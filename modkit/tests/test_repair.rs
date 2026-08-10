use std::path::Path;

use mod_kit::mod_bam::MN_TAG;
use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::Read;

use crate::common::{run_modkit, run_simple_summary};

mod common;

#[test]
fn test_repair_help() {
    run_modkit(&["repair", "--help"]).unwrap();
}

#[test]
fn test_repair_regression() {
    let temp_dir = tempfile::tempdir().unwrap();
    let out_bam = temp_dir.path().join("test_repair_regression.bam");
    let donor_fp = "../tests/resources/donor_read_sort.bam";
    let acceptor_fp = "../tests/resources/trimmed_read_sort.mapped.bam";
    run_modkit(&[
        "repair",
        "--donor",
        donor_fp,
        "--acceptor",
        acceptor_fp,
        "-o",
        out_bam.to_str().unwrap(),
    ])
    .unwrap();

    // dbg!(out_bam.to_str().unwrap());
    let mut test_bam = bam::Reader::from_path(out_bam).unwrap();
    // this BAM was hand-checked
    let mut ref_bam = bam::Reader::from_path(
        "../tests/resources/trimmed_read_sort_mods.mapped.bam",
    )
    .unwrap();
    let mut test_records = test_bam
        .records()
        .map(|r| r.unwrap())
        .map(|mut record| {
            // todo consider removing this later, added MN tag by default but
            // the old  test data doesn't have it.
            record.remove_aux(MN_TAG.as_bytes()).expect("should remove MN tag");
            record
        })
        .collect::<Vec<_>>();

    let mut expected_records =
        ref_bam.records().map(|r| r.unwrap()).collect::<Vec<_>>();

    // The hand-checked legacy reference predates ordered repair output. Sort
    // both complete vectors for the content oracle; the synthetic regression
    // below independently checks exact acceptor-relative output order.
    test_records.sort_by(|left, right| left.qname().cmp(right.qname()));
    expected_records.sort_by(|left, right| left.qname().cmp(right.qname()));
    assert_eq!(
        test_records.len(),
        expected_records.len(),
        "repair output cardinality changed"
    );
    for (index, (actual, expected)) in
        test_records.iter().zip(expected_records.iter()).enumerate()
    {
        assert_eq!(actual.qname(), expected.qname(), "record {index} missing");
        assert_eq!(actual, expected, "record {index} not the same");
    }
}

fn query_name_header(sub_sort: Option<&str>) -> bam::Header {
    let mut header = bam::Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6").push_tag(b"SO", "queryname");
    if let Some(sub_sort) = sub_sort {
        hd.push_tag(b"SS", sub_sort);
    }
    header.push_record(&hd);
    header
}

fn write_repair_input(
    path: &Path,
    records: &[(&str, u32)],
    donor_probabilities: Option<&[u8]>,
    sub_sort: Option<&str>,
) {
    if let Some(probabilities) = donor_probabilities {
        assert_eq!(probabilities.len(), records.len());
    }
    let header = query_name_header(sub_sort);
    let mut writer =
        bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
    for (index, (name, marker)) in records.iter().enumerate() {
        let sequence =
            if donor_probabilities.is_some() { "TACGT" } else { "ACG" };
        let mut record = bam::Record::new();
        record.set(
            name.as_bytes(),
            None,
            sequence.as_bytes(),
            &vec![255; sequence.len()],
        );
        record.push_aux(b"XI", Aux::U32(*marker)).unwrap();
        if let Some(probabilities) = donor_probabilities {
            record.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
            let probability = [probabilities[index]];
            let probability_array: AuxArray<u8> = (&probability[..]).into();
            record.push_aux(b"ML", Aux::ArrayU8(probability_array)).unwrap();
        }
        writer.write(&record).unwrap();
    }
}

fn repaired_identities(path: &Path) -> Vec<(String, u32, u8)> {
    let mut reader = bam::Reader::from_path(path).unwrap();
    reader
        .records()
        .map(|record| {
            let record = record.unwrap();
            let marker = match record.aux(b"XI").unwrap() {
                Aux::U32(marker) => marker,
                other => panic!("unexpected marker tag: {other:?}"),
            };
            let probability = match record.aux(b"ML").unwrap() {
                Aux::ArrayU8(probabilities) => {
                    let probabilities =
                        probabilities.iter().collect::<Vec<_>>();
                    assert_eq!(probabilities.len(), 1);
                    probabilities[0]
                }
                other => panic!("unexpected ML tag: {other:?}"),
            };
            (
                String::from_utf8_lossy(record.qname()).into_owned(),
                marker,
                probability,
            )
        })
        .collect()
}

fn repair_error(donor: &Path, acceptor: &Path, output: &Path) -> String {
    let executable = Path::new(env!("CARGO_BIN_EXE_modkit"));
    let result = std::process::Command::new(executable)
        .args([
            "repair",
            "--donor",
            donor.to_str().unwrap(),
            "--acceptor",
            acceptor.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--threads",
            "2",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success(), "repair unexpectedly succeeded");
    String::from_utf8_lossy(&result.stderr).into_owned()
}

#[test]
fn test_repair_preserves_acceptor_order_and_exact_cardinality() {
    let temp_dir = tempfile::tempdir().unwrap();
    let donor = temp_dir.path().join("donor.bam");
    let acceptor = temp_dir.path().join("acceptor.bam");
    let output = temp_dir.path().join("repaired.bam");
    write_repair_input(
        &donor,
        &[
            ("read01", 90),
            ("read1", 91),
            ("read001", 92),
            ("read3", 93),
            ("read10", 100),
        ],
        Some(&[22, 11, 33, 44, 55]),
        Some("queryname:natural"),
    );
    write_repair_input(
        &acceptor,
        &[
            ("read1", 0),
            ("read1", 1),
            ("read01", 2),
            ("read001", 7),
            ("read0001", 8),
            ("read2", 3),
            ("read3", 4),
            ("read3", 5),
            ("read10", 6),
        ],
        None,
        Some("queryname:natural"),
    );

    run_modkit(&[
        "repair",
        "--donor",
        donor.to_str().unwrap(),
        "--acceptor",
        acceptor.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--threads",
        "2",
    ])
    .unwrap();

    let reader = bam::Reader::from_path(&output).unwrap();
    let header = String::from_utf8_lossy(reader.header().as_bytes());
    assert!(header.contains("SO:queryname"));
    assert!(header.contains("SS:queryname:natural"));
    drop(reader);
    let actual = repaired_identities(&output);
    assert_eq!(
        actual,
        vec![
            ("read1".to_string(), 0, 11),
            ("read1".to_string(), 1, 11),
            ("read01".to_string(), 2, 22),
            ("read001".to_string(), 7, 33),
            ("read3".to_string(), 4, 44),
            ("read3".to_string(), 5, 44),
            ("read10".to_string(), 6, 55),
        ]
    );
}

#[test]
fn test_repair_rejects_mixed_natural_headers_in_both_directions() {
    let temp_dir = tempfile::tempdir().unwrap();
    for (index, donor_uses_explicit_natural) in
        [true, false].into_iter().enumerate()
    {
        let donor = temp_dir.path().join(format!("donor-{index}.bam"));
        let acceptor = temp_dir.path().join(format!("acceptor-{index}.bam"));
        let output = temp_dir.path().join(format!("repaired-{index}.bam"));
        if donor_uses_explicit_natural {
            write_repair_input(
                &donor,
                &[("read1a1", 1), ("read01a2", 2)],
                Some(&[11, 22]),
                Some("queryname:natural"),
            );
            write_repair_input(
                &acceptor,
                &[("read01a2", 2), ("read1a1", 1)],
                None,
                None,
            );
        } else {
            write_repair_input(
                &donor,
                &[("read01a2", 2), ("read1a1", 1)],
                Some(&[22, 11]),
                None,
            );
            write_repair_input(
                &acceptor,
                &[("read1a1", 1), ("read01a2", 2)],
                None,
                Some("queryname:natural"),
            );
        }

        let expected_diagnostic = if donor_uses_explicit_natural {
            "donor and acceptor BAMs use incompatible query-name sorting: \
             donor is natural, acceptor is natural (SO-only legacy fallback)"
        } else {
            "donor and acceptor BAMs use incompatible query-name sorting: \
             donor is natural (SO-only legacy fallback), acceptor is natural"
        };
        let error = repair_error(&donor, &acceptor, &output);
        assert!(
            error.contains(expected_diagnostic),
            "unexpected mixed-sort error: {error}"
        );
    }
}

#[test]
fn test_repair_rejects_ambiguous_or_incompatibly_sorted_inputs() {
    let temp_dir = tempfile::tempdir().unwrap();
    let donor = temp_dir.path().join("donor.bam");
    let acceptor = temp_dir.path().join("acceptor.bam");
    let output = temp_dir.path().join("repaired.bam");
    write_repair_input(
        &donor,
        &[("read1", 0), ("read1", 1)],
        Some(&[100, 101]),
        Some("queryname:natural"),
    );
    write_repair_input(
        &acceptor,
        &[("read1", 0)],
        None,
        Some("queryname:natural"),
    );
    let error = repair_error(&donor, &acceptor, &output);
    assert!(
        error.contains("multiple primary records for query name read1"),
        "unexpected duplicate-donor error: {error}"
    );

    write_repair_input(
        &donor,
        &[("read2", 0), ("read1", 1)],
        Some(&[100, 101]),
        Some("queryname:natural"),
    );
    let error = repair_error(&donor, &acceptor, &output);
    assert!(
        error.contains("donor BAM is not natural query-name sorted"),
        "unexpected donor-inversion error: {error}"
    );

    write_repair_input(
        &donor,
        &[("read1", 0), ("read2", 1)],
        Some(&[100, 101]),
        Some("queryname:natural"),
    );
    write_repair_input(
        &acceptor,
        &[("read2", 0), ("read1", 1)],
        None,
        Some("queryname:natural"),
    );
    let error = repair_error(&donor, &acceptor, &output);
    assert!(
        error.contains("acceptor BAM is not natural query-name sorted"),
        "unexpected acceptor-inversion error: {error}"
    );

    write_repair_input(
        &donor,
        &[("read1", 0)],
        Some(&[100]),
        Some("queryname:lexicographical"),
    );
    write_repair_input(
        &acceptor,
        &[("read1", 0)],
        None,
        Some("queryname:natural"),
    );
    let error = repair_error(&donor, &acceptor, &output);
    assert!(
        error.contains("incompatible query-name sorting"),
        "unexpected incompatible-sort error: {error}"
    );
}

#[test]
fn test_repair_mn_tag() {
    let temp_dir = tempfile::tempdir().unwrap();
    let out_bam = temp_dir.path().join("test_repair_mn_tag.bam");
    let donor_fp = "../tests/resources/donor_read_sort_mn_tag.bam";
    let acceptor_fp = "../tests/resources/trimmed_read_sort_mn_tag.mapped.bam";
    run_modkit(&[
        "repair",
        "--donor",
        donor_fp,
        "--acceptor",
        acceptor_fp,
        "-o",
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    let summary = run_simple_summary(out_bam.to_str().unwrap(), 25)
        .expect("should run summary");
    assert_eq!(summary.total_reads_used, 10);
}
