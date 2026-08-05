use crate::common::run_modkit;
use anyhow::Context;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{Format, Header, Record, Writer};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use tempfile::tempdir;

mod common;

fn run_validate(
    bam_bed_pairs: &[(&Path, &Path)],
    output_file: &Path,
) -> (String, HashMap<String, String>) {
    let mut args = vec!["validate".to_string()];
    for (bam, bed) in bam_bed_pairs {
        args.push("--bam-and-bed".to_string());
        args.push(bam.to_str().unwrap().to_string());
        args.push(bed.to_str().unwrap().to_string());
    }
    args.extend([
        "--out-filepath".to_string(),
        output_file.to_str().unwrap().to_string(),
        "--filter-threshold".to_string(),
        "0".to_string(),
        "--suppress-progress".to_string(),
        "--threads".to_string(),
        "1".to_string(),
    ]);
    let args = args.iter().map(String::as_str).collect::<Vec<_>>();
    run_modkit(&args).unwrap();

    let output = fs::read_to_string(output_file).unwrap();
    let fields = output
        .lines()
        .map(|line| {
            let (key, value) = line.split_once(": ").unwrap();
            (key.to_string(), value.to_string())
        })
        .collect();
    (output, fields)
}

#[test]
fn test_validate_help() {
    run_modkit(&["validate", "--help"])
        .context("validate --help failed")
        .unwrap();
}

#[test]
fn test_validate_expected() {
    let output_file = std::env::temp_dir().join("test_validate_output.tsv");
    run_modkit(&[
        "validate",
        "--bam-and-bed",
        "../tests/resources/input_5mC.bam",
        "../tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed",
        "--bam-and-bed",
        "../tests/resources/input_C.bam",
        "../tests/resources/CGI_ladder_3.6kb_ref_CG_C.bed",
        "--out-filepath",
        output_file.to_str().unwrap(),
    ])
    .context("should run validate")
    .unwrap();

    let reader = BufReader::new(File::open(output_file).unwrap());
    for line in reader.lines() {
        let line_content = line.unwrap();
        if line_content.starts_with("raw_accuracy") {
            let accuracy: f32 = line_content
                .split_whitespace()
                .skip(1)
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or_default();
            assert_eq!(accuracy, 85.305214);
        }
        if line_content.starts_with("filtered_accuracy") {
            let accuracy: f32 = line_content
                .split_whitespace()
                .skip(1)
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or_default();
            assert_eq!(accuracy, 89.00287);
        }
    }
}

#[test]
fn test_validate_reopens_bam_for_each_truth_bed() {
    let temp_dir = tempdir().unwrap();
    let bam = Path::new("../tests/resources/input_5mC.bam");
    let truth = Path::new("../tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed");
    let truth_lines = fs::read_to_string(truth)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect::<Vec<_>>();
    assert_eq!(truth_lines.len(), 48);

    let mut split_a_lines = Vec::new();
    let mut split_b_lines = Vec::new();
    for (idx, line) in truth_lines.iter().enumerate() {
        if idx % 2 == 0 {
            split_a_lines.push(line.clone());
        } else {
            split_b_lines.push(line.clone());
        }
    }
    let truth_set = truth_lines.iter().collect::<HashSet<_>>();
    let split_a_set = split_a_lines.iter().collect::<HashSet<_>>();
    let split_b_set = split_b_lines.iter().collect::<HashSet<_>>();
    assert_eq!(truth_set.len(), truth_lines.len());
    assert!(split_a_set.is_disjoint(&split_b_set));
    assert_eq!(
        split_a_set.union(&split_b_set).copied().collect::<HashSet<_>>(),
        truth_set
    );

    let split_a = temp_dir.path().join("truth-a.bed");
    let split_b = temp_dir.path().join("truth-b.bed");
    fs::write(&split_a, format!("{}\n", split_a_lines.join("\n"))).unwrap();
    fs::write(&split_b, format!("{}\n", split_b_lines.join("\n"))).unwrap();

    let (union_output, union) = run_validate(
        &[(bam, truth)],
        &temp_dir.path().join("union-output.tsv"),
    );
    assert_eq!(
        union["full_contingency_table"],
        "[[\"ground_truth_label\",\"C\",\"h\",\"m\",\"No Call\",\"A\",\"G\",\"T\",\"Deletion\"],[\"m\",16625,20432,120926,9069,4826,8337,859,4608]]"
    );
    assert_eq!(
        union["raw_contingency_table"],
        "[[\"ground_truth_label\",\"C\",\"h\",\"m\"],[\"m\",16625,20432,120926]]"
    );

    let (split_ab_output, split_ab) = run_validate(
        &[(bam, split_a.as_path()), (bam, split_b.as_path())],
        &temp_dir.path().join("split-ab-output.tsv"),
    );
    assert_eq!(
        split_ab["full_contingency_table"],
        union["full_contingency_table"]
    );
    assert_eq!(split_ab, union);
    assert_eq!(split_ab_output, union_output);

    let (split_ba_output, split_ba) = run_validate(
        &[(bam, split_b.as_path()), (bam, split_a.as_path())],
        &temp_dir.path().join("split-ba-output.tsv"),
    );
    assert_eq!(
        split_ba["full_contingency_table"],
        union["full_contingency_table"]
    );
    assert_eq!(split_ba, union);
    assert_eq!(split_ba_output, union_output);
}

fn synthetic_record(
    name: &str,
    sequence: &str,
    cigar: Vec<Cigar>,
    mm: &str,
    ml: &[u8],
    nm: u32,
) -> Record {
    let mut record = Record::new();
    let cigar = CigarString(cigar);
    record.set(
        name.as_bytes(),
        Some(&cigar),
        sequence.as_bytes(),
        &vec![30; sequence.len()],
    );
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    record.push_aux(b"MM", Aux::String(mm)).unwrap();
    record.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
    record.push_aux(b"MN", Aux::U32(sequence.len() as u32)).unwrap();
    record.push_aux(b"NM", Aux::U32(nm)).unwrap();
    record
}

fn synthetic_called(name: &str, mod_strand: char) -> Record {
    let (sequence, fundamental_base) = match mod_strand {
        '+' => ("CC", 'C'),
        '-' => ("GG", 'G'),
        _ => panic!("invalid modification strand"),
    };
    synthetic_record(
        name,
        sequence,
        vec![Cigar::Match(2)],
        &format!("{fundamental_base}{mod_strand}m?,1;"),
        &[255],
        0,
    )
}

fn synthetic_uncalled(name: &str, sequence: &str, mod_strand: char) -> Record {
    let (expected_sequence, fundamental_base) = match mod_strand {
        '+' => ("CC", 'C'),
        '-' => ("GG", 'G'),
        _ => panic!("invalid modification strand"),
    };
    synthetic_record(
        name,
        sequence,
        vec![Cigar::Match(2)],
        &format!("{fundamental_base}{mod_strand}m?,0;"),
        &[255],
        u32::from(sequence != expected_sequence),
    )
}

fn write_synthetic_bam(path: &Path, records: Vec<Record>) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 2);
    header.push_record(&sq);
    let mut writer = Writer::from_path(path, &header, Format::Bam).unwrap();
    for record in records {
        writer.write(&record).unwrap();
    }
}

fn assert_validate_fixture(
    root: &Path,
    name: &str,
    records: Vec<Record>,
    bed: &str,
    expected_full_table: &str,
) {
    let bam = root.join(format!("{name}.bam"));
    let bed_path = root.join(format!("{name}.bed"));
    let output = root.join(format!("{name}.tsv"));
    let check_tags_dir = root.join(format!("{name}-check-tags"));
    write_synthetic_bam(&bam, records);
    let mut bed_file = File::create(&bed_path).unwrap();
    bed_file.write_all(bed.as_bytes()).unwrap();

    run_modkit(&[
        "modbam",
        "check-tags",
        bam.to_str().unwrap(),
        "--ignore-index",
        "--suppress-progress",
        "--out-dir",
        check_tags_dir.to_str().unwrap(),
    ])
    .with_context(|| {
        format!("{name}: synthetic records should pass check-tags")
    })
    .unwrap();

    run_modkit(&[
        "validate",
        "--bam-and-bed",
        bam.to_str().unwrap(),
        bed_path.to_str().unwrap(),
        "--canonical-base",
        "C",
        "--filter-threshold",
        "0",
        "--threads",
        "1",
        "--suppress-progress",
        "--out-filepath",
        output.to_str().unwrap(),
    ])
    .with_context(|| format!("{name}: validate should succeed"))
    .unwrap();

    let first_line = BufReader::new(File::open(output).unwrap())
        .lines()
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(
        first_line,
        format!("full_contingency_table: {expected_full_table}"),
        "{name}"
    );
}

#[test]
fn test_validate_counts_observations_without_truth_overlapping_seed_calls() {
    let temp_dir = tempfile::tempdir().unwrap();
    let root = temp_dir.path();

    let gt_m_records =
        (0..3)
            .map(|i| synthetic_called(&format!("called-{i}"), '+'))
            .chain((0..6).map(|i| {
                synthetic_uncalled(&format!("uncalled-{i}"), "CC", '+')
            }))
            .collect();
    assert_validate_fixture(
        root,
        "gt-m",
        gt_m_records,
        "chr1\t1\t2\tm\t0\t+\n",
        "[[\"ground_truth_label\",\"m\",\"No Call\"],[\"m\",3,6]]",
    );

    let mirror_records =
        (0..2)
            .map(|i| synthetic_called(&format!("called-{i}"), '+'))
            .chain((0..6).map(|i| {
                synthetic_uncalled(&format!("uncalled-{i}"), "CC", '+')
            }))
            .chain((0..2).map(|i| {
                synthetic_uncalled(&format!("mismatch-{i}"), "CA", '+')
            }))
            .chain((0..2).map(|i| {
                synthetic_record(
                    &format!("deletion-{i}"),
                    "C",
                    vec![Cigar::Match(1), Cigar::Del(1)],
                    "C+m?,0;",
                    &[255],
                    1,
                )
            }))
            .collect();
    assert_validate_fixture(
        root,
        "mirror",
        mirror_records,
        "chr1\t1\t2\tm\t0\t+\n",
        "[[\"ground_truth_label\",\"m\",\"No Call\",\"A\",\"Deletion\"],[\"m\",2,6,2,2]]",
    );

    let both_strands_records =
        std::iter::once(synthetic_called("positive-called", '+'))
            .chain((0..5).map(|i| {
                synthetic_uncalled(&format!("positive-uncalled-{i}"), "CC", '+')
            }))
            .chain(std::iter::once(synthetic_called("negative-called", '-')))
            .chain((0..5).map(|i| {
                synthetic_uncalled(&format!("negative-uncalled-{i}"), "GG", '-')
            }))
            .collect();
    assert_validate_fixture(
        root,
        "both-strands",
        both_strands_records,
        "chr1\t1\t2\tm\t0\t+\nchr1\t1\t2\tm\t0\t-\n",
        "[[\"ground_truth_label\",\"m\",\"No Call\"],[\"m\",2,10]]",
    );

    let empty_descriptor_records = vec![
        synthetic_called("called-positive", '+'),
        synthetic_record(
            "empty-positive",
            "CC",
            vec![Cigar::Match(2)],
            "C+m?;",
            &[],
            0,
        ),
        synthetic_record(
            "empty-negative",
            "GG",
            vec![Cigar::Match(2)],
            "G-m?;",
            &[],
            0,
        ),
    ];
    assert_validate_fixture(
        root,
        "empty-descriptors",
        empty_descriptor_records,
        "chr1\t1\t2\tm\t0\t+\nchr1\t1\t2\tm\t0\t-\n",
        "[[\"ground_truth_label\",\"m\",\"No Call\"],[\"m\",1,2]]",
    );

    let implicit_records = (0..3)
        .map(|i| {
            synthetic_record(
                &format!("implicit-{i}"),
                "CC",
                vec![Cigar::Match(2)],
                "C+m.,0;",
                &[255],
                0,
            )
        })
        .collect();
    assert_validate_fixture(
        root,
        "implicit",
        implicit_records,
        "chr1\t1\t2\t-\t0\t+\n",
        "[[\"ground_truth_label\",\"C\"],[\"C\",3]]",
    );
}
