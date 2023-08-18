use std::collections::HashMap;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use crate::common::run_modkit;

mod common;

#[test]
fn test_repair_help() {
    run_modkit(&["repair", "--help"]).unwrap();
}

#[test]
fn test_repair_regression() {
    let out_bam = std::env::temp_dir().join("test_repair_regression.bam");
    let donor_fp = "tests/resources/donor_read_sort.bam";
    let acceptor_fp = "tests/resources/trimmed_read_sort.mapped.bam";
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
        "tests/resources/trimmed_read_sort_mods.mapped.bam",
    )
    .unwrap();
    let test_records = test_bam
        .records()
        .map(|r| r.unwrap())
        .map(|record| {
            let qname = record
                .qname()
                .iter()
                .map(|&b| b as char)
                .collect::<String>();
            (qname, record)
        })
        .collect::<HashMap<_, _>>();

    let expected_records = ref_bam
        .records()
        .map(|r| r.unwrap())
        .map(|record| {
            let qname = record
                .qname()
                .iter()
                .map(|&b| b as char)
                .collect::<String>();
            (qname, record)
        })
        .collect::<HashMap<_, _>>();

    for (q, r) in test_records.iter() {
        assert_eq!(
            expected_records.get(q).unwrap(),
            r,
            "record {q} not the same"
        );
    }
}
