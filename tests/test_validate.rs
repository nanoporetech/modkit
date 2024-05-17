use crate::common::run_modkit;
use anyhow::Context;
use std::fs::File;
use std::io::{BufRead, BufReader};

mod common;

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
        "tests/resources/input_5mC.bam",
        "tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed",
        "--bam-and-bed",
        "tests/resources/input_C.bam",
        "tests/resources/CGI_ladder_3.6kb_ref_CG_C.bed",
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
