use crate::common::run_modkit;
use anyhow::Context;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::Command;

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
fn test_validate_bed_errors_preserve_output_and_report_line() {
    let temp_dir = tempfile::tempdir().unwrap();
    let cases = [
        ("field-count", "not-a-bed-row", "Invalid number of fields"),
        (
            "coordinate",
            "chr1\t5\t5\tm\t.\t+",
            "BED end must be greater than start",
        ),
        ("strand", "chr1\t5\t6\tm\t.\t+junk", "expected `+` or `-`"),
        ("mod-code", "chr1\t5\t6\t.\t.\t+", "failed to parse mod code"),
        (
            "conflicting-label",
            "chr1\t4\t6\th\t.\t+",
            "conflicting ground truth labels",
        ),
    ];

    for (case_name, invalid_line, expected_error) in cases {
        let bed_path = temp_dir.path().join(format!("{case_name}.bed"));
        let output_path = temp_dir.path().join(format!("{case_name}.tsv"));
        std::fs::write(
            &bed_path,
            format!("chr1\t4\t5\tm\t.\t+\n{invalid_line}\n"),
        )
        .unwrap();
        std::fs::write(&output_path, "sentinel\n").unwrap();

        let output = Command::new(Path::new(env!("CARGO_BIN_EXE_modkit")))
            .arg("validate")
            .arg("--bam-and-bed")
            .arg("../tests/resources/input_5mC.bam")
            .arg(&bed_path)
            .arg("--canonical-base")
            .arg("C")
            .arg("--out-filepath")
            .arg(&output_path)
            .arg("--suppress-progress")
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);

        assert!(!output.status.success(), "{case_name}: {stderr}");
        assert!(
            stderr.contains(&bed_path.display().to_string()),
            "{case_name}: {stderr}"
        );
        assert!(stderr.contains("line 2"), "{case_name}: {stderr}");
        assert!(stderr.contains(expected_error), "{case_name}: {stderr}");
        assert_eq!(
            std::fs::read(&output_path).unwrap(),
            b"sentinel\n",
            "{case_name}"
        );
    }
}
