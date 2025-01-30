use crate::common::{check_against_expected_text_file, run_modkit};

mod common;

#[test]
fn test_modbam_helps() {
    let _ = run_modkit(&["modbam", "--help"])
        .expect("failed to run modkit modbam help");
    let _ = run_modkit(&["modbam", "check-tags", "--help"])
        .expect("failed to run modkit modbam check-tags help");
    let _ = run_modkit(&["mb", "--help"])
        .expect("failed to run modkit modbam help");
    let _ = run_modkit(&["mb", "check-tags", "--help"])
        .expect("failed to run modkit modbam check-tags help");
}

#[test]
fn test_modbam_check_tags_expected_valid_reads_output() {
    let tmp_dir = std::env::temp_dir().join("test_expected_valid_reads_output");
    let _ = run_modkit(&[
        "modbam",
        "check-tags",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        "--interval-size",
        "20",
        "--force",
        "--out-dir",
        tmp_dir.to_str().unwrap(),
    ])
    .unwrap();
    check_against_expected_text_file(
        tmp_dir.join("modified_bases.tsv").to_str().unwrap(),
        "tests/resources/modified_bases.tsv",
    );
    check_against_expected_text_file(
        tmp_dir.join("valid_mm_headers.tsv").to_str().unwrap(),
        "tests/resources/valid_mm_headers.tsv",
    );
}
