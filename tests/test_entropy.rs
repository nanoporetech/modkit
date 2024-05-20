use crate::common::{check_against_expected_text_file, run_modkit};

mod common;

#[test]
fn test_entropy_help() {
    run_modkit(&["entropy", "--help"]).expect("entropy help");
}

#[test]
fn test_entropy_regression() {
    let td = std::env::temp_dir().join("test_entropy_regression");
    std::fs::create_dir_all(&td).expect("should make temp dir");
    run_modkit(&[
        "entropy",
        "-s",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        "-o",
        td.to_str().unwrap(),
        "--min-coverage",
        "1",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--regions",
        "tests/resources/entropy_test_regions.bed",
        "--prefix",
        "prefix_test",
        "--cpg",
        "--force",
    ])
    .expect("should run entropy on regions");
    let regions = td.join("prefix_test_regions.bed");
    let windows = td.join("prefix_test_windows.bedgraph");
    assert!(regions.exists());
    assert!(windows.exists());
    // todo too much wiggle in the calculation, make an assert_approx on the
    // scores check_against_expected_text_file(regions.to_str().unwrap(),
    // "tests/resources/expected_entropy_regions.bed");
    // check_against_expected_text_file(windows.to_str().unwrap(),
    // "tests/resources/expected_entropy_windows.bed");
}
