use std::{fs, path::Path, process::Command};

use crate::common::run_modkit;

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
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        "-o",
        td.to_str().unwrap(),
        "--min-coverage",
        "1",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--regions",
        "../tests/resources/entropy_test_regions.bed",
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
    // "../tests/resources/expected_entropy_regions.bed");
    // check_against_expected_text_file(windows.to_str().unwrap(),
    // "../tests/resources/expected_entropy_windows.bed");
}

#[test]
fn no_filtering_bypasses_sampling_and_matches_zero_threshold() {
    let temp_dir = tempfile::tempdir().unwrap();
    let no_filter_output = temp_dir.path().join("no_filter.bed");
    let zero_threshold_output = temp_dir.path().join("zero_threshold.bed");
    let executable = Path::new(env!("CARGO_BIN_EXE_modkit"));

    let run = |output: &Path, threshold_args: &[&str]| {
        Command::new(executable)
            .args([
                "entropy",
                "--in-bam",
                "../tests/resources/bc_anchored_10_reads.sorted.bam",
                "--out-bed",
                output.to_str().unwrap(),
                "--ref",
                "../tests/resources/CGI_ladder_3.6kb_ref.fa",
                "--cpg",
                "--min-coverage",
                "1",
                "--threads",
                "1",
                "--io-threads",
                "1",
                "--num-reads",
                "0",
                "--suppress-progress",
            ])
            .args(threshold_args)
            .output()
            .unwrap()
    };

    let no_filter = run(&no_filter_output, &["--no-filtering"]);
    assert!(
        no_filter.status.success(),
        "no-filtering failed:\n{}",
        String::from_utf8_lossy(&no_filter.stderr)
    );
    let no_filter_stderr = String::from_utf8_lossy(&no_filter.stderr);
    assert!(no_filter_stderr.contains("not performing filtering"));
    assert!(!no_filter_stderr.contains("calculated thresholds:"));

    let zero_threshold =
        run(&zero_threshold_output, &["--filter-threshold", "0"]);
    assert!(
        zero_threshold.status.success(),
        "zero threshold failed:\n{}",
        String::from_utf8_lossy(&zero_threshold.stderr)
    );

    let no_filter_text = fs::read_to_string(no_filter_output).unwrap();
    let zero_threshold_text =
        fs::read_to_string(zero_threshold_output).unwrap();
    let no_filter_rows = no_filter_text.lines().collect::<Vec<_>>();
    let zero_threshold_rows = zero_threshold_text.lines().collect::<Vec<_>>();
    assert_eq!(no_filter_rows.len(), zero_threshold_rows.len());
    for (no_filter_row, zero_threshold_row) in
        no_filter_rows.into_iter().zip(zero_threshold_rows)
    {
        let no_filter_fields = no_filter_row.split('\t').collect::<Vec<_>>();
        let zero_threshold_fields =
            zero_threshold_row.split('\t').collect::<Vec<_>>();
        assert_eq!(&no_filter_fields[..3], &zero_threshold_fields[..3]);
        assert_eq!(&no_filter_fields[4..], &zero_threshold_fields[4..]);
        let no_filter_entropy = no_filter_fields[3].parse::<f64>().unwrap();
        let zero_threshold_entropy =
            zero_threshold_fields[3].parse::<f64>().unwrap();
        assert!(
            (no_filter_entropy - zero_threshold_entropy).abs() <= 1e-6,
            "entropy differed: {no_filter_entropy} versus \
             {zero_threshold_entropy}"
        );
    }
}

#[test]
fn no_filtering_conflicts_with_per_mod_thresholds_before_output() {
    let temp_dir = tempfile::tempdir().unwrap();
    let absent_output = temp_dir.path().join("absent.bed");
    let sentinel_output = temp_dir.path().join("sentinel.bed");
    fs::write(&sentinel_output, "keep-sentinel\n").unwrap();
    let executable = Path::new(env!("CARGO_BIN_EXE_modkit"));

    for (threshold_option, threshold_value, output) in [
        ("--mod-threshold", "m:0.5", &absent_output),
        ("--mod-thresholds", "not-valid", &sentinel_output),
    ] {
        let result = Command::new(executable)
            .args([
                "entropy",
                "--in-bam",
                "../tests/resources/bc_anchored_10_reads.sorted.bam",
                "--out-bed",
                output.to_str().unwrap(),
                "--ref",
                "../tests/resources/CGI_ladder_3.6kb_ref.fa",
                "--cpg",
                "--no-filtering",
                threshold_option,
                threshold_value,
                "--suppress-progress",
            ])
            .output()
            .unwrap();

        assert_eq!(result.status.code(), Some(2));
        assert!(
            String::from_utf8_lossy(&result.stderr)
                .contains("cannot be used with"),
            "unexpected stderr:\n{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    assert!(!absent_output.exists());
    assert_eq!(fs::read_to_string(sentinel_output).unwrap(), "keep-sentinel\n");
}
