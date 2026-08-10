use std::{fs, path::Path, process::Command};

use crate::common::run_modkit;

mod common;

fn entropy_command() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "entropy",
        "--in-bam",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--cpg",
        "--min-coverage",
        "1",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--suppress-progress",
    ]);
    command
}

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

#[test]
fn threshold_failure_does_not_create_or_truncate_entropy_output() {
    let temp_dir = tempfile::tempdir().unwrap();
    let absent_output = temp_dir.path().join("absent.bed");
    let sentinel_output = temp_dir.path().join("sentinel.bed");
    let sentinel = b"existing scientific output\n";
    fs::write(&sentinel_output, sentinel).unwrap();

    for output_path in [&absent_output, &sentinel_output] {
        let output = entropy_command()
            .args([
                "--out-bed",
                output_path.to_str().unwrap(),
                "--num-reads",
                "0",
                "--header",
                "--force",
            ])
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "empty automatic-threshold setup unexpectedly succeeded"
        );
    }

    assert!(!absent_output.exists());
    assert_eq!(fs::read(sentinel_output).unwrap(), sentinel);
}

#[test]
fn entropy_file_output_requires_force() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_path = temp_dir.path().join("entropy.bed");
    let sentinel = b"existing scientific output\n";
    fs::write(&output_path, sentinel).unwrap();

    let rejected = entropy_command()
        .args([
            "--out-bed",
            output_path.to_str().unwrap(),
            "--filter-threshold",
            "0",
        ])
        .output()
        .unwrap();
    assert!(!rejected.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), sentinel);

    let overwritten = entropy_command()
        .args([
            "--out-bed",
            output_path.to_str().unwrap(),
            "--filter-threshold",
            "0",
            "--force",
        ])
        .output()
        .unwrap();
    assert!(
        overwritten.status.success(),
        "forced output failed:\n{}",
        String::from_utf8_lossy(&overwritten.stderr)
    );
    assert_ne!(fs::read(output_path).unwrap(), sentinel);
}

#[test]
fn entropy_regions_preflight_both_prefixed_targets() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_dir = temp_dir.path().join("entropy");
    fs::create_dir(&output_dir).unwrap();
    let unrelated = output_dir.join("keep.txt");
    let unrelated_contents = b"unrelated contents\n";
    fs::write(&unrelated, unrelated_contents).unwrap();

    let regions_output = output_dir.join("sample_regions.bed");
    let windows_output = output_dir.join("sample_windows.bedgraph");
    let sentinel = b"existing windows output\n";
    fs::write(&windows_output, sentinel).unwrap();

    let run_regions = |prefix: &str, force: bool| {
        let mut command = entropy_command();
        command.args([
            "--out-bed",
            output_dir.to_str().unwrap(),
            "--regions",
            "../tests/resources/entropy_test_regions.bed",
            "--prefix",
            prefix,
            "--filter-threshold",
            "0",
        ]);
        if force {
            command.arg("--force");
        }
        command.output().unwrap()
    };

    let rejected = run_regions("sample", false);
    assert!(!rejected.status.success());
    assert!(!regions_output.exists());
    assert_eq!(fs::read(&windows_output).unwrap(), sentinel);
    assert_eq!(fs::read(&unrelated).unwrap(), unrelated_contents);

    let fresh = run_regions("fresh", false);
    assert!(
        fresh.status.success(),
        "fresh prefixed output failed:\n{}",
        String::from_utf8_lossy(&fresh.stderr)
    );
    assert!(output_dir.join("fresh_regions.bed").exists());
    assert!(output_dir.join("fresh_windows.bedgraph").exists());
    assert_eq!(fs::read(&unrelated).unwrap(), unrelated_contents);

    let forced = run_regions("sample", true);
    assert!(
        forced.status.success(),
        "forced prefixed output failed:\n{}",
        String::from_utf8_lossy(&forced.stderr)
    );
    assert!(regions_output.exists());
    assert_ne!(fs::read(windows_output).unwrap(), sentinel);
    assert_eq!(fs::read(unrelated).unwrap(), unrelated_contents);
}

#[test]
fn entropy_regions_rejects_escaping_prefixes_before_forced_mutation() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_dir = temp_dir.path().join("entropy");
    fs::create_dir(&output_dir).unwrap();
    fs::create_dir(output_dir.join("a")).unwrap();
    let unrelated = output_dir.join("keep.txt");
    let unrelated_contents = b"unrelated contents\n";
    fs::write(&unrelated, unrelated_contents).unwrap();

    let absolute_prefix = temp_dir.path().join("absolute-x");
    let absolute_prefix_string = absolute_prefix.to_str().unwrap().to_string();
    let sentinel = b"outside scientific output\n";
    let invalid_cases = [
        (
            "../x".to_string(),
            [
                temp_dir.path().join("x_regions.bed"),
                temp_dir.path().join("x_windows.bedgraph"),
            ],
        ),
        (
            "a/b".to_string(),
            [
                output_dir.join("a/b_regions.bed"),
                output_dir.join("a/b_windows.bedgraph"),
            ],
        ),
        (
            absolute_prefix_string.clone(),
            [
                format!("{absolute_prefix_string}_regions.bed").into(),
                format!("{absolute_prefix_string}_windows.bedgraph").into(),
            ],
        ),
    ];

    for (prefix, targets) in invalid_cases {
        for target in &targets {
            fs::write(target, sentinel).unwrap();
        }
        let rejected = entropy_command()
            .args([
                "--out-bed",
                output_dir.to_str().unwrap(),
                "--regions",
                "../tests/resources/entropy_test_regions.bed",
                "--prefix",
                &prefix,
                "--filter-threshold",
                "0",
                "--force",
            ])
            .output()
            .unwrap();
        assert!(!rejected.status.success(), "prefix {prefix:?} was accepted");
        assert!(
            String::from_utf8_lossy(&rejected.stderr)
                .contains("exactly one non-empty filename component"),
            "unexpected stderr for {prefix:?}:\n{}",
            String::from_utf8_lossy(&rejected.stderr)
        );
        for target in targets {
            assert_eq!(fs::read(target).unwrap(), sentinel);
        }
        assert_eq!(fs::read(&unrelated).unwrap(), unrelated_contents);
    }

    let ordinary = entropy_command()
        .args([
            "--out-bed",
            output_dir.to_str().unwrap(),
            "--regions",
            "../tests/resources/entropy_test_regions.bed",
            "--prefix",
            "ordinary",
            "--filter-threshold",
            "0",
            "--force",
        ])
        .output()
        .unwrap();
    assert!(
        ordinary.status.success(),
        "ordinary prefix failed:\n{}",
        String::from_utf8_lossy(&ordinary.stderr)
    );
    assert!(output_dir.join("ordinary_regions.bed").is_file());
    assert!(output_dir.join("ordinary_windows.bedgraph").is_file());
    assert_eq!(fs::read(unrelated).unwrap(), unrelated_contents);
}

#[test]
fn entropy_regions_second_open_failure_preserves_first_target() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_dir = temp_dir.path().join("entropy");
    fs::create_dir(&output_dir).unwrap();

    // NAME_MAX is 255 bytes on the supported Unix platforms. The first
    // suffix brings this component to exactly 255 bytes, while the second
    // makes it 260 bytes and must fail before the first target is truncated.
    let prefix = "x".repeat(243);
    let regions_output = output_dir.join(format!("{prefix}_regions.bed"));
    let sentinel = b"existing scientific output\n";
    fs::write(&regions_output, sentinel).unwrap();

    let rejected = entropy_command()
        .args([
            "--out-bed",
            output_dir.to_str().unwrap(),
            "--regions",
            "../tests/resources/entropy_test_regions.bed",
            "--prefix",
            &prefix,
            "--filter-threshold",
            "0",
            "--force",
        ])
        .output()
        .unwrap();

    assert!(!rejected.status.success());
    assert_eq!(fs::read(regions_output).unwrap(), sentinel);
}

#[cfg(unix)]
#[test]
fn entropy_file_force_rejects_existing_and_dangling_symlinks() {
    use std::os::unix::fs::symlink;

    let temp_dir = tempfile::tempdir().unwrap();
    let sentinel = b"outside scientific output\n";
    let outside_output = temp_dir.path().join("outside.bed");
    fs::write(&outside_output, sentinel).unwrap();
    let linked_output = temp_dir.path().join("linked.bed");
    symlink(&outside_output, &linked_output).unwrap();

    let run = |output_path: &Path| {
        entropy_command()
            .args([
                "--out-bed",
                output_path.to_str().unwrap(),
                "--filter-threshold",
                "0",
                "--force",
            ])
            .output()
            .unwrap()
    };

    let linked = run(&linked_output);
    assert!(!linked.status.success());
    assert!(
        String::from_utf8_lossy(&linked.stderr).contains("regular file"),
        "unexpected stderr:\n{}",
        String::from_utf8_lossy(&linked.stderr)
    );
    assert_eq!(fs::read(&outside_output).unwrap(), sentinel);
    assert!(fs::symlink_metadata(&linked_output)
        .unwrap()
        .file_type()
        .is_symlink());

    let missing_output = temp_dir.path().join("missing.bed");
    let dangling_output = temp_dir.path().join("dangling.bed");
    symlink(&missing_output, &dangling_output).unwrap();
    let dangling = run(&dangling_output);
    assert!(!dangling.status.success());
    assert!(
        String::from_utf8_lossy(&dangling.stderr).contains("regular file"),
        "unexpected stderr:\n{}",
        String::from_utf8_lossy(&dangling.stderr)
    );
    assert!(!missing_output.exists());
    assert!(fs::symlink_metadata(dangling_output)
        .unwrap()
        .file_type()
        .is_symlink());
}

#[cfg(unix)]
#[test]
fn entropy_regions_force_rejects_symlinks_before_pair_mutation() {
    use std::os::unix::fs::symlink;

    let temp_dir = tempfile::tempdir().unwrap();
    let output_dir = temp_dir.path().join("entropy");
    fs::create_dir(&output_dir).unwrap();
    let sentinel = b"outside scientific output\n";

    let run = |prefix: &str| {
        entropy_command()
            .args([
                "--out-bed",
                output_dir.to_str().unwrap(),
                "--regions",
                "../tests/resources/entropy_test_regions.bed",
                "--prefix",
                prefix,
                "--filter-threshold",
                "0",
                "--force",
            ])
            .output()
            .unwrap()
    };

    let outside_regions = temp_dir.path().join("outside-regions.bed");
    let outside_windows = temp_dir.path().join("outside-windows.bedgraph");
    fs::write(&outside_regions, sentinel).unwrap();
    fs::write(&outside_windows, sentinel).unwrap();
    let linked_regions = output_dir.join("linked_regions.bed");
    let linked_windows = output_dir.join("linked_windows.bedgraph");
    symlink(&outside_regions, &linked_regions).unwrap();
    symlink(&outside_windows, &linked_windows).unwrap();

    let linked = run("linked");
    assert!(!linked.status.success());
    assert!(String::from_utf8_lossy(&linked.stderr).contains("regular file"));
    assert_eq!(fs::read(&outside_regions).unwrap(), sentinel);
    assert_eq!(fs::read(&outside_windows).unwrap(), sentinel);
    assert!(fs::symlink_metadata(&linked_regions)
        .unwrap()
        .file_type()
        .is_symlink());
    assert!(fs::symlink_metadata(&linked_windows)
        .unwrap()
        .file_type()
        .is_symlink());

    let first_output = output_dir.join("second_regions.bed");
    let outside_second = temp_dir.path().join("outside-second.bedgraph");
    fs::write(&first_output, sentinel).unwrap();
    fs::write(&outside_second, sentinel).unwrap();
    symlink(&outside_second, output_dir.join("second_windows.bedgraph"))
        .unwrap();
    let second = run("second");
    assert!(!second.status.success());
    assert_eq!(fs::read(&first_output).unwrap(), sentinel);
    assert_eq!(fs::read(&outside_second).unwrap(), sentinel);

    let dangling_first = output_dir.join("dangling_regions.bed");
    let missing_second = temp_dir.path().join("missing-second.bedgraph");
    fs::write(&dangling_first, sentinel).unwrap();
    symlink(&missing_second, output_dir.join("dangling_windows.bedgraph"))
        .unwrap();
    let dangling = run("dangling");
    assert!(!dangling.status.success());
    assert_eq!(fs::read(dangling_first).unwrap(), sentinel);
    assert!(!missing_second.exists());
}

#[test]
fn entropy_stdout_output_remains_available() {
    let output =
        entropy_command().args(["--filter-threshold", "0"]).output().unwrap();
    assert!(
        output.status.success(),
        "stdout output failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!output.stdout.is_empty());
}
