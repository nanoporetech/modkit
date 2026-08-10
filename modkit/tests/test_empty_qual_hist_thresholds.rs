use std::fs;
use std::path::Path;
use std::process::{Command, Output};

const MOD_BAM: &str = "../tests/resources/bc_anchored_10_reads.sorted.bam";
const EMPTY_TAGS_BAM: &str = "../tests/resources/empty-tags.sorted.bam";
const REFERENCE: &str = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
const EMPTY_THRESHOLD_ERROR: &str =
    "cannot calculate automatic thresholds because no modification \
     probabilities were sampled";

fn run_modkit(args: &[&str]) -> Output {
    Command::new(Path::new(env!("CARGO_BIN_EXE_modkit")))
        .args(args)
        .output()
        .unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_empty_auto_failure(output: &Output) {
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains(EMPTY_THRESHOLD_ERROR),
        "unexpected stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn pileup_rejects_empty_automatic_observations_but_controls_succeed() {
    let temp_dir = tempfile::tempdir().unwrap();
    let include_bed = temp_dir.path().join("empty_site.bed");
    fs::write(&include_bed, "oligo_1512_adapters\t0\t1\n").unwrap();

    let run = |output_path: &Path, threshold_args: &[&str]| {
        Command::new(Path::new(env!("CARGO_BIN_EXE_modkit")))
            .args([
                "pileup",
                MOD_BAM,
                output_path.to_str().unwrap(),
                "--include-bed",
                include_bed.to_str().unwrap(),
                "--modified-bases",
                "C:m",
                "--reference",
                REFERENCE,
                "--threads",
                "1",
                "--sampling-threads",
                "1",
                "--io-threads",
                "1",
                "--suppress-progress",
            ])
            .args(threshold_args)
            .output()
            .unwrap()
    };

    let automatic_path = temp_dir.path().join("automatic.bed");
    let automatic = run(&automatic_path, &[]);
    assert_empty_auto_failure(&automatic);

    for (name, threshold_args) in [
        ("explicit.bed", vec!["--filter-threshold", "0"]),
        ("no_filter.bed", vec!["--no-filtering"]),
    ] {
        let output_path = temp_dir.path().join(name);
        let output = run(&output_path, &threshold_args);
        assert_success(&output);
        assert_eq!(fs::metadata(output_path).unwrap().len(), 0);
    }
}

#[test]
fn summary_rejects_empty_automatic_observations_but_explicit_succeeds() {
    let automatic = run_modkit(&[
        "summary",
        EMPTY_TAGS_BAM,
        "--num-reads",
        "10",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--tsv",
        "--suppress-progress",
    ]);
    assert_empty_auto_failure(&automatic);
    assert!(automatic.stdout.is_empty());

    let explicit = run_modkit(&[
        "summary",
        EMPTY_TAGS_BAM,
        "--num-reads",
        "10",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--filter-threshold",
        "0",
        "--tsv",
        "--suppress-progress",
    ]);
    assert_success(&explicit);
    assert!(String::from_utf8_lossy(&explicit.stdout)
        .contains("total_reads_used\t0"));
}

#[test]
fn extract_rejects_empty_automatic_observations_but_controls_succeed() {
    let temp_dir = tempfile::tempdir().unwrap();
    let include_bed = temp_dir.path().join("empty_site.bed");
    fs::write(&include_bed, "oligo_1512_adapters\t0\t1\n").unwrap();

    let run = |output_path: &Path, threshold_args: &[&str]| {
        Command::new(Path::new(env!("CARGO_BIN_EXE_modkit")))
            .args([
                "extract",
                "calls",
                MOD_BAM,
                output_path.to_str().unwrap(),
                "--include-bed",
                include_bed.to_str().unwrap(),
                "--threads",
                "1",
                "--io-threads",
                "1",
                "--sample-num-reads",
                "10",
                "--suppress-progress",
                "--no-headers",
            ])
            .args(threshold_args)
            .output()
            .unwrap()
    };

    let automatic_path = temp_dir.path().join("automatic.tsv");
    let automatic = run(&automatic_path, &[]);
    assert_empty_auto_failure(&automatic);
    assert!(!automatic_path.exists());

    for (name, threshold_args) in [
        ("explicit.tsv", vec!["--filter-threshold", "0"]),
        ("no_filter.tsv", vec!["--no-filtering"]),
    ] {
        let output_path = temp_dir.path().join(name);
        let output = run(&output_path, &threshold_args);
        assert_success(&output);
        assert_eq!(fs::metadata(output_path).unwrap().len(), 0);
    }
}
