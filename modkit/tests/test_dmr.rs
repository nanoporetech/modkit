use crate::common::{
    check_against_expected_text_file, check_legal_csv, run_modkit,
};
use std::fs;
use std::path::Path;
use std::process::{Command, Output};

mod common;

const ZERO_COVERAGE_BED: &str =
    "../tests/resources/dmr_zero_coverage_chr1.bed.gz";
const CANONICAL_ONLY_BED: &str =
    "../tests/resources/dmr_canonical_only_chr1.bed.gz";
const MODIFIED_10X_BED: &str =
    "../tests/resources/dmr_10x_modified_chr1.bed.gz";
const HALF_MODIFIED_10X_BED: &str =
    "../tests/resources/dmr_10x_half_modified_chr1.bed.gz";
const CANONICAL_10X_BED: &str =
    "../tests/resources/dmr_10x_canonical_chr1.bed.gz";
const EMPTY_BED: &str = "../tests/resources/dmr_empty_chr1.bed.gz";
const ZERO_COVERAGE_REF: &str = "../tests/resources/dmr_zero_coverage_chr1.fa";
const ZERO_PERCENTILE_BED: &str =
    "../tests/resources/dmr_zero_percentile_21_sites_chr1.bed.gz";
const ZERO_PERCENTILE_REF: &str =
    "../tests/resources/dmr_zero_percentile_21_sites_chr1.fa";

#[test]
fn test_dmr_helps() {
    let _ = run_modkit(&["dmr", "pair", "--help"])
        .expect("failed to run modkit dmr pair help");
    let _ = run_modkit(&["dmr", "multi", "--help"])
        .expect("failed to run modkit dmr multi help");
}

#[test]
fn test_dmr_regression() {
    let out_bed = std::env::temp_dir().join("test_dmr_regression.bed");
    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "-f",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );

    let out_bed =
        std::env::temp_dir().join("foo").join("test_dmr_regression_2.bed");

    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "-f",
        "--header",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );
}

fn run_dmr_with_prior(output: &Path, alpha: &str, beta: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "dmr",
            "pair",
            "-a",
            "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "-b",
            "../tests/resources/lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "-o",
            output.to_str().unwrap(),
            "--ref",
            "../tests/resources/GRCh38_chr20.fa",
            "--base",
            "C",
            "--prior",
            alpha,
            beta,
            "--delta",
            "1",
            "--max-coverages",
            "100",
            "100",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
            "--force",
        ])
        .output()
        .unwrap()
}

#[test]
fn dmr_prior_cli_accepts_boundary_and_interior_but_rejects_invalid_inputs() {
    let temp_dir = tempfile::tempdir().unwrap();

    for (label, alpha, beta) in
        [("boundary", "0.5", "0.5"), ("interior", "0.55", "0.55")]
    {
        let output_path = temp_dir.path().join(format!("{label}.bed"));
        let output = run_dmr_with_prior(&output_path, alpha, beta);
        assert!(
            output.status.success(),
            "{label} prior ({alpha}, {beta}) was rejected: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let below_boundary = run_dmr_with_prior(
        &temp_dir.path().join("below-boundary.bed"),
        "0.4",
        "0.5",
    );
    assert!(!below_boundary.status.success());
    assert!(String::from_utf8_lossy(&below_boundary.stderr)
        .contains("alpha + beta must be >= 1.0 for numerical stability"));

    let non_positive =
        run_dmr_with_prior(&temp_dir.path().join("non-positive.bed"), "0", "1");
    assert!(!non_positive.status.success());
    assert!(String::from_utf8_lossy(&non_positive.stderr)
        .contains("invalid beta parameters 0, 1"));
}

fn run_dmr_with_max_coverages(
    output: &Path,
    max_coverages: [usize; 2],
    prior: Option<(&str, &str)>,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "dmr",
        "pair",
        "-a",
        MODIFIED_10X_BED,
        "-b",
        CANONICAL_10X_BED,
        "-o",
        output.to_str().unwrap(),
        "--ref",
        ZERO_COVERAGE_REF,
        "--base",
        "C",
        "--header",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--suppress-progress",
        "--force",
    ]);
    command
        .arg("--max-coverages")
        .arg(max_coverages[0].to_string())
        .arg(max_coverages[1].to_string());
    if let Some((alpha, beta)) = prior {
        command.args(["--prior", alpha, beta]);
    }
    command.output().unwrap()
}

fn assert_no_scientific_dmr_rows(output_path: &Path) {
    if output_path.exists() {
        let output = fs::read_to_string(output_path).unwrap();
        assert!(
            output.lines().all(|line| line.is_empty() || line.starts_with('#')),
            "unexpected scientific DMR output:\n{output}"
        );
    }
}

#[test]
fn explicit_zero_max_coverage_fails_clearly_without_scientific_rows() {
    let temp_dir = tempfile::tempdir().unwrap();

    for (prior_label, prior) in
        [("default", None), ("boundary", Some(("0.5", "0.5")))]
    {
        for max_coverages in [[0, 0], [0, 10], [10, 0]] {
            let output_path = temp_dir.path().join(format!(
                "{prior_label}-{}-{}.bed",
                max_coverages[0], max_coverages[1]
            ));
            let output =
                run_dmr_with_max_coverages(&output_path, max_coverages, prior);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(!output.status.success(), "{stderr}");
            assert!(
                stderr.contains(&format!(
                    "resolved maximum coverage must be greater than zero for \
                     both conditions, got control {} and experiment {}",
                    max_coverages[0], max_coverages[1]
                )),
                "{stderr}"
            );
            assert!(!stderr.contains("NaN"), "{stderr}");
            assert_no_scientific_dmr_rows(&output_path);
        }
    }
}

#[test]
fn explicit_positive_max_coverage_still_produces_exact_finite_output() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_path = temp_dir.path().join("positive-10-10.bed");

    let output = run_dmr_with_max_coverages(&output_path, [10, 10], None);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("processed 1 sites successfully, 0 failed"),
        "{stderr}"
    );

    let output = fs::read_to_string(output_path).unwrap();
    assert!(!output.contains("NaN"), "{output}");
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2, "{output}");
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let row = lines[1].split('\t').collect::<Vec<_>>();
    let field = |name| header.iter().position(|field| *field == name).unwrap();
    assert_eq!(row[field("a_counts")], "m:10");
    assert_eq!(row[field("a_total")], "10");
    assert_eq!(row[field("b_counts")], "m:0");
    assert_eq!(row[field("b_total")], "10");
    assert_eq!(row[field("a_pct_modified")], "1");
    assert_eq!(row[field("b_pct_modified")], "0");
    assert_eq!(row[field("effect_size")], "1");
    assert_eq!(row[field("map_pvalue")], "0.0000006230948043897833");
}

fn run_automatic_zero_percentile_dmr(
    output: &Path,
    threads: &str,
    io_threads: &str,
    max_coverages: Option<[&str; 2]>,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "dmr",
        "pair",
        "-a",
        ZERO_PERCENTILE_BED,
        "-b",
        ZERO_PERCENTILE_BED,
        "-o",
        output.to_str().unwrap(),
        "--ref",
        ZERO_PERCENTILE_REF,
        "--base",
        "C",
        "--header",
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ]);
    if let Some([control, experiment]) = max_coverages {
        command.args(["--max-coverages", control, experiment]);
    } else {
        command.args(["-N", "21", "--interval-size", "4"]);
    }
    command.output().unwrap()
}

#[test]
fn automatic_zero_percentile_cap_fails_and_explicit_positive_cap_recovers() {
    let temp_dir = tempfile::tempdir().unwrap();
    let mut error_lines = Vec::new();

    for (label, threads, io_threads) in
        [("one", "1", "1"), ("four-a", "4", "2"), ("four-b", "4", "2")]
    {
        let output_path = temp_dir.path().join(format!("{label}.bed"));
        let output = run_automatic_zero_percentile_dmr(
            &output_path,
            threads,
            io_threads,
            None,
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(!output.status.success(), "{stderr}");
        assert!(!stderr.contains("NaN"), "{stderr}");
        let error_line = stderr
            .lines()
            .find(|line| {
                line.contains(
                    "resolved maximum coverage must be greater than zero",
                )
            })
            .unwrap_or_else(|| {
                panic!("missing maximum-coverage error: {stderr}")
            });
        error_lines.push(error_line.to_string());
        assert_no_scientific_dmr_rows(&output_path);
    }

    assert!(error_lines.iter().all(|line| line == &error_lines[0]));
    assert!(error_lines[0].contains("control 0 and experiment 0"));

    let output_path = temp_dir.path().join("explicit-positive.bed");
    let output = run_automatic_zero_percentile_dmr(
        &output_path,
        "1",
        "1",
        Some(["10", "10"]),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let beta_error_line = stderr
        .lines()
        .find(|line| line.contains("beta-diff-calc-error"))
        .unwrap_or_else(|| panic!("missing beta-diff error count: {stderr}"));
    assert!(beta_error_line.contains("20"), "{stderr}");
    assert!(
        stderr.contains("processed 1 sites successfully, 20 failed"),
        "{stderr}"
    );
    assert!(!stderr.contains("NaN"), "{stderr}");

    let output = fs::read_to_string(output_path).unwrap();
    assert!(!output.contains("NaN"), "{output}");
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2, "{output}");
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let row = lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), row.len());
    let field = |name| header.iter().position(|field| *field == name).unwrap();
    assert_eq!(row[field("#chrom")], "chr1");
    assert_eq!(row[field("start")], "20");
    assert_eq!(row[field("end")], "21");
    assert_eq!(row[field("a_counts")], "C:0");
    assert_eq!(row[field("a_total")], "1");
    assert_eq!(row[field("b_counts")], "C:0");
    assert_eq!(row[field("b_total")], "1");
    assert_eq!(row[field("map_pvalue")], "1");
    assert_eq!(row[field("effect_size")], "0");
}

fn run_zero_coverage_dmr(
    output: &Path,
    prior: Option<(&str, &str)>,
    delta: Option<&str>,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "dmr",
        "pair",
        "-a",
        ZERO_COVERAGE_BED,
        "-b",
        CANONICAL_ONLY_BED,
        "-o",
        output.to_str().unwrap(),
        "--ref",
        ZERO_COVERAGE_REF,
        "--base",
        "C",
        "--max-coverages",
        "1",
        "1",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--suppress-progress",
        "--force",
    ]);
    if let Some((alpha, beta)) = prior {
        command.args(["--prior", alpha, beta]);
    }
    if let Some(delta) = delta {
        command.args(["--delta", delta]);
    }
    command.output().unwrap()
}

fn assert_zero_coverage_is_failed(output: Output, output_path: &Path) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "zero-coverage DMR site should be a recoverable site failure: {stderr}"
    );
    assert!(stderr.contains("beta-diff-calc-error"), "{stderr}");
    assert!(
        stderr.contains("processed 0 sites successfully, 1 failed"),
        "{stderr}"
    );

    let output_bytes = fs::read(output_path).unwrap();
    assert!(
        output_bytes.is_empty(),
        "zero-coverage site produced output: {}",
        String::from_utf8_lossy(&output_bytes)
    );
}

#[test]
fn zero_coverage_site_is_failed_without_nan_under_default_prior() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_path = temp_dir.path().join("default-prior.bed");

    let output = run_zero_coverage_dmr(&output_path, None, None);

    assert_zero_coverage_is_failed(output, &output_path);
}

#[test]
fn zero_coverage_site_does_not_abort_at_boundary_prior() {
    let temp_dir = tempfile::tempdir().unwrap();

    for (label, delta) in [("default-delta", None), ("delta-one", Some("1"))] {
        let output_path = temp_dir.path().join(format!("{label}.bed"));
        let output =
            run_zero_coverage_dmr(&output_path, Some(("0.5", "0.5")), delta);

        assert_zero_coverage_is_failed(output, &output_path);
    }
}

fn run_replicate_dmr(
    output_path: &Path,
    samples_a: &[&str],
    samples_b: &[&str],
    min_coverage: &str,
    threads: &str,
    io_threads: &str,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args(["dmr", "pair"]);
    for sample in samples_a {
        command.args(["-a", sample]);
    }
    for sample in samples_b {
        command.args(["-b", sample]);
    }
    command.args([
        "-o",
        output_path.to_str().unwrap(),
        "--ref",
        ZERO_COVERAGE_REF,
        "--base",
        "C",
        "--header",
        "--max-coverages",
        "10",
        "10",
        "--min-valid-coverage",
        min_coverage,
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ]);
    command.output().unwrap()
}

fn read_single_successful_dmr_row(
    output: Output,
    output_path: &Path,
) -> String {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("processed 1 sites successfully, 0 failed"),
        "{stderr}"
    );
    let output = fs::read_to_string(output_path).unwrap();
    assert_eq!(output.lines().count(), 2, "{output}");
    output
}

fn dmr_field<'a>(output: &'a str, field_name: &str) -> &'a str {
    let mut lines = output.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let row = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), row.len(), "{output}");
    let index = header.iter().position(|field| *field == field_name).unwrap();
    row[index]
}

#[test]
fn complete_matched_replicates_follow_cli_order_and_independent_oracles() {
    let temp_dir = tempfile::tempdir().unwrap();
    let samples_a =
        [MODIFIED_10X_BED, CANONICAL_10X_BED, HALF_MODIFIED_10X_BED];
    let samples_b =
        [CANONICAL_10X_BED, HALF_MODIFIED_10X_BED, MODIFIED_10X_BED];
    let expected_maps = [
        "0.0000006230948043897833",
        "0.020888505198542493",
        "0.01995939262165838",
    ];
    let expected_effects = ["1", "-0.5", "-0.5"];

    for (pair_index, ((sample_a, sample_b), (map, effect))) in samples_a
        .iter()
        .zip(samples_b.iter())
        .zip(expected_maps.iter().zip(expected_effects.iter()))
        .enumerate()
    {
        let output_path =
            temp_dir.path().join(format!("independent-{pair_index}.bed"));
        let output = run_replicate_dmr(
            &output_path,
            &[*sample_a],
            &[*sample_b],
            "0",
            "1",
            "1",
        );
        let output = read_single_successful_dmr_row(output, &output_path);
        assert_eq!(dmr_field(&output, "map_pvalue"), *map);
        assert_eq!(dmr_field(&output, "effect_size"), *effect);
    }

    let expected_maps = expected_maps.join(",");
    let expected_effects = expected_effects.join(",");
    let mut outputs = Vec::new();
    for (threads, io_threads) in [("1", "1"), ("4", "2")] {
        let output_path = temp_dir
            .path()
            .join(format!("complete-{threads}-{io_threads}.bed"));
        let output = run_replicate_dmr(
            &output_path,
            &samples_a,
            &samples_b,
            "0",
            threads,
            io_threads,
        );
        let output = read_single_successful_dmr_row(output, &output_path);
        assert_eq!(dmr_field(&output, "replicate_map_pvalues"), expected_maps);
        assert_eq!(
            dmr_field(&output, "replicate_effect_sizes"),
            expected_effects
        );
        for (field, expected) in [
            ("a_counts", "m:15"),
            ("a_total", "30"),
            ("b_counts", "m:15"),
            ("b_total", "30"),
            ("a_pct_modified", "0.5"),
            ("b_pct_modified", "0.5"),
            ("map_pvalue", "1"),
            ("effect_size", "0"),
            ("balanced_map_pvalue", "1"),
            ("balanced_effect_size", "0"),
            ("pct_a_samples", "100"),
            ("pct_b_samples", "100"),
        ] {
            assert_eq!(dmr_field(&output, field), expected);
        }
        outputs.push(output);
    }
    assert_eq!(outputs[0].as_bytes(), outputs[1].as_bytes());
}

#[test]
fn incomplete_matched_replicates_suppress_both_fields_and_keep_group_row() {
    let temp_dir = tempfile::tempdir().unwrap();
    let complete_b =
        [CANONICAL_10X_BED, HALF_MODIFIED_10X_BED, MODIFIED_10X_BED];
    let scenarios = [
        (
            "aligned-missing",
            [MODIFIED_10X_BED, EMPTY_BED, HALF_MODIFIED_10X_BED],
            [CANONICAL_10X_BED, EMPTY_BED, MODIFIED_10X_BED],
            "0",
        ),
        (
            "crossed-missing",
            [EMPTY_BED, CANONICAL_10X_BED, HALF_MODIFIED_10X_BED],
            [CANONICAL_10X_BED, EMPTY_BED, MODIFIED_10X_BED],
            "0",
        ),
        (
            "zero-member",
            [MODIFIED_10X_BED, ZERO_COVERAGE_BED, HALF_MODIFIED_10X_BED],
            complete_b,
            "0",
        ),
        (
            "threshold-filtered",
            [MODIFIED_10X_BED, CANONICAL_ONLY_BED, HALF_MODIFIED_10X_BED],
            complete_b,
            "2",
        ),
    ];

    for (label, samples_a, samples_b, min_coverage) in scenarios {
        let mut outputs = Vec::new();
        for (threads, io_threads) in [("1", "1"), ("4", "2")] {
            let output_path = temp_dir
                .path()
                .join(format!("{label}-{threads}-{io_threads}.bed"));
            let output = run_replicate_dmr(
                &output_path,
                &samples_a,
                &samples_b,
                min_coverage,
                threads,
                io_threads,
            );
            let output = read_single_successful_dmr_row(output, &output_path);
            assert_eq!(
                dmr_field(&output, "replicate_map_pvalues"),
                "-",
                "{label}: {output}"
            );
            assert_eq!(
                dmr_field(&output, "replicate_effect_sizes"),
                "-",
                "{label}: {output}"
            );
            outputs.push(output);
        }
        assert_eq!(outputs[0].as_bytes(), outputs[1].as_bytes(), "{label}");
    }
}

#[test]
fn positive_coverage_oracle_is_filter_and_thread_invariant() {
    let temp_dir = tempfile::tempdir().unwrap();
    let scenarios = [
        ("zero", [ZERO_COVERAGE_BED, MODIFIED_10X_BED], "0"),
        ("min-one", [ZERO_COVERAGE_BED, MODIFIED_10X_BED], "1"),
        ("absent", [EMPTY_BED, MODIFIED_10X_BED], "0"),
    ];
    let mut outputs = Vec::new();

    for (label, samples_a, min_coverage) in scenarios {
        for (threads, io_threads) in [("1", "1"), ("4", "2")] {
            let output_path = temp_dir
                .path()
                .join(format!("{label}-{threads}-{io_threads}.bed"));
            let output = run_replicate_dmr(
                &output_path,
                &samples_a,
                &[CANONICAL_10X_BED],
                min_coverage,
                threads,
                io_threads,
            );
            let output = read_single_successful_dmr_row(output, &output_path);
            for (field, expected) in [
                ("a_counts", "m:10"),
                ("a_total", "10"),
                ("b_counts", "m:0"),
                ("b_total", "10"),
                ("map_pvalue", "0.0000006230948043897833"),
                ("effect_size", "1"),
                ("balanced_map_pvalue", "0.0000006230948043897833"),
                ("balanced_effect_size", "1"),
                ("pct_a_samples", "50"),
                ("pct_b_samples", "100"),
            ] {
                assert_eq!(dmr_field(&output, field), expected);
            }
            assert!(!output.contains("h:0"), "{label}: {output}");
            outputs.push(output);
        }
    }

    for output in &outputs[1..] {
        assert_eq!(&outputs[0], output);
    }
}

fn run_segmented_dmr(
    output: &Path,
    segments: &Path,
    threads: &str,
    io_threads: &str,
) {
    run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        output.to_str().unwrap(),
        "--segment",
        segments.to_str().unwrap(),
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "--base",
        "C",
        "--max-coverages",
        "100",
        "100",
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ])
    .expect("segmented DMR run should succeed");
}

#[test]
fn segmentation_includes_last_site_and_is_thread_deterministic() {
    let temp_dir = tempfile::tempdir().unwrap();
    let sites_one = temp_dir.path().join("sites-threads-1.bed");
    let segments_one = temp_dir.path().join("segments-threads-1.bed");
    let sites_four = temp_dir.path().join("sites-threads-4.bed");
    let segments_four = temp_dir.path().join("segments-threads-4.bed");

    run_segmented_dmr(&sites_one, &segments_one, "1", "1");
    run_segmented_dmr(&sites_four, &segments_four, "4", "2");

    let sites_one = fs::read_to_string(sites_one).unwrap();
    let sites_four = fs::read_to_string(sites_four).unwrap();
    let segments_one = fs::read_to_string(segments_one).unwrap();
    let segments_four = fs::read_to_string(segments_four).unwrap();
    assert_eq!(sites_one.as_bytes(), sites_four.as_bytes());
    assert_eq!(segments_one.as_bytes(), segments_four.as_bytes());

    let site_rows = sites_one.lines().skip(1).collect::<Vec<_>>();
    let segment_rows = segments_one.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(site_rows.len(), 17_271);

    let segment_site_total = segment_rows
        .iter()
        .map(|row| row.split('\t').nth(5).unwrap().parse::<usize>().unwrap())
        .sum::<usize>();
    assert_eq!(segment_site_total, site_rows.len());

    let last_site_end = site_rows.last().unwrap().split('\t').nth(2).unwrap();
    let last_segment_end =
        segment_rows.last().unwrap().split('\t').nth(2).unwrap();
    assert_eq!(last_site_end, "10804378");
    assert_eq!(last_segment_end, last_site_end);
}

// todo
//  test pair with explicit index
//  test multi
