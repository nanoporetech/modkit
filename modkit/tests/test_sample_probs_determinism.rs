use std::{collections::BTreeMap, fs, path::Path, process::Command};

use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};

fn run_sample_probs(
    out_dir: &Path,
    prefix: &str,
    input: &str,
    reference: Option<&str>,
    fraction: Option<&str>,
    seed: u64,
    threads: usize,
    interval_size: u32,
) -> (Vec<u8>, Vec<u8>) {
    let mut args = vec![
        "sample-probs".to_string(),
        input.to_string(),
        "--threads".to_string(),
        threads.to_string(),
        "--interval-size".to_string(),
        interval_size.to_string(),
        "--hist".to_string(),
        "--out-dir".to_string(),
        out_dir.to_str().unwrap().to_string(),
        "--prefix".to_string(),
        prefix.to_string(),
        "--suppress-progress".to_string(),
    ];
    if let Some(reference) = reference {
        args.extend(["--reference".to_string(), reference.to_string()]);
    }
    if let Some(fraction) = fraction {
        args.extend([
            "--sample-frac".to_string(),
            fraction.to_string(),
            "--seed".to_string(),
            seed.to_string(),
        ]);
    } else {
        args.push("--no-sampling".to_string());
    }
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args(args)
        .output()
        .expect("failed to run modkit sample-probs");
    assert!(
        output.status.success(),
        "sample-probs failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let thresholds =
        fs::read(out_dir.join(format!("{prefix}_thresholds.tsv"))).unwrap();
    let probabilities =
        fs::read(out_dir.join(format!("{prefix}_probabilities.tsv"))).unwrap();
    (thresholds, probabilities)
}

fn run_summary(threads: usize, interval_size: u32) -> BTreeMap<String, String> {
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "summary",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
            "--sampling-frac",
            "0.5",
            "--seed",
            "7",
            "--threads",
            &threads.to_string(),
            "--interval-size",
            &interval_size.to_string(),
            "--tsv",
            "--suppress-progress",
        ])
        .output()
        .expect("failed to run modkit summary");
    assert!(
        output.status.success(),
        "summary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    // Summary currently has a separate, known nondeterministic row-order issue
    // (U46), so compare its complete key/value row map while rejecting duplicate
    // keys instead of weakening this test to a sorted bag of lines.
    let stdout = String::from_utf8(output.stdout).unwrap();
    let mut rows = BTreeMap::new();
    for line in stdout.lines() {
        let (key, value) = line
            .split_once('\t')
            .unwrap_or_else(|| panic!("invalid summary TSV row: {line}"));
        assert!(
            rows.insert(key.to_string(), value.to_string()).is_none(),
            "duplicate summary TSV key: {key}"
        );
    }
    for required_key in ["mod_bases", "total_reads_used"] {
        assert!(
            rows.contains_key(required_key),
            "summary TSV is missing metadata key: {required_key}"
        );
    }
    rows
}

fn run_pileup(
    out_dir: &Path,
    name: &str,
    threads: usize,
    sampling_interval_size: u32,
) -> Vec<u8> {
    let output_path = out_dir.join(format!("{name}.bed"));
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "pileup",
            "-i",
            "25",
            "--sampling-frac",
            "0.5",
            "--filter-percentile",
            "0.25",
            "--seed",
            "7",
            "--threads",
            &threads.to_string(),
            "--sampling-threads",
            &threads.to_string(),
            "--sampling-interval-size",
            &sampling_interval_size.to_string(),
            "--suppress-progress",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run modkit pileup");
    assert!(
        output.status.success(),
        "pileup failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read(output_path).unwrap()
}

fn help_text(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args(args)
        .output()
        .expect("failed to run modkit help");
    assert!(
        output.status.success(),
        "help command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

#[test]
fn seeded_indexed_sampling_is_stable_across_worker_and_interval_matrix() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input = "../tests/resources/bc_anchored_10_reads.sorted.bam";
    let mut reader = bam::Reader::from_path(input).unwrap();
    assert!(reader.records().any(|result| {
        let record = result.unwrap();
        record.reference_end() - record.pos() > 20
    }));
    let baseline = run_sample_probs(
        temp_dir.path(),
        "w1_i20",
        input,
        None,
        Some("0.5"),
        7,
        1,
        20,
    );

    for (prefix, threads, interval_size) in [
        ("w1_i20_repeat", 1, 20),
        ("w2_i20", 2, 20),
        ("w3_i20", 3, 20),
        ("w8_i20", 8, 20),
        ("w8_i20_repeat", 8, 20),
        ("w1_i100", 1, 100),
        ("w3_i100", 3, 100),
        ("w8_i1000", 8, 1_000),
    ] {
        let observed = run_sample_probs(
            temp_dir.path(),
            prefix,
            input,
            None,
            Some("0.5"),
            7,
            threads,
            interval_size,
        );
        assert!(
            observed == baseline,
            "seeded threshold or probability output changed for {threads} \
             workers and interval size {interval_size}"
        );
    }

    let distinct_seed = run_sample_probs(
        temp_dir.path(),
        "seed8",
        input,
        None,
        Some("0.5"),
        8,
        3,
        20,
    );
    assert!(
        distinct_seed != baseline,
        "different explicit seeds unexpectedly produced identical outputs"
    );
}

#[test]
fn seeded_indexed_sampling_matches_bam_and_selective_cram_decoding() {
    let temp_dir = tempfile::tempdir().unwrap();
    let bam = run_sample_probs(
        temp_dir.path(),
        "bam",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        None,
        Some("0.5"),
        7,
        3,
        20,
    );
    let cram = run_sample_probs(
        temp_dir.path(),
        "cram",
        "../tests/resources/bc_anchored_10_reads.sorted.cram",
        Some("../tests/resources/CGI_ladder_3.6kb_ref.fa"),
        Some("0.5"),
        7,
        3,
        20,
    );

    assert!(cram == bam, "seeded BAM and CRAM outputs differed");
}

#[test]
fn boundary_fractions_exclude_none_or_all_exactly() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input = "../tests/resources/bc_anchored_10_reads.sorted.bam";
    let all = run_sample_probs(
        temp_dir.path(),
        "all",
        input,
        None,
        Some("1"),
        7,
        1,
        20,
    );
    let no_sampling = run_sample_probs(
        temp_dir.path(),
        "no_sampling",
        input,
        None,
        None,
        0,
        8,
        1_000,
    );
    assert!(
        all == no_sampling,
        "fraction 1 did not match the all-records control"
    );

    let none_serial = run_sample_probs(
        temp_dir.path(),
        "none_serial",
        input,
        None,
        Some("0"),
        7,
        1,
        20,
    );
    let none_parallel = run_sample_probs(
        temp_dir.path(),
        "none_parallel",
        input,
        None,
        Some("0"),
        7,
        8,
        1_000,
    );
    assert!(
        none_serial == none_parallel,
        "fraction 0 changed across scheduling"
    );
    assert_eq!(
        String::from_utf8(none_serial.1).unwrap().lines().count(),
        1,
        "fraction 0 emitted probability observations"
    );
}

#[test]
fn seeded_summary_and_pileup_consumers_are_schedule_independent() {
    assert_eq!(run_summary(1, 20), run_summary(8, 1_000));

    let temp_dir = tempfile::tempdir().unwrap();
    assert_eq!(
        run_pileup(temp_dir.path(), "serial", 1, 20),
        run_pileup(temp_dir.path(), "parallel", 8, 1_000)
    );
}

#[test]
fn indexed_sampling_seed_help_describes_repeatable_decisions() {
    for args in [["sample-probs", "--help"], ["call-mods", "--help"]] {
        let help = help_text(&args);
        assert!(help.contains(
            "Provide a seed to make fractional read sampling decisions \
             repeatable for indexed and unindexed inputs"
        ));
        assert!(!help.contains("only used when no BAM index is provided"));
    }
}
