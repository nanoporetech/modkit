use std::{fs, path::Path, process::Command};

fn sam_record_body(path: &Path) -> Vec<u8> {
    let sam = String::from_utf8(fs::read(path).unwrap()).unwrap();
    sam.lines()
        .filter(|line| !line.starts_with('@'))
        .collect::<Vec<_>>()
        .join("\n")
        .into_bytes()
}

fn run_call_mods(
    out_dir: &Path,
    name: &str,
    input: &str,
    seed: u64,
    threads: usize,
    sampling_interval_size: u32,
) -> Vec<u8> {
    let output_path = out_dir.join(format!("{name}.sam"));
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "call-mods",
            input,
            output_path.to_str().unwrap(),
            "--output-sam",
            "--sampling-frac",
            "0.5",
            "--seed",
            &seed.to_string(),
            "--threads",
            &threads.to_string(),
            "--sampling-interval-size",
            &sampling_interval_size.to_string(),
            "--filter-percentile",
            "0.25",
            "--suppress-progress",
        ])
        .output()
        .expect("failed to run modkit call-mods");
    assert!(
        output.status.success(),
        "call-mods failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    sam_record_body(&output_path)
}

fn run_pileup_hemi(
    out_dir: &Path,
    name: &str,
    threads: usize,
    sampling_interval_size: u32,
) -> (Vec<u8>, String) {
    let output_path = out_dir.join(format!("{name}.bed"));
    let log_path = out_dir.join(format!("{name}.log"));
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "pileup-hemi",
            "../tests/resources/duplex_modcalls_sort.bam",
            "--out-bed",
            output_path.to_str().unwrap(),
            "--ref",
            "../tests/resources/GRCh38_chr20.fa",
            "--cpg",
            "--region",
            "chr20:22,613,835-22,640,468",
            "--sampling-frac",
            "0.5",
            "--seed",
            "7",
            "--threads",
            &threads.to_string(),
            "--sampling-interval-size",
            &sampling_interval_size.to_string(),
            "--suppress-progress",
            "--log-filepath",
            log_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run modkit pileup-hemi");
    assert!(
        output.status.success(),
        "pileup-hemi failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let log = fs::read_to_string(log_path).unwrap();
    let threshold = log
        .lines()
        .find_map(|line| {
            line.find("estimated pass threshold")
                .map(|start| line[start..].to_string())
        })
        .expect("pileup-hemi log did not contain an estimated threshold");
    (fs::read(output_path).unwrap(), threshold)
}

#[test]
fn call_mods_seeded_fraction_is_stable_across_workers_and_intervals() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input = "../tests/resources/bc_anchored_10_reads.sorted.bam";
    let baseline = run_call_mods(temp_dir.path(), "baseline", input, 7, 1, 20);

    for (name, threads, interval_size) in [
        ("repeat", 1, 20),
        ("more_workers", 8, 20),
        ("larger_intervals", 8, 1_000),
    ] {
        assert_eq!(
            run_call_mods(
                temp_dir.path(),
                name,
                input,
                7,
                threads,
                interval_size
            ),
            baseline,
            "call-mods output changed for {threads} workers and interval size \
             {interval_size}"
        );
    }

    assert_ne!(
        run_call_mods(temp_dir.path(), "different_seed", input, 8, 3, 20),
        baseline,
        "call-mods ignored the explicit fractional-sampling seed"
    );
}

#[test]
fn call_mods_unindexed_seeded_streaming_behavior_remains_repeatable() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input = temp_dir.path().join("unindexed.bam");
    fs::copy("../tests/resources/bc_anchored_10_reads.sorted.bam", &input)
        .unwrap();
    let input = input.to_str().unwrap();
    let baseline =
        run_call_mods(temp_dir.path(), "unindexed_a", input, 7, 1, 20);

    assert_eq!(
        run_call_mods(temp_dir.path(), "unindexed_b", input, 7, 8, 1_000),
        baseline
    );
    assert_ne!(
        run_call_mods(temp_dir.path(), "unindexed_seed8", input, 8, 3, 20),
        baseline
    );
}

#[test]
fn duplex_threshold_sampling_is_stable_across_workers_and_intervals() {
    let temp_dir = tempfile::tempdir().unwrap();
    let baseline = run_pileup_hemi(temp_dir.path(), "baseline", 1, 1_000);

    for (name, threads, interval_size) in
        [("repeat", 1, 1_000), ("matrix", 8, 100_000)]
    {
        assert_eq!(
            run_pileup_hemi(temp_dir.path(), name, threads, interval_size),
            baseline,
            "pileup-hemi threshold or output changed for {threads} workers \
             and interval size {interval_size}"
        );
    }
}
