use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use common::regional_query::{
    make_indexed_bedmethyl, OUTPUT_SENTINEL, REGIONS,
};

mod common;

fn run_stats(
    bedmethyl: &Path,
    regions: &Path,
    output: &Path,
    threads: usize,
    force: bool,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command
        .arg("stats")
        .arg(bedmethyl)
        .arg("--regions")
        .arg(regions)
        .arg("--out-table")
        .arg(output)
        .args(["--threads", &threads.to_string(), "--io-threads", "1"]);
    if force {
        command.arg("--force");
    }
    command.output().unwrap()
}

fn assert_query_failure(output: &Output) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "stats unexpectedly succeeded: {stderr}");
    assert!(
        stderr.contains("failed to calculate stats for region chr1:")
            && stderr.contains("invalid-bedmethyl-data"),
        "stats error lacks region/decode context: {stderr}"
    );
}

#[test]
fn stats_propagates_regional_query_errors_without_mutating_output() {
    for threads in [1usize, 4] {
        let temp_dir = tempfile::tempdir().unwrap();
        let bedmethyl = make_indexed_bedmethyl(temp_dir.path(), true);
        let regions = temp_dir.path().join("regions.bed");
        fs::write(&regions, REGIONS).unwrap();

        let absent_output = temp_dir.path().join("absent.tsv");
        assert_query_failure(&run_stats(
            &bedmethyl,
            &regions,
            &absent_output,
            threads,
            false,
        ));
        assert!(
            !absent_output.exists(),
            "failed stats run created output with {threads} threads"
        );

        let existing_output = temp_dir.path().join("existing.tsv");
        fs::write(&existing_output, OUTPUT_SENTINEL).unwrap();
        assert_query_failure(&run_stats(
            &bedmethyl,
            &regions,
            &existing_output,
            threads,
            true,
        ));
        assert_eq!(fs::read(existing_output).unwrap(), OUTPUT_SENTINEL);
    }
}

#[test]
fn stats_valid_output_is_byte_identical_across_thread_counts() {
    let temp_dir = tempfile::tempdir().unwrap();
    let bedmethyl = make_indexed_bedmethyl(temp_dir.path(), false);
    let regions = temp_dir.path().join("regions.bed");
    fs::write(&regions, REGIONS).unwrap();

    let one_thread = temp_dir.path().join("one.tsv");
    let many_threads = temp_dir.path().join("many.tsv");
    for (threads, output) in [(1usize, &one_thread), (4, &many_threads)] {
        let result = run_stats(&bedmethyl, &regions, output, threads, false);
        assert!(
            result.status.success(),
            "valid stats run failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    let expected = b"#chrom\tstart\tend\tname\tstrand\tcount_m\tcount_valid_m\tpercent_m\n\
chr1\t100\t101\t.\t.\t1\t4\t25\n\
chr1\t200\t201\t.\t.\t3\t6\t50\n\
chr1\t300\t301\t.\t.\t4\t8\t50\n\
chr1\t400\t401\t.\t.\t2\t10\t20\n";
    assert_eq!(fs::read(one_thread).unwrap(), expected);
    assert_eq!(fs::read(many_threads).unwrap(), expected);
}
