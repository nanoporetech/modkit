use crate::common::regional_query::{
    make_indexed_bedmethyl, GENOME_SIZES, OUTPUT_SENTINEL, REGIONS,
};
use crate::common::run_modkit;
use std::fs::{read_to_string, write};
use std::path::Path;
use std::process::{Command, Output};

mod common;

fn run_localize(
    regions: &Path,
    genome_sizes: &Path,
    output: &Path,
    force: bool,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "localize",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "--regions",
        regions.to_str().unwrap(),
        "--genome-sizes",
        genome_sizes.to_str().unwrap(),
        "--window",
        "10",
        "--min-coverage",
        "1",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--out-file",
        output.to_str().unwrap(),
    ]);
    if force {
        command.arg("--force");
    }
    command.output().unwrap()
}

fn run_localize_query(
    bedmethyl: &Path,
    regions: &Path,
    genome_sizes: &Path,
    output: &Path,
    threads: usize,
    force: bool,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command
        .arg("localize")
        .arg(bedmethyl)
        .arg("--regions")
        .arg(regions)
        .arg("--genome-sizes")
        .arg(genome_sizes)
        .arg("--out-file")
        .arg(output)
        .args([
            "--window",
            "1",
            "--min-coverage",
            "1",
            "--threads",
            &threads.to_string(),
            "--io-threads",
            "1",
        ]);
    if force {
        command.arg("--force");
    }
    command.output().unwrap()
}

fn assert_query_failure(output: &Output) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !output.status.success(),
        "localize unexpectedly succeeded: {stderr}"
    );
    assert!(
        stderr.contains("failed to localize region chr1:")
            && stderr.contains("invalid-bedmethyl-data"),
        "localize error lacks region/decode context: {stderr}"
    );
}

#[test]
fn test_localise_helps() {
    let _ = run_modkit(&["localize", "--help"])
        .expect("failed to run modkit localize help");
    let _ = run_modkit(&["localise", "--help"])
        .expect("failed to run modkit localise help");
}

#[test]
fn test_localize_rejects_mixed_regions_without_mutating_output() {
    let temp_dir = tempfile::tempdir().unwrap();
    let regions = temp_dir.path().join("mixed-regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    let new_output = temp_dir.path().join("new-localize.tsv");
    let existing_output = temp_dir.path().join("existing-localize.tsv");
    write(
        &regions,
        "# test regions\n\
         chr20\t9681998\t9681999\n\
         not-a-valid-bed-row\n",
    )
    .unwrap();
    write(&genome_sizes, "chr20\t100000000\n").unwrap();

    let create_result =
        run_localize(&regions, &genome_sizes, &new_output, false);
    let sentinel = "keep existing output content\n";
    write(&existing_output, sentinel).unwrap();
    let force_result =
        run_localize(&regions, &genome_sizes, &existing_output, true);

    let create_stderr = String::from_utf8_lossy(&create_result.stderr);
    assert!(!create_result.status.success(), "{create_stderr}");
    assert!(
        create_stderr.contains(regions.to_string_lossy().as_ref()),
        "failure did not reach region parsing: {create_stderr}"
    );
    assert!(create_stderr.contains("line 3"), "{create_stderr}");
    assert!(
        create_stderr.contains("delimiter mode changed"),
        "{create_stderr}"
    );
    assert!(!force_result.status.success());
    assert!(
        !new_output.exists()
            && read_to_string(&existing_output).unwrap() == sentinel,
        "failed region parsing mutated output paths: new_exists={}, \
         existing_content={:?}",
        new_output.exists(),
        read_to_string(existing_output).unwrap()
    );
}

#[test]
fn test_localize_accepts_valid_regions() {
    let temp_dir = tempfile::tempdir().unwrap();
    let regions = temp_dir.path().join("valid-regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    let output = temp_dir.path().join("localize.tsv");
    write(
        &regions,
        "# test regions\n\
         chr20\t9681998\t9681999\n\
         chr20\t9682013\t9682014\n",
    )
    .unwrap();
    write(&genome_sizes, "chr20\t100000000\n").unwrap();

    let result = run_localize(&regions, &genome_sizes, &output, false);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(read_to_string(output)
        .unwrap()
        .starts_with("mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n"));
}

#[test]
fn localize_propagates_regional_query_errors_without_mutating_output() {
    for threads in [1usize, 4] {
        let temp_dir = tempfile::tempdir().unwrap();
        let bedmethyl = make_indexed_bedmethyl(temp_dir.path(), true);
        let regions = temp_dir.path().join("regions.bed");
        let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
        write(&regions, REGIONS).unwrap();
        write(&genome_sizes, GENOME_SIZES).unwrap();

        let absent_output = temp_dir.path().join("absent.tsv");
        assert_query_failure(&run_localize_query(
            &bedmethyl,
            &regions,
            &genome_sizes,
            &absent_output,
            threads,
            false,
        ));
        assert!(
            !absent_output.exists(),
            "failed localize run created output with {threads} threads"
        );

        let existing_output = temp_dir.path().join("existing.tsv");
        write(&existing_output, OUTPUT_SENTINEL).unwrap();
        assert_query_failure(&run_localize_query(
            &bedmethyl,
            &regions,
            &genome_sizes,
            &existing_output,
            threads,
            true,
        ));
        assert_eq!(std::fs::read(existing_output).unwrap(), OUTPUT_SENTINEL);
    }
}

#[test]
fn localize_valid_output_is_byte_identical_across_thread_counts() {
    let temp_dir = tempfile::tempdir().unwrap();
    let bedmethyl = make_indexed_bedmethyl(temp_dir.path(), false);
    let regions = temp_dir.path().join("regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    write(&regions, REGIONS).unwrap();
    write(&genome_sizes, GENOME_SIZES).unwrap();

    let one_thread = temp_dir.path().join("one.tsv");
    let many_threads = temp_dir.path().join("many.tsv");
    for (threads, output) in [(1usize, &one_thread), (4, &many_threads)] {
        let result = run_localize_query(
            &bedmethyl,
            &regions,
            &genome_sizes,
            output,
            threads,
            false,
        );
        assert!(
            result.status.success(),
            "valid localize run failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    let expected = b"mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n\
m\t-1\t28\t10\t35.714287\n";
    assert_eq!(std::fs::read(one_thread).unwrap(), expected);
    assert_eq!(std::fs::read(many_threads).unwrap(), expected);
}

#[test]
fn test_localize_zero_window_includes_anchor() {
    let temp_dir = tempfile::tempdir().unwrap();
    let regions = temp_dir.path().join("regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    let output = temp_dir.path().join("localize.tsv");

    write(&regions, "chr20\t9681998\t9681999\n").unwrap();
    write(&genome_sizes, "chr20\t100000000\n").unwrap();

    run_modkit(&[
        "localize",
        "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "--regions",
        regions.to_str().unwrap(),
        "--genome-sizes",
        genome_sizes.to_str().unwrap(),
        "--window",
        "0",
        "--min-coverage",
        "1",
        "--out-file",
        output.to_str().unwrap(),
    ])
    .expect("failed to run modkit localize");

    assert_eq!(
        read_to_string(output).unwrap(),
        "mod_code\toffset\tn_valid\tn_mod\tpercent_modified\nC\t0\t1\t1\t100\n"
    );
}

#[test]
fn test_localize_offsets_use_reference_axis_for_all_feature_strands() {
    let temp_dir = tempfile::tempdir().unwrap();
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    write(&genome_sizes, "chr20\t100000000\n").unwrap();

    let expected = concat!(
        "mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n",
        "C\t-7\t1\t1\t100\n",
        "C\t8\t1\t1\t100\n",
    );
    let feature_strands = [
        ("positive", "chr20\t9682005\t9682006\tfeature\t0\t+\n"),
        ("negative", "chr20\t9682005\t9682006\tfeature\t0\t-\n"),
        ("both", "chr20\t9682005\t9682006\n"),
    ];

    for (label, region_line) in feature_strands {
        let regions = temp_dir.path().join(format!("regions-{label}.bed"));
        let output = temp_dir.path().join(format!("localize-{label}.tsv"));
        write(&regions, region_line).unwrap();

        run_modkit(&[
            "localize",
            "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "--regions",
            regions.to_str().unwrap(),
            "--genome-sizes",
            genome_sizes.to_str().unwrap(),
            "--window",
            "8",
            "--stranded-features",
            "both",
            "--min-coverage",
            "1",
            "--out-file",
            output.to_str().unwrap(),
        ])
        .expect("failed to run modkit localize");

        assert_eq!(read_to_string(output).unwrap(), expected, "{label}");
    }
}
