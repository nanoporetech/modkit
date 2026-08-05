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
