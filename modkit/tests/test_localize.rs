use crate::common::run_modkit;
use std::fs;
use tempfile::tempdir;

mod common;

#[test]
fn test_localise_helps() {
    let _ = run_modkit(&["localize", "--help"])
        .expect("failed to run modkit localize help");
    let _ = run_modkit(&["localise", "--help"])
        .expect("failed to run modkit localise help");
}

#[test]
fn test_localize_min_coverage_filters_before_aggregation() {
    let temp_dir = tempdir().unwrap();
    let regions = temp_dir.path().join("regions.bed");
    let genome_sizes = temp_dir.path().join("genome-sizes.tsv");
    fs::write(&regions, "chr20\t9681998\t9681999\nchr20\t9838537\t9838538\n")
        .unwrap();
    fs::write(&genome_sizes, "chr20\t64444167\n").unwrap();

    let run = |min_coverage: u64| {
        let output = temp_dir.path().join(format!("min-{min_coverage}.tsv"));
        let min_coverage = min_coverage.to_string();
        run_modkit(&[
            "localize",
            "../tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
            "--regions",
            regions.to_str().unwrap(),
            "--genome-sizes",
            genome_sizes.to_str().unwrap(),
            "--window",
            "1",
            "--min-coverage",
            &min_coverage,
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--out-file",
            output.to_str().unwrap(),
        ])
        .unwrap();
        fs::read_to_string(output).unwrap().replace("\r\n", "\n")
    };

    assert_eq!(
        run(1),
        "mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n\
         C\t-1\t24\t3\t12.5\n"
    );
    assert_eq!(
        run(3),
        "mod_code\toffset\tn_valid\tn_mod\tpercent_modified\n\
         C\t-1\t23\t2\t8.695652\n"
    );
}
