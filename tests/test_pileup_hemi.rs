use crate::common::{check_against_expected_text_file, run_modkit};

mod common;

#[test]
fn test_pileup_hemi_help() {
    let pileup_help_args = ["pileup-hemi", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

#[test]
fn test_pileup_hemi_hm() {
    let temp_file = std::env::temp_dir().join("test_pileup_hemi_hm.bed");
    let args = [
        "pileup-hemi",
        "tests/resources/duplex_modcalls_sort.bam",
        "-o",
        temp_file.to_str().unwrap(),
        "-r",
        "tests/resources/GRCh38_chr20.fa",
        "--motif",
        "CG",
        "0",
        "--region",
        "chr20:22,613,835-22,640,468",
        "--no-filtering",
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/duplex_hemi_nofilt.bed",
    );
}

#[test]
fn test_pileup_hemi_preset() {
    let temp_file = std::env::temp_dir().join("test_pileup_hemi_preset.bed");
    let args = [
        "pileup-hemi",
        "tests/resources/duplex_modcalls_sort.bam",
        "-o",
        temp_file.to_str().unwrap(),
        "-r",
        "tests/resources/GRCh38_chr20.fa",
        "--cpg",
        "--region",
        "chr20:22,613,835-22,640,468",
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/duplex_hemi.bed",
    );
}

// todo test with combine mods
