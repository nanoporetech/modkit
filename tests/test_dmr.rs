use crate::common::{check_against_expected_text_file, run_modkit};

mod common;

#[test]
fn test_dmr_regression() {
    let out_bed = std::env::temp_dir().join("test_dmr_regression.bed");
    let _ = run_modkit( &[
        "dmr",
        "tests/resources/lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "tests/resources/lung_00733-M_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        out_bed.to_str().unwrap(),
        "-r",
        "tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "tests/resources/grch38_chr20.fa",
        "-f",
        "--base", "C",
        ]).expect("failed to run modkit dmr");

    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "/Users/art.rand/projects/mod_flatten/data/dmr-dev/test_output_chr20-2.bed",
    )
}
