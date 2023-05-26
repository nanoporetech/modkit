use rust_htslib::bam;
use std::fs::File;
use std::io::{BufRead, BufReader};

use common::{check_against_expected_text_file, run_modkit};

mod common;

#[test]
fn test_help() {
    let pileup_help_args = ["pileup", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

#[test]
fn test_mod_pileup_no_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_nofilt.methyl.bed",
    );
}

#[test]
fn test_mod_pileup_with_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_withfilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-f",
        "1.0",
        "-p",
        "0.25",
        "--only-tabs",
        "--seed",
        "42",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_filt025.methyl.bed",
    );
}

#[test]
fn test_mod_pileup_combine() {
    let test_adjusted_bam = std::env::temp_dir().join("test_combined.bed");
    let pileup_args = [
        "pileup",
        "--combine-mods",
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_adjusted_bam.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_adjusted_bam.exists());

    check_against_expected_text_file(
        test_adjusted_bam.to_str().unwrap(),
        "tests/resources/modbam.modpileup_combined.methyl.bed",
    );
}

#[test]
fn test_mod_pileup_collapse() {
    let test_collapsed_bam = std::env::temp_dir().join("test_collapsed.bam");
    let test_collapsed_bed = std::env::temp_dir().join("test_collapsed.bed");
    let test_restricted_bed = std::env::temp_dir().join("test_restricted.bed");

    let collapse_args = [
        "adjust-mods",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_collapsed_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();
    assert!(test_collapsed_bam.exists());
    bam::index::build(
        test_collapsed_bam.clone(),
        None,
        bam::index::Type::Bai,
        1,
    )
    .unwrap();

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        test_collapsed_bam.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_collapsed_bed.exists());

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--ignore",
        "h",
        "--no-filtering",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_restricted_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_restricted_bed.exists());
    check_against_expected_text_file(
        test_restricted_bed.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_no_mod_calls() {
    let empty_bedfile =
        std::env::temp_dir().join("test_pileup_no_mod_calls_outbed.bed");
    let args = [
        "pileup",
        "--no-filtering",
        "tests/resources/empty-tags.sorted.bam",
        empty_bedfile.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    let reader = BufReader::new(File::open(empty_bedfile).unwrap());
    let lines = reader.lines().collect::<Vec<Result<String, _>>>();
    assert_eq!(lines.len(), 0);
}

#[test]
fn test_pileup_old_tags() {
    let updated_file =
        std::env::temp_dir().join("test_pileup_old_tags_updated.bam");
    run_modkit(&[
        "update-tags",
        "tests/resources/HG002_small.ch20._other.sorted.bam",
        "--mode",
        "ambiguous",
        updated_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(updated_file.exists());
    bam::index::build(updated_file.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    let out_file = std::env::temp_dir().join("test_pileup_old_tags.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        "--only-tabs",
        updated_file.to_str().unwrap(),
        out_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(out_file.exists());
    check_against_expected_text_file(
        out_file.to_str().unwrap(),
        "tests/resources/pileup-old-tags-regressiontest.methyl.bed",
    );
}

#[test]
fn test_pileup_with_region() {
    let temp_file = std::env::temp_dir().join("test_pileup_with_region.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--region",
        "oligo_1512_adapters:0-50",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_nofilt_oligo_1512_adapters_10_50.bed",
    );
}

#[test]
fn test_pileup_duplex_reads() {
    let temp_file = std::env::temp_dir().join("test_pileup_duplex_reads.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/duplex_modbam.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/duplex_modbam_pileup_nofilt.bed",
    );
}

#[test]
fn test_cpg_motif_filtering() {
    let temp_file = std::env::temp_dir().join("test_cpg_motif_filtering.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_cpg_motif_filtering_strand_combine() {
    let temp_file = std::env::temp_dir()
        .join("test_cpg_motif_filtering_strand_combine.bed");
    for interval_size in
        ["10", "88", "89", "90", "91", "92", "93", "94", "10000"]
    {
        run_modkit(&[
            "pileup",
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            temp_file.to_str().unwrap(),
            "--no-filtering",
            "-i",
            interval_size,
            "--cpg",
            "--combine-strands",
            "--ref",
            "tests/resources/CGI_ladder_3.6kb_ref.fa",
        ])
        .unwrap();
        check_against_expected_text_file(
            temp_file.to_str().unwrap(),
            "tests/resources/bc_anchored_10_reads_nofilt_cg_motif_strand_combine.bed",
        );
    }
}

#[test]
fn test_presets_traditional_same_as_options() {
    let preset_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options.bed");
    let options_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options2.bed");

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        preset_temp_file.to_str().unwrap(),
        "--no-filtering",
        "--preset",
        "traditional",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        options_temp_file.to_str().unwrap(),
        "--cpg",
        "--no-filtering",
        "--ignore",
        "h",
        "--combine-strands",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        preset_temp_file.to_str().unwrap(),
        options_temp_file.to_str().unwrap(),
    );
}

#[test]
fn test_duplicated_reads_ignored() {
    let control_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_control.bed");
    let test_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_marked.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        control_fp.to_str().unwrap(),
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/duplicated.marked.fixed.bam",
        test_fp.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        control_fp.to_str().unwrap(),
        test_fp.to_str().unwrap(),
    );
}
