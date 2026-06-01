use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
};

use common::run_modkit;
use mod_kit::dmr::bedmethyl::BedMethylLine;

mod common;

#[test]
fn test_merge_bedmethyl_help() {
    run_modkit(&["bedmethyl", "--help"]).unwrap();
    run_modkit(&["bedmethyl", "merge", "--help"]).unwrap();
    run_modkit(&["bedmethyl", "tobigwig", "--help"]).unwrap();
}

#[test]
fn test_bedmethyl_merge() {
    let sizes = "chr20\t64444167\n";
    let sizes_fp = std::env::temp_dir().join("test_merge_bedmethyl_sizes.tsv");
    {
        let w = File::create(&sizes_fp).unwrap();
        let mut w = BufWriter::new(w);
        w.write(sizes.as_bytes()).unwrap();
        drop(w);
    };
    assert!(sizes_fp.exists());
    let bed_fp = "../tests/resources/\
                  lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.\
                  gz";
    let out_bed = std::env::temp_dir().join("test_merge_bedmethyl.bed");

    run_modkit(&[
        "bedmethyl",
        "merge",
        bed_fp,
        bed_fp,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_bed.to_str().unwrap(),
        "--force",
    ])
    .unwrap();
    assert!(out_bed.exists());
    // dbg!(&out_bed);

    let mut buff = vec![];
    let _ = rust_htslib::bgzf::Reader::from_path(bed_fp)
        .unwrap()
        .read_to_end(&mut buff);
    let input_records = buff
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .collect::<Vec<BedMethylLine>>();
    let merged_records = BufReader::new(File::open(out_bed).unwrap())
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .collect::<Vec<BedMethylLine>>();
    assert_eq!(input_records.len(), merged_records.len());

    for (x, y) in input_records.into_iter().zip(merged_records) {
        assert_eq!(x.chrom, y.chrom);
        assert_eq!(x.start(), y.start());
        assert_eq!(x.raw_mod_code, y.raw_mod_code);
        assert_eq!(x.strand, y.strand);
        assert_eq!(x.count_methylated * 2, y.count_methylated);
        assert_eq!(x.valid_coverage * 2, y.valid_coverage);
        assert_eq!(x.count_canonical * 2, y.count_canonical);
        assert_eq!(x.count_other * 2, y.count_other);
        assert_eq!(x.count_delete * 2, y.count_delete);
        assert_eq!(x.count_fail * 2, y.count_fail);
        assert_eq!(x.count_diff * 2, y.count_diff);
        assert_eq!(x.count_nocall * 2, y.count_nocall);
    }
}

// Shared setup helpers for the --min-samples / --min-sample-coverage tests.
// The committed chr20 pileup resource is passed as multiple inputs: passing it
// N times makes every position present in exactly N inputs, which lets us
// exercise the inner-join threshold without building new fixtures.
fn min_samples_test_sizes() -> std::path::PathBuf {
    let sizes = "chr20\t64444167\n";
    let sizes_fp =
        std::env::temp_dir().join("test_merge_min_samples_sizes.tsv");
    let mut w = BufWriter::new(File::create(&sizes_fp).unwrap());
    w.write(sizes.as_bytes()).unwrap();
    drop(w);
    sizes_fp
}

const MIN_SAMPLES_BED_FP: &str = "../tests/resources/\
    lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz";

fn count_bed_records(fp: &std::path::Path) -> usize {
    BufReader::new(File::open(fp).unwrap())
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .count()
}

#[test]
fn test_bedmethyl_merge_min_samples() {
    let sizes_fp = min_samples_test_sizes();

    // Decode the input once to know how many positions it has.
    let mut buff = vec![];
    rust_htslib::bgzf::Reader::from_path(MIN_SAMPLES_BED_FP)
        .unwrap()
        .read_to_end(&mut buff)
        .expect("failed to read MIN_SAMPLES_BED_FP fixture");
    let n_positions = buff
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .count();

    // --min-samples 2 with the file passed twice: every position is present in
    // both inputs, so all positions are kept (an inner join that retains all).
    let out_keep = std::env::temp_dir().join("test_merge_min_samples_keep.bed");
    run_modkit(&[
        "bedmethyl",
        "merge",
        MIN_SAMPLES_BED_FP,
        MIN_SAMPLES_BED_FP,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_keep.to_str().unwrap(),
        "--force",
        "--min-samples",
        "2",
    ])
    .unwrap();
    assert_eq!(count_bed_records(&out_keep), n_positions);

    // --min-samples all resolves to the number of inputs (2 here), so it is
    // equivalent to --min-samples 2 and keeps every shared position.
    let out_all = std::env::temp_dir().join("test_merge_min_samples_all.bed");
    run_modkit(&[
        "bedmethyl",
        "merge",
        MIN_SAMPLES_BED_FP,
        MIN_SAMPLES_BED_FP,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_all.to_str().unwrap(),
        "--force",
        "--min-samples",
        "all",
    ])
    .unwrap();
    assert_eq!(count_bed_records(&out_all), n_positions);

    // --min-samples greater than the number of inputs is an error.
    let out_err = std::env::temp_dir().join("test_merge_min_samples_err.bed");
    assert!(run_modkit(&[
        "bedmethyl",
        "merge",
        MIN_SAMPLES_BED_FP,
        MIN_SAMPLES_BED_FP,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_err.to_str().unwrap(),
        "--force",
        "--min-samples",
        "3",
    ])
    .is_err());
}

#[test]
fn test_bedmethyl_merge_min_sample_coverage() {
    let sizes_fp = min_samples_test_sizes();

    // A very high per-sample valid-coverage floor excludes every record, so no
    // input contributes to any position and the output is empty.
    let out_bed =
        std::env::temp_dir().join("test_merge_min_sample_coverage.bed");
    run_modkit(&[
        "bedmethyl",
        "merge",
        MIN_SAMPLES_BED_FP,
        MIN_SAMPLES_BED_FP,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_bed.to_str().unwrap(),
        "--force",
        "--min-sample-coverage",
        "1000000000",
    ])
    .unwrap();
    assert_eq!(count_bed_records(&out_bed), 0);
}
