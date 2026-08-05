use std::{
    ffi::CString,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
    process::{Command, Output},
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
    let expected_output = input_records
        .iter()
        .cloned()
        .map(|mut line| {
            line.count_methylated *= 2;
            line.valid_coverage *= 2;
            line.count_canonical *= 2;
            line.count_other *= 2;
            line.count_delete *= 2;
            line.count_fail *= 2;
            line.count_diff *= 2;
            line.count_nocall *= 2;
            line.to_string()
        })
        .collect::<String>();
    assert_eq!(
        fs::read(&out_bed).unwrap(),
        expected_output.as_bytes(),
        "valid merge output must remain byte-identical"
    );
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

const PREFLIGHT_OUTPUT_SENTINEL: &[u8] = b"existing output must survive\n";

fn copy_bedmethyl_fixture(data_dst: &Path) {
    fs::copy(MIN_SAMPLES_BED_FP, data_dst).unwrap();
}

fn run_merge_for_preflight(
    inputs: &[&Path],
    sizes_fp: &Path,
    out_fp: &Path,
    extra_args: &[&str],
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args(["bedmethyl", "merge"]);
    for input in inputs {
        command.arg(input);
    }
    command
        .arg("-g")
        .arg(sizes_fp)
        .arg("-o")
        .arg(out_fp)
        .arg("--force")
        .args(extra_args)
        .output()
        .unwrap()
}

fn assert_preflight_failure(output: Output, bad_input: &Path, out_fp: &Path) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !output.status.success(),
        "merge unexpectedly accepted invalid input {}\nstderr: {stderr}",
        bad_input.display()
    );
    let bad_name = bad_input.file_name().unwrap().to_string_lossy();
    assert!(
        stderr.contains(bad_name.as_ref()),
        "error must identify invalid input {}\nstderr: {stderr}",
        bad_input.display()
    );
    assert_eq!(
        fs::read(out_fp).unwrap(),
        PREFLIGHT_OUTPUT_SENTINEL,
        "input preflight failure must precede output truncation"
    );
}

fn preflight_case_paths(
    tempdir: &tempfile::TempDir,
) -> (std::path::PathBuf, std::path::PathBuf) {
    let sizes_fp = tempdir.path().join("sizes.tsv");
    fs::write(&sizes_fp, b"chr20\t64444167\n").unwrap();
    let out_fp = tempdir.path().join("merged.bed");
    fs::write(&out_fp, PREFLIGHT_OUTPUT_SENTINEL).unwrap();
    (sizes_fp, out_fp)
}

#[test]
fn test_bedmethyl_merge_rejects_missing_input_before_output() {
    let tempdir = tempfile::tempdir().unwrap();
    let (sizes_fp, out_fp) = preflight_case_paths(&tempdir);
    let missing = tempdir.path().join("missing.bed.gz");
    let output = run_merge_for_preflight(
        &[Path::new(MIN_SAMPLES_BED_FP), &missing],
        &sizes_fp,
        &out_fp,
        &[],
    );
    assert_preflight_failure(output, &missing, &out_fp);
}

#[test]
fn test_bedmethyl_merge_rejects_absent_index_before_output() {
    let tempdir = tempfile::tempdir().unwrap();
    let (sizes_fp, out_fp) = preflight_case_paths(&tempdir);
    let unindexed = tempdir.path().join("unindexed.bed.gz");
    copy_bedmethyl_fixture(&unindexed);
    let output = run_merge_for_preflight(
        &[Path::new(MIN_SAMPLES_BED_FP), &unindexed],
        &sizes_fp,
        &out_fp,
        &[],
    );
    assert_preflight_failure(output, &unindexed, &out_fp);
}

#[test]
fn test_bedmethyl_merge_rejects_corrupt_index_before_output() {
    let tempdir = tempfile::tempdir().unwrap();
    let (sizes_fp, out_fp) = preflight_case_paths(&tempdir);
    let corrupt = tempdir.path().join("corrupt-index.bed.gz");
    copy_bedmethyl_fixture(&corrupt);
    fs::write(format!("{}.tbi", corrupt.display()), b"not a tabix index")
        .unwrap();
    let output = run_merge_for_preflight(
        &[Path::new(MIN_SAMPLES_BED_FP), &corrupt],
        &sizes_fp,
        &out_fp,
        &[],
    );
    assert_preflight_failure(output, &corrupt, &out_fp);
}

#[test]
fn test_bedmethyl_merge_rejects_one_bad_input_with_min_samples_all() {
    let tempdir = tempfile::tempdir().unwrap();
    let (sizes_fp, out_fp) = preflight_case_paths(&tempdir);
    let missing = tempdir.path().join("missing-middle.bed.gz");
    let valid = Path::new(MIN_SAMPLES_BED_FP);
    let output = run_merge_for_preflight(
        &[valid, &missing, valid],
        &sizes_fp,
        &out_fp,
        &["--min-samples", "all"],
    );
    assert_preflight_failure(output, &missing, &out_fp);
}

#[test]
fn test_bedmethyl_merge_accepts_csi_index() {
    let tempdir = tempfile::tempdir().unwrap();
    let (sizes_fp, out_fp) = preflight_case_paths(&tempdir);
    let csi_input = tempdir.path().join("csi-indexed.bed.gz");
    copy_bedmethyl_fixture(&csi_input);

    let csi_input_cstr = CString::new(csi_input.to_str().unwrap()).unwrap();
    let result = unsafe {
        rust_htslib::htslib::tbx_index_build(
            csi_input_cstr.as_ptr(),
            14,
            &rust_htslib::htslib::tbx_conf_bed,
        )
    };
    assert_eq!(result, 0, "failed to build CSI fixture");
    assert!(Path::new(&format!("{}.csi", csi_input.display())).exists());

    let output = run_merge_for_preflight(
        &[&csi_input, &csi_input],
        &sizes_fp,
        &out_fp,
        &[],
    );
    assert!(
        output.status.success(),
        "merge rejected a valid CSI-backed input: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(count_bed_records(&out_fp) > 0);
}

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
