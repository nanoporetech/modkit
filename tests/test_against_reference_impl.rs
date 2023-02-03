use std::fs::File;
use std::io::{BufRead, BufReader, Read as StdRead};
use std::process::Output;

use mod_kit::mod_bam::parse_raw_mod_tags;
use mod_kit::mod_base_code::{DnaBase, ModCode};
use mod_kit::summarize::summarize_modbam;
use rust_htslib::bam;
use rust_htslib::bam::Read;

#[test]
fn test_help() {
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let help = std::process::Command::new(exe)
        .arg("adjust-mods")
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();
    assert!(help.status.success());
}

fn run_modkit(args: &[&str]) -> Output {
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let output = std::process::Command::new(exe)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();
    assert!(output.status.success(), "failed to run {:?}", args);
    output
}

fn test_adjust_output(
    input_path: &str,
    output_path: &str,
    check_file_path: &str,
) {
    let temp_file = std::env::temp_dir().join(output_path);
    let args = ["adjust-mods", input_path, temp_file.to_str().unwrap()];
    run_modkit(&args);
    assert!(temp_file.exists());

    let mut test_bam = bam::Reader::from_path(temp_file).unwrap();
    let mut ref_bam = bam::Reader::from_path(check_file_path).unwrap();
    for (test_res, ref_res) in test_bam.records().zip(ref_bam.records()) {
        let test_record = test_res.unwrap();
        let ref_record = ref_res.unwrap();
        assert_eq!(ref_record, test_record);
    }
}

#[test]
fn test_adjust_canonical() {
    test_adjust_output(
        "tests/resources/input_C.bam",
        "test_C.bam",
        "tests/resources/ref_out_C_auto.bam",
    );
}

#[test]
fn test_adjust_methyl() {
    test_adjust_output(
        "tests/resources/input_5mC.bam",
        "test_5mC.bam",
        "tests/resources/ref_out_5mC_auto.bam",
    );
}

#[test]
fn test_adjust_no_tags() {
    let temp_file = std::env::temp_dir().join("test_out_no_tags.bam");
    run_modkit(&[
        "adjust-mods",
        "tests/resources/input_C_no_tags.bam",
        temp_file.to_str().unwrap(),
    ]);
}

fn check_against_expected_text_file(output_fp: &str, expected_fp: &str) {
    let test = {
        let mut fh = std::fs::File::open(output_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };
    let expected = {
        // this file was hand-checked for correctness.
        let mut fh = std::fs::File::open(expected_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };

    similar_asserts::assert_eq!(test, expected);
}

#[test]
fn test_mod_pileup_no_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--output-bed",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args);

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_nofilt.bed",
    );
}

#[test]
fn test_mod_pileup_with_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_withfilt.bed");
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-f",
        "1.0",
        "-p",
        "0.25",
        "--seed",
        "42",
        "--output-bed",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args);

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_filt025.bed",
    );
}

#[test]
fn test_mod_pileup_combine() {
    let test_adjusted_bam = std::env::temp_dir().join("test_combined.bam");
    let pileup_args = [
        "pileup",
        "--combine",
        "--no-filtering",
        "--output-bed",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_adjusted_bam.to_str().unwrap(),
    ];
    run_modkit(&pileup_args);
    assert!(test_adjusted_bam.exists());

    check_against_expected_text_file(
        test_adjusted_bam.to_str().unwrap(),
        "tests/resources/modbam.modpileup_combined.bed",
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
    run_modkit(&collapse_args);
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
    run_modkit(&pileup_args);
    assert!(test_collapsed_bed.exists());

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--collapse",
        "h",
        "--method",
        "norm",
        "--no-filtering",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_restricted_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args);
    assert!(test_restricted_bed.exists());
    check_against_expected_text_file(
        test_restricted_bed.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    );
}
#[test]
fn test_adjust_to_no_mods() {
    let test_ignore_h_bam =
        std::env::temp_dir().join("test_adjust_to_no_mods_ignore_h.bam");
    let test_both_bam =
        std::env::temp_dir().join("test_adjust_to_no_mods_ignore_both.bam");
    let first_adjust_args = [
        "adjust-mods",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_ignore_h_bam.to_str().unwrap(),
    ];
    run_modkit(&first_adjust_args);
    let mut reader =
        bam::Reader::from_path(test_ignore_h_bam.to_str().unwrap()).unwrap();
    for record in reader.records().map(|r| r.expect("should parse record")) {
        let raw_mod_tags = parse_raw_mod_tags(&record).unwrap().unwrap();
        let mm = raw_mod_tags.get_raw_mm();
        assert!(mm.starts_with("C+m?"));
    }
    let second_adjust_args = [
        "adjust-mods",
        "--ignore",
        "m",
        test_ignore_h_bam.to_str().unwrap(),
        test_both_bam.to_str().unwrap(),
    ];
    run_modkit(&second_adjust_args);
    let mut reader =
        bam::Reader::from_path(test_both_bam.to_str().unwrap()).unwrap();
    for record in reader.records().map(|r| r.expect("should parse record")) {
        let raw_mod_tags = parse_raw_mod_tags(&record).unwrap().unwrap();
        let mm = raw_mod_tags.get_raw_mm();
        assert!(mm.starts_with("C+C?"));
    }
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

    run_modkit(&args);

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
    ]);
    assert!(updated_file.exists());
    bam::index::build(updated_file.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    let out_file = std::env::temp_dir().join("test_pileup_old_tags.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        updated_file.to_str().unwrap(),
        out_file.to_str().unwrap(),
    ]);
    check_against_expected_text_file(
        out_file.to_str().unwrap(),
        "tests/resources/pileup-old-tags-regressiontest.bed",
    );
}

#[test]
fn test_adjust_convert_old_tags() {
    let out_file =
        std::env::temp_dir().join("test_adjust_convert_old_tags.bam");
    let args = [
        "adjust-mods",
        "--convert",
        "m",
        "C",
        "tests/resources/HG002_small.ch20._other.sorted.bam",
        out_file.to_str().unwrap(),
    ];

    run_modkit(&args);
    let mut reader =
        bam::Reader::from_path(out_file.to_str().unwrap()).unwrap();
    for record in reader.records().map(|r| r.expect("should parse record")) {
        let raw_mod_tags = parse_raw_mod_tags(&record).unwrap().unwrap();
        assert!(!raw_mod_tags.mm_is_new_style());
        assert!(!raw_mod_tags.ml_is_new_style());
        let mm = raw_mod_tags.get_raw_mm();
        if !mm.is_empty() {
            assert!(mm.starts_with("C+C,"), "wrong: {mm}");
        }
    }
}

#[test]
fn test_mod_adjust_convert_sum_probs() {
    let test_convered_bam =
        std::env::temp_dir().join("test_convert_sum_probs.bam");

    let initial_mod_summary =
        summarize_modbam("tests/resources/bc_anchored_10_reads.sorted.bam", 1)
            .unwrap();

    let collapse_args = [
        "adjust-mods",
        "--convert",
        "h",
        "m",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_convered_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args);

    let converted_mod_summary =
        summarize_modbam(test_convered_bam.to_str().unwrap(), 1).unwrap();

    let initial_m_calls = initial_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::m))
        .unwrap();
    let initial_h_calls = initial_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::h))
        .unwrap();

    let converted_m_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::m))
        .unwrap();
    assert_eq!(*converted_m_calls, initial_m_calls + initial_h_calls);
    let converted_h_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::h));
    assert!(converted_h_calls.is_none());
}

#[test]
fn test_mod_adjust_convert_rename() {
    let test_convered_bam =
        std::env::temp_dir().join("test_convert_convert_rename.bam");

    let initial_mod_summary =
        summarize_modbam("tests/resources/bc_anchored_10_reads.sorted.bam", 1)
            .unwrap();

    let collapse_args = [
        "adjust-mods",
        "--convert",
        "h",
        "C",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_convered_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args);

    let converted_mod_summary =
        summarize_modbam(test_convered_bam.to_str().unwrap(), 1).unwrap();

    let initial_h_calls = initial_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::h))
        .unwrap();
    let converted_any_c_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::anyC))
        .unwrap();
    assert_eq!(initial_h_calls, converted_any_c_calls);
}

#[test]
fn test_mod_adjust_convert_sum_probs_rename() {
    let test_convered_bam =
        std::env::temp_dir().join("test_convert_sum_probs_rename.bam");

    let initial_mod_summary =
        summarize_modbam("tests/resources/bc_anchored_10_reads.sorted.bam", 1)
            .unwrap();

    let collapse_args = [
        "adjust-mods",
        "--convert",
        "h",
        "C",
        "--convert",
        "m",
        "C",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_convered_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args);

    let converted_mod_summary =
        summarize_modbam(test_convered_bam.to_str().unwrap(), 1).unwrap();

    let initial_m_calls = initial_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::m))
        .unwrap();
    let initial_h_calls = initial_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::h))
        .unwrap();

    let converted_any_c_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::anyC))
        .unwrap();
    assert_eq!(*converted_any_c_calls, initial_m_calls + initial_h_calls);
    let converted_h_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::h));
    assert!(converted_h_calls.is_none());
    let converted_m_calls = converted_mod_summary
        .mod_call_counts
        .get(&DnaBase::C)
        .and_then(|counts| counts.get(&ModCode::m));
    assert!(converted_m_calls.is_none());
}
