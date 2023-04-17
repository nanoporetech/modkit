use rust_htslib::{bam, bam::Read};

use crate::common::run_simple_summary;
use common::run_modkit;
use mod_kit::mod_bam::parse_raw_mod_tags;
use mod_kit::mod_base_code::{DnaBase, ModCode};

mod common;

#[test]
fn test_help() {
    let pileup_help_args = ["adjust-mods", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

fn tests_adjust_output(
    input_path: &str,
    output_path: &str,
    check_file_path: &str,
) -> Result<(), (bam::Record, bam::Record)> {
    let temp_file = std::env::temp_dir().join(output_path);
    let args = [
        "adjust-mods",
        "--ignore",
        "h",
        input_path,
        temp_file.to_str().unwrap(),
    ];
    run_modkit(&args).unwrap();
    assert!(temp_file.exists());

    let mut test_bam = bam::Reader::from_path(temp_file).unwrap();
    let mut ref_bam = bam::Reader::from_path(check_file_path).unwrap();
    for (test_res, ref_res) in test_bam.records().zip(ref_bam.records()) {
        let test_record = test_res.unwrap();
        let ref_record = ref_res.unwrap();
        if test_record == ref_record {
            continue;
        } else {
            return Err((test_record, ref_record));
        }
    }
    Ok(())
}

#[test]
fn test_adjust_canonical() {
    tests_adjust_output(
        "tests/resources/input_C.bam",
        "test_C.bam",
        "tests/resources/ref_out_C_auto.bam",
    )
    .unwrap();
}

#[test]
fn test_adjust_methyl() {
    tests_adjust_output(
        "tests/resources/input_5mC.bam",
        "test_5mC.bam",
        "tests/resources/ref_out_5mC_auto.bam",
    )
    .unwrap();
}

#[test]
fn test_adjust_no_tags() {
    let temp_file = std::env::temp_dir().join("test_out_no_tags.bam");
    run_modkit(&[
        "adjust-mods",
        "tests/resources/input_C_no_tags.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();
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

    run_modkit(&args).unwrap();
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

    let initial_mod_summary = run_simple_summary(
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        25,
    )
    .unwrap();

    let collapse_args = [
        "adjust-mods",
        "--convert",
        "h",
        "m",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_convered_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();

    let converted_mod_summary =
        run_simple_summary(test_convered_bam.to_str().unwrap(), 25).unwrap();

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

    let initial_mod_summary = run_simple_summary(
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        25,
    )
    .unwrap();

    let collapse_args = [
        "adjust-mods",
        "--convert",
        "h",
        "C",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_convered_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();

    let converted_mod_summary =
        run_simple_summary(test_convered_bam.to_str().unwrap(), 25).unwrap();

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
    let test_converted_bam =
        std::env::temp_dir().join("test_convert_sum_probs_rename.bam");

    let initial_mod_summary = run_simple_summary(
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        25,
    )
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
        test_converted_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();

    let converted_mod_summary =
        run_simple_summary(test_converted_bam.to_str().unwrap(), 25).unwrap();

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
    run_modkit(&first_adjust_args).unwrap();
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
    run_modkit(&second_adjust_args).unwrap();
    let mut reader =
        bam::Reader::from_path(test_both_bam.to_str().unwrap()).unwrap();
    for record in reader.records().map(|r| r.expect("should parse record")) {
        let raw_mod_tags = parse_raw_mod_tags(&record).unwrap().unwrap();
        let mm = raw_mod_tags.get_raw_mm();
        assert!(mm.starts_with("C+C?"));
    }
}

#[test]
fn test_adjust_out_of_spec_codes() {
    let adjusted_bam =
        std::env::temp_dir().join("test_adjust_out_of_spec_codes.bam");
    run_modkit(&[
        "adjust-mods",
        "tests/resources/bc_anchored_10_reads_old_tags.bam",
        adjusted_bam.to_str().unwrap(),
        "--convert",
        "Z",
        "m",
        "--convert",
        "Y",
        "h",
    ])
    .expect("should run adjust");

    let expected_summary = summarize_modbam(
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        1,
        0f32,
        None,
    )
    .expect("should get expected summary");
    let adjusted_summary =
        summarize_modbam(adjusted_bam.to_str().unwrap(), 1, 0f32, None)
            .expect("should get adjusted summary");
    assert_eq!(expected_summary, adjusted_summary);
}
