use crate::common::{check_against_expected_text_file, run_modkit};
use anyhow::{anyhow, Context};
use mod_kit::mod_bam::{parse_raw_mod_tags, ModBaseInfo};
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

mod common;

pub fn check_two_bams_mod_probs_are_the_same(
    observed: &str,
    expected: &str,
) -> anyhow::Result<()> {
    let mut test_bam = bam::Reader::from_path(observed).unwrap();
    let mut ref_bam = bam::Reader::from_path(expected).unwrap();
    for (test_res, ref_res) in test_bam.records().zip(ref_bam.records()) {
        let test_record = test_res.unwrap();
        let ref_record = ref_res.unwrap();
        let test_modbase_info =
            ModBaseInfo::new_from_record(&test_record).unwrap();
        let ref_modbase_info =
            ModBaseInfo::new_from_record(&ref_record).unwrap();

        let ref_mod_probs = ref_modbase_info
            .into_iter_base_mod_probs()
            .1
            .map(|(base, strand, probs)| ((base, strand), probs))
            .collect::<HashMap<_, _>>();
        for (base, strand, seq_pos_probs) in
            test_modbase_info.iter_seq_base_mod_probs()
        {
            let ref_probs = ref_mod_probs.get(&(*base, strand)).unwrap();
            if ref_probs != seq_pos_probs {
                let test_record_id =
                    String::from_utf8(test_record.qname().to_vec()).unwrap();
                let expected_record_id =
                    String::from_utf8(ref_record.qname().to_vec()).unwrap();
                return Err(anyhow!(
                    "difference at test record id {test_record_id} =/= \
                     {expected_record_id}"
                ));
            }
        }

        if test_record == ref_record {
            continue;
        } else {
        }
    }
    Ok(())
}

#[test]
fn test_call_mods_basic_regression() {
    // Tests BAM against one checked by eye. Canary test, there has been a
    // change in the algorithm if this tests fails, but not necessarily
    // because it's broken
    let mod_call_out_bam =
        std::env::temp_dir().join("test_call_mods_same_positions_mod_call.bam");
    run_modkit(&[
        "call-mods",
        "tests/resources/ecoli_reg.sorted.bam",
        mod_call_out_bam.to_str().unwrap(),
        "--filter-threshold",
        "A:0.65",
        "--mod-threshold",
        "a:0.95",
        "--filter-threshold",
        "C:0.85",
        "--mod-threshold",
        "m:0.95",
    ])
    .unwrap();
    check_two_bams_mod_probs_are_the_same(
        mod_call_out_bam.to_str().unwrap(),
        "tests/resources/ecoli_reg.call_mods.bam",
    )
    .unwrap();
}

#[test]
fn test_call_mods_keeps_all_mod_calls() {
    let n_lines_in_file = |fp: &PathBuf| -> anyhow::Result<usize> {
        let reader = BufReader::new(File::open(fp)?);
        Ok(reader.lines().map(|_| 1).sum::<usize>())
    };
    let extract_control_fp = std::env::temp_dir()
        .join("test_call_mods_keeps_all_mod_calls_control.tsv");
    let extract_call_mods_fp = std::env::temp_dir()
        .join("test_call_mods_keeps_all_mod_calls_call_mods.bam");
    let extract_call_mods_fp_tsv = std::env::temp_dir()
        .join("test_call_mods_keeps_all_mod_calls_call_mods.tsv");

    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        extract_control_fp.to_str().unwrap(),
        "--force",
    ])
    .unwrap();

    run_modkit(&[
        "call-mods",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        extract_call_mods_fp.to_str().unwrap(),
        "--no-filtering",
    ])
    .unwrap();
    run_modkit(&[
        "extract",
        extract_call_mods_fp.to_str().unwrap(),
        extract_call_mods_fp_tsv.to_str().unwrap(),
        "--force",
    ])
    .unwrap();
    let before_call_mods = n_lines_in_file(&extract_control_fp).unwrap();
    let after_call_mods = n_lines_in_file(&extract_call_mods_fp_tsv).unwrap();
    assert_eq!(before_call_mods, after_call_mods);
}

#[test]
fn test_call_mods_same_pileup() {
    // tests that the pileup generated from a pre-thresholded BAM is equivalent
    // to one where the thresholds are used _during_ pileup
    let mod_call_out_bam =
        std::env::temp_dir().join("test_call_mods_same_pileup.bam");
    run_modkit(&[
        "call-mods",
        "tests/resources/ecoli_reg.sorted.bam",
        mod_call_out_bam.to_str().unwrap(),
        "--filter-threshold",
        "A:0.65",
        "--mod-threshold",
        "a:0.95",
        "--filter-threshold",
        "C:0.85",
        "--mod-threshold",
        "m:0.95",
    ])
    .unwrap();
    bam::index::build(mod_call_out_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();
    let mod_called_pileup =
        std::env::temp_dir().join("test_call_mods_same_pileup-1.bed");
    run_modkit(&[
        "pileup",
        mod_call_out_bam.to_str().unwrap(),
        mod_called_pileup.to_str().unwrap(),
        "--no-filtering",
    ])
    .unwrap();
    let in_situ_threshold_pileup =
        std::env::temp_dir().join("test_call_mods_same_pileup-2.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/ecoli_reg.sorted.bam",
        in_situ_threshold_pileup.to_str().unwrap(),
        "--filter-threshold",
        "A:0.65",
        "--mod-threshold",
        "a:0.95",
        "--filter-threshold",
        "C:0.85",
        "--mod-threshold",
        "m:0.95",
    ])
    .unwrap();
    check_against_expected_text_file(
        mod_called_pileup.to_str().unwrap(),
        in_situ_threshold_pileup.to_str().unwrap(),
    );
}

#[test]
fn test_call_mods_supplementary_secondary() {
    fn check(bam_fp: &PathBuf) {
        let mut reader = bam::Reader::from_path(&bam_fp).unwrap();
        let n_records = reader
            .records()
            .map(|r| r.unwrap())
            .map(|record| parse_raw_mod_tags(&record).unwrap())
            .count();
        assert_eq!(n_records, 3);
    }

    let out_bam =
        std::env::temp_dir().join("test_call_mods_supplementary_secondary.bam");
    run_modkit(&[
        "adjust-mods",
        "tests/resources/test_supplementary_secondary.bam",
        out_bam.to_str().unwrap(),
        "--ignore",
        "h",
        "--ff",
    ])
    .context(format!("failed to run adjust"))
    .unwrap();
    check(&out_bam);
}
