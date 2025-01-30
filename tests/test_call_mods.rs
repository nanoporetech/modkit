use crate::common::run_modkit;
use anyhow::{anyhow, Context};
use mod_kit::dmr::bedmethyl::BedMethylLine;
use mod_kit::errs::MkError;
use mod_kit::mod_bam::{
    BaseModCall, ModBaseInfo, RawModTags, SeqPosBaseModProbs,
};
use mod_kit::mod_base_code::{
    DnaBase, ModCodeRepr, ParseChar, METHYL_CYTOSINE, SIX_METHYL_ADENINE,
};
use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
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
            let ref_probs = ref_mod_probs.get(&(base, strand)).unwrap();
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

fn get_mod_probs(fp: &str) -> HashMap<String, ModBaseInfo> {
    let mut reader = bam::Reader::from_path(fp).unwrap();
    reader
        .records()
        .map(|r| r.unwrap())
        .filter_map(|rec| match ModBaseInfo::new_from_record(&rec) {
            Err(MkError::NonPrimaryMissingMn) => None,
            Err(e) => panic!("should extract modbase info, {e}"),
            Ok(modbase_info) => Some((
                rec.qname().iter().map(|b| *b as char).collect::<String>(),
                modbase_info,
            )),
        })
        .collect()
}

#[test]
fn test_call_mods_thresholds_correctly() {
    // Tests BAM against one checked by eye. Canary test, there has been a
    // change in the algorithm if this tests fails, but not necessarily
    // because it's broken
    let mod_call_out_bam =
        std::env::temp_dir().join("test_call_mods_same_positions_mod_call.bam");
    let test_bam_fp = "tests/resources/ecoli_reg.sorted.bam";
    let uncalled_mod_probs = get_mod_probs(test_bam_fp);
    run_modkit(&[
        "call-mods",
        &test_bam_fp,
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
    let called_mod_probs = get_mod_probs(mod_call_out_bam.to_str().unwrap());
    let canonical_thresholds =
        HashMap::from([(DnaBase::C, 0.85), (DnaBase::A, 0.65)]);
    let per_mod_thresholds =
        HashMap::from([(SIX_METHYL_ADENINE, 0.95), (METHYL_CYTOSINE, 0.95)]);

    let caller = MultipleThresholdModCaller::new(
        canonical_thresholds.clone(),
        per_mod_thresholds.clone(),
        0.0,
    );

    fn check_probs(
        called: &SeqPosBaseModProbs,
        uncalled: &SeqPosBaseModProbs,
        caller: &MultipleThresholdModCaller,
        primary_base: &DnaBase,
        canonical_threshold: f32,
        mod_thresholds: &HashMap<ModCodeRepr, f32>,
    ) -> bool {
        uncalled.pos_to_base_mod_probs.iter().all(
            |(position, base_mod_probs)| {
                let base_mod_call = called
                    .pos_to_base_mod_probs
                    .get(position)
                    .map(|x| x.argmax_base_mod_call());
                match caller.call(primary_base, base_mod_probs) {
                    BaseModCall::Canonical(p) => {
                        p >= canonical_threshold
                            && base_mod_call.unwrap().is_canonical()
                    }
                    BaseModCall::Modified(p, code) => {
                        let mod_threshold = mod_thresholds.get(&code).unwrap();
                        p >= *mod_threshold
                            && base_mod_call.unwrap().is_match_modcall(&code)
                    }
                    BaseModCall::Filtered => base_mod_call.is_none(),
                }
            },
        )
    }

    for (read_id, mod_base_info) in called_mod_probs {
        let uncalled_info = uncalled_mod_probs.get(&read_id).unwrap();

        assert!(mod_base_info.neg_seq_base_mod_probs.is_empty());
        assert!(uncalled_info.neg_seq_base_mod_probs.is_empty());

        for (dna_base, uncalled_probs) in
            uncalled_info.pos_seq_base_mod_probs.iter()
        {
            let called_probs =
                mod_base_info.pos_seq_base_mod_probs.get(&dna_base).unwrap();
            let can_thresh = *canonical_thresholds.get(&dna_base).unwrap();
            assert!(check_probs(
                &called_probs,
                &uncalled_probs,
                &caller,
                &dna_base,
                can_thresh,
                &per_mod_thresholds
            ))
        }
    }
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
        "full",
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
        "full",
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
    let update_tags_bam = std::env::temp_dir()
        .join("test_call_mods_same_pileup_updated_tags.bam");
    run_modkit(&[
        "update-tags",
        "tests/resources/ecoli_reg.sorted.bam",
        update_tags_bam.to_str().unwrap(),
        "--no-implicit-probs",
        "--mode",
        "explicit",
    ])
    .unwrap();
    bam::index::build(update_tags_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    let mod_call_out_bam =
        std::env::temp_dir().join("test_call_mods_same_pileup_called_mods.bam");
    run_modkit(&[
        "call-mods",
        update_tags_bam.to_str().unwrap(),
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
        update_tags_bam.to_str().unwrap(),
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

    let parse_bedmethyl_fp = |fp: &PathBuf| -> Vec<BedMethylLine> {
        let reader = BufReader::new(File::open(fp).unwrap());
        reader
            .lines()
            .map(|l| BedMethylLine::parse(l.unwrap().as_str()).unwrap())
            .collect()
    };

    let called_records = parse_bedmethyl_fp(&mod_called_pileup);
    let in_situ_records = parse_bedmethyl_fp(&in_situ_threshold_pileup);
    assert_eq!(called_records.len(), in_situ_records.len());
    for (x, y) in called_records.into_iter().zip(in_situ_records) {
        assert_eq!(x.chrom, y.chrom);
        assert_eq!(x.start(), y.start());
        assert_eq!(x.raw_mod_code, y.raw_mod_code);
        assert_eq!(x.strand, y.strand);
        assert_eq!(x.count_methylated, y.count_methylated);
        assert_eq!(x.valid_coverage, y.valid_coverage);
        assert_eq!(x.count_canonical, y.count_canonical);
        assert_eq!(x.count_other, y.count_other);
        assert_eq!(
            x.count_diff + x.count_nocall,
            y.count_fail + y.count_diff + y.count_nocall,
            "{x:?}\n{y:?}"
        )
    }
}

#[test]
fn test_call_mods_supplementary_secondary() {
    fn check(bam_fp: &PathBuf) {
        let mut reader = bam::Reader::from_path(&bam_fp).unwrap();
        let n_records = reader
            .records()
            .map(|r| r.unwrap())
            .map(|record| RawModTags::new_from_record(&record).unwrap())
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
