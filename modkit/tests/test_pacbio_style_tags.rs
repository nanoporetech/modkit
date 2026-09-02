//! Tests with PacBio Jasmine-style tags (`A+a.;C+h?;C+m?;T-a.`), see
//! tests/make_pacbio_style_tags.py. The 5mC and 5hmC probabilities at the
//! first C of each record sum to more than one, one record has no MM/ML tags
//! and 6mA is called on both strands.
use std::path::PathBuf;

use common::run_modkit;
use rust_htslib::{bam, bam::Read};

mod common;

const INPUT: &str = "../tests/resources/pacbio_style_tags.bam";
const INPUT_NO_UNTAGGED: &str =
    "../tests/resources/pacbio_style_tags_no_untagged.bam";
const REFERENCE: &str = "../tests/resources/CGI_ladder_3.6kb_ref.fa";

fn count_records(fp: &PathBuf) -> usize {
    let mut reader = bam::Reader::from_path(fp).unwrap();
    reader.records().count()
}

fn read_lines(fp: &PathBuf) -> Vec<String> {
    std::fs::read_to_string(fp)
        .unwrap()
        .lines()
        .map(|l| l.to_string())
        .collect()
}

#[test]
fn test_adjust_mods_keeps_records_with_conflicting_probs() {
    let out_bam = std::env::temp_dir().join("test_pacbio_style_adjust.bam");
    run_modkit(&[
        "adjust-mods",
        "--ignore",
        "h",
        INPUT,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    // 10 records in, the record without MM/ML tags is dropped, the 9 records
    // with a position where P(5mC) + P(5hmC) > 1 are all kept.
    assert_eq!(count_records(&out_bam), 9);
    // the output has no conflicts left, check-tags exits 0
    run_modkit(&[
        "modbam",
        "check-tags",
        out_bam.to_str().unwrap(),
        "--interval-size",
        "20",
    ])
    .unwrap();

    // converting h to m also keeps every record
    let out_bam = std::env::temp_dir().join("test_pacbio_style_convert.bam");
    run_modkit(&[
        "adjust-mods",
        "--convert",
        "h",
        "m",
        INPUT,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    assert_eq!(count_records(&out_bam), 9);
}

#[test]
fn test_check_tags_reports_conflicting_probs() {
    // check-tags stays strict, the conflicting records make it exit 1
    assert!(run_modkit(&[
        "modbam",
        "check-tags",
        INPUT,
        "--interval-size",
        "20"
    ])
    .is_err());

    let tmp_dir = std::env::temp_dir().join("test_pacbio_style_check_tags");
    run_modkit(&[
        "modbam",
        "check-tags",
        INPUT,
        "--interval-size",
        "20",
        "--permissive",
        "--force",
        "--out-dir",
        tmp_dir.to_str().unwrap(),
    ])
    .unwrap();
    let error_counts = read_lines(&tmp_dir.join("error_counts.tsv"));
    let get_count = |error: &str| -> usize {
        error_counts
            .iter()
            .find_map(|line| {
                let parts = line.split('\t').collect::<Vec<&str>>();
                if parts[0] == error {
                    Some(parts[1].parse::<usize>().unwrap())
                } else {
                    None
                }
            })
            .unwrap_or_else(|| panic!("{error} not in {error_counts:?}"))
    };
    assert_eq!(get_count("conflict-explicit-prob-greater-than-one"), 9);
    assert_eq!(get_count("MM-tag-missing"), 1);
}

#[test]
fn test_pileup_pacbio_style_tags() {
    // The general workers, a record without MM/ML tags must not abort the
    // interval.
    let out_general =
        std::env::temp_dir().join("test_pacbio_style_pileup_general.bed");
    run_modkit(&[
        "pileup",
        INPUT,
        out_general.to_str().unwrap(),
        "--motif",
        "CG",
        "0",
        "--ref",
        REFERENCE,
        "--no-filtering",
    ])
    .unwrap();
    let general_lines = read_lines(&out_general);
    assert!(!general_lines.is_empty());

    // Same input without the record lacking MM/ML tags gives the same output.
    let out_no_untagged =
        std::env::temp_dir().join("test_pacbio_style_pileup_no_untagged.bed");
    run_modkit(&[
        "pileup",
        INPUT_NO_UNTAGGED,
        out_no_untagged.to_str().unwrap(),
        "--motif",
        "CG",
        "0",
        "--ref",
        REFERENCE,
        "--no-filtering",
    ])
    .unwrap();
    assert_eq!(general_lines, read_lines(&out_no_untagged));

    // Requesting a preset (--cpg --modified-bases) would use the optimized
    // workers, which do not support the opposite-strand `T-a.` calls; pileup
    // must detect this and fall back to the general workers. With
    // --modified-bases 5mC the other modifications are counted as N_other, so
    // only the 5mC and canonical counts are compared with the run above.
    let run_preset = |input: &str, out_fn: &str| -> Vec<String> {
        let out_preset = std::env::temp_dir().join(out_fn);
        run_modkit(&[
            "pileup",
            input,
            out_preset.to_str().unwrap(),
            "--cpg",
            "--ref",
            REFERENCE,
            "--modified-bases",
            "5mC",
            "--no-filtering",
        ])
        .unwrap();
        read_lines(&out_preset)
    };
    let preset_lines = run_preset(INPUT, "test_pacbio_style_pileup_preset.bed");
    assert!(!preset_lines.is_empty());
    let m_counts = |lines: &[String]| -> Vec<Vec<String>> {
        lines
            .iter()
            .map(|l| l.split('\t').map(|x| x.to_string()).collect::<Vec<_>>())
            .filter(|parts| parts[3] == "m")
            // chrom, start, strand, N_mod, N_canonical
            .map(|parts| {
                vec![
                    parts[0].clone(),
                    parts[1].clone(),
                    parts[5].clone(),
                    parts[11].clone(),
                    parts[12].clone(),
                ]
            })
            .collect()
    };
    assert!(!m_counts(&preset_lines).is_empty());
    assert_eq!(m_counts(&general_lines), m_counts(&preset_lines));
    // and the record without MM/ML tags does not change the output
    assert_eq!(
        preset_lines,
        run_preset(
            INPUT_NO_UNTAGGED,
            "test_pacbio_style_pileup_preset_no_untagged.bed"
        )
    );
}

#[test]
fn test_extract_pacbio_style_tags() {
    let out_tsv = std::env::temp_dir().join("test_pacbio_style_extract.tsv");
    run_modkit(&[
        "extract",
        "full",
        INPUT,
        out_tsv.to_str().unwrap(),
        "--force",
    ])
    .unwrap();
    let lines = read_lines(&out_tsv);
    // header plus rows from all 9 tagged records
    let mut read_ids = lines
        .iter()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect::<Vec<String>>();
    read_ids.sort();
    read_ids.dedup();
    assert_eq!(read_ids.len(), 9);
}
