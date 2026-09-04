//! Tests with PacBio Jasmine-style tags (`A+a.;C+h?;C+m?;T-a.`), see
//! tests/make_pacbio_style_tags.py. The 5mC and 5hmC probabilities at the
//! first C of each record sum to more than one, the second C carries a
//! regular 5hmC call, one record has no MM/ML tags, 6mA is called on both
//! strands and the records are haplotagged (HP=1, HP=2 or untagged).
use std::{collections::HashMap, path::PathBuf};

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

/// bedMethyl rows split into columns: [3] mod code, [5] strand,
/// [9] N_valid_cov, [11] N_mod, [12] N_canonical, [13] N_other_mod.
fn parse_rows(lines: &[String]) -> Vec<Vec<String>> {
    lines
        .iter()
        .map(|l| l.split('\t').map(|x| x.to_string()).collect())
        .collect()
}

type RowKey = (String, String, String);

fn row_key(row: &[String]) -> RowKey {
    (row[0].clone(), row[1].clone(), row[5].clone())
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
    // regression: identical to the output of the general workers before the
    // --phased and --modified-bases changes (snapshot made with 697de7b)
    assert_eq!(
        general_lines,
        read_lines(&PathBuf::from(
            "../tests/resources/pacbio_style_tags_pileup_general.bed"
        ))
    );

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

#[test]
fn test_pileup_pacbio_style_modified_bases_filters_rows() {
    // The general workers (opposite-strand 6mA calls). Without
    // --modified-bases every code in the BAM gets a row, with
    // --modified-bases 5mC only the m rows are written and the h calls are
    // counted in N_other of the m rows.
    let run = |extra: &[&str], out_fn: &str| -> Vec<Vec<String>> {
        let out = std::env::temp_dir().join(out_fn);
        let mut args = vec![
            "pileup",
            INPUT,
            out.to_str().unwrap(),
            "--motif",
            "CG",
            "0",
            "--ref",
            REFERENCE,
            "--no-filtering",
        ];
        args.extend_from_slice(extra);
        run_modkit(&args).unwrap();
        parse_rows(&read_lines(&out))
    };
    let all_rows = run(&[], "test_pacbio_style_pileup_all_codes.bed");
    let m_rows = run(
        &["--modified-bases", "5mC"],
        "test_pacbio_style_pileup_5mc_only.bed",
    );
    assert!(all_rows.iter().any(|r| r[3] == "h"));
    assert!(all_rows.iter().any(|r| r[3] == "m"));
    assert!(!m_rows.is_empty());
    assert!(m_rows.iter().all(|r| r[3] == "m"));
    // the m rows themselves are not changed by the filter
    let expected_m_rows = all_rows
        .iter()
        .filter(|r| r[3] == "m")
        .cloned()
        .collect::<Vec<Vec<String>>>();
    assert_eq!(m_rows, expected_m_rows);
    // N_other of every m row is N_mod of the h row at the same position
    let h_n_mod = all_rows
        .iter()
        .filter(|r| r[3] == "h")
        .map(|r| (row_key(r), r[11].parse::<usize>().unwrap()))
        .collect::<HashMap<RowKey, usize>>();
    let mut n_other_total = 0usize;
    for row in &m_rows {
        let n_other = row[13].parse::<usize>().unwrap();
        assert_eq!(Some(&n_other), h_n_mod.get(&row_key(row)), "{row:?}");
        n_other_total += n_other;
    }
    // the fixture has a regular 5hmC call on the second C of each record
    assert!(n_other_total > 0);

    // --combine-mods sums all codes into one row per position instead. (With
    // --combine-mods a single-base `C` motif is added next to `CG`, so the
    // rows carry multiple-motif labels; only the CG rows are compared.)
    let combined_rows = run(
        &["--modified-bases", "5mC", "--combine-mods"],
        "test_pacbio_style_pileup_combine_mods.bed",
    )
    .into_iter()
    .filter(|r| r[3] == "C,CG,0")
    .collect::<Vec<Vec<String>>>();
    assert_eq!(combined_rows.len(), m_rows.len());
    for (combined, m) in combined_rows.iter().zip(m_rows.iter()) {
        assert_eq!(row_key(combined), row_key(m));
        let n_mod_combined = combined[11].parse::<usize>().unwrap();
        let n_mod = m[11].parse::<usize>().unwrap();
        let n_other = m[13].parse::<usize>().unwrap();
        assert_eq!(n_mod_combined, n_mod + n_other);
    }
}

#[test]
fn test_pileup_pacbio_style_phased() {
    // The general workers (opposite-strand 6mA calls) with --phased: the
    // counts are partitioned on the HP tag, untagged records only count
    // towards the combined output.
    let out_dir = std::env::temp_dir().join("test_pacbio_style_pileup_phased");
    run_modkit(&[
        "pileup",
        INPUT,
        out_dir.to_str().unwrap(),
        "--cpg",
        "--ref",
        REFERENCE,
        "--modified-bases",
        "5mC",
        "--no-filtering",
        "--phased",
        "--prefix",
        "pb",
    ])
    .unwrap();
    let combined =
        parse_rows(&read_lines(&out_dir.join("pb_combined.bedmethyl")));
    let hp1 = parse_rows(&read_lines(&out_dir.join("pb_hp1.bedmethyl")));
    let hp2 = parse_rows(&read_lines(&out_dir.join("pb_hp2.bedmethyl")));
    assert!(!hp1.is_empty());
    assert!(!hp2.is_empty());
    assert_ne!(hp1, hp2);
    for rows in [&combined, &hp1, &hp2] {
        assert!(rows.iter().all(|r| r[3] == "m"));
    }

    // the combined output is the same as without --phased
    let out_unphased =
        std::env::temp_dir().join("test_pacbio_style_pileup_unphased.bed");
    run_modkit(&[
        "pileup",
        INPUT,
        out_unphased.to_str().unwrap(),
        "--cpg",
        "--ref",
        REFERENCE,
        "--modified-bases",
        "5mC",
        "--no-filtering",
    ])
    .unwrap();
    assert_eq!(combined, parse_rows(&read_lines(&out_unphased)));

    // per position the haplotype coverages sum to at most the combined
    // coverage, and in total to less because of the untagged records
    let coverage = |rows: &[Vec<String>]| -> HashMap<RowKey, usize> {
        rows.iter()
            .map(|r| (row_key(r), r[9].parse::<usize>().unwrap()))
            .collect()
    };
    let combined_cov = coverage(&combined);
    let hp1_cov = coverage(&hp1);
    let hp2_cov = coverage(&hp2);
    for key in hp1_cov.keys().chain(hp2_cov.keys()) {
        assert!(combined_cov.contains_key(key), "{key:?}");
    }
    let mut total_haplotagged = 0usize;
    let mut total_combined = 0usize;
    for (key, cov) in &combined_cov {
        let haplotagged = hp1_cov.get(key).copied().unwrap_or(0)
            + hp2_cov.get(key).copied().unwrap_or(0);
        assert!(haplotagged <= *cov, "{key:?}");
        total_haplotagged += haplotagged;
        total_combined += cov;
    }
    assert!(total_haplotagged > 0);
    assert!(total_haplotagged < total_combined);
}
