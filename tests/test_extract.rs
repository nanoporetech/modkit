use crate::common::{
    check_against_expected_text_file, parse_mod_profile, ModData,
};
use anyhow::{anyhow, Context};
use common::run_modkit;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::{Path, PathBuf};

mod common;

fn parse_bed_file(fp: &PathBuf) -> HashMap<String, HashSet<(i64, char)>> {
    let reader = BufReader::new(File::open(fp).unwrap());
    reader
        .lines()
        .map(|l| l.unwrap())
        .map(|line| {
            let parts = line.split_ascii_whitespace().collect::<Vec<&str>>();
            let contig = parts[0].to_owned();
            let pos = parts[1].parse::<i64>().unwrap();
            let strand = parts[5].parse::<char>().unwrap();
            (contig, (pos, strand))
        })
        .fold(HashMap::new(), |mut acc, (contig, (pos, strand))| {
            acc.entry(contig)
                .or_insert(HashSet::new())
                .insert((pos, strand));
            acc
        })
}

#[test]
fn test_mod_data_ord() {
    let mod_data1 =
        ModData::new(0, 1, 'm', '+', 100, "".to_string(), "".to_string());
    let mod_data2 =
        ModData::new(0, 1, 'h', '+', 100, "".to_string(), "".to_string());
    let mod_data3 =
        ModData::new(1, 1, 'h', '+', 100, "".to_string(), "".to_string());
    assert!(mod_data2 < mod_data1);
    assert!(mod_data1 < mod_data3);
}

fn check_mod_profiles_same(
    output_fp: &PathBuf,
    expected_fp: &PathBuf,
) -> anyhow::Result<()> {
    let output_profile = parse_mod_profile(output_fp).unwrap();
    let expected_profile = parse_mod_profile(expected_fp).unwrap();
    if output_profile == expected_profile {
        Ok(())
    } else {
        for (read, profile) in expected_profile.iter() {
            if let Some(obs) = output_profile.get(read) {
                if obs != profile {
                    return Err(anyhow!(
                        "read {read}'s profile is different expected, test fp: {:?}, expected fp {:?}",
                        &output_fp, &expected_fp
                    ));
                }
            } else {
                return Err(anyhow!("read {read} was missing from output"));
            }
        }
        Err(anyhow!("they were different"))
    }
}

#[test]
fn test_extract_correct_output() {
    let out_fp = std::env::temp_dir().join("test_extract_correct_output.tsv");
    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "-i",
        "25",
        "--force",
    ])
    .unwrap();
    check_mod_profiles_same(
        &out_fp,
        &Path::new(
            "tests/resources/bc_anchored_10_reads.sorted.methylprofile.tsv",
        )
        .to_path_buf(),
    )
    .context("test_extract_correct_output, failed")
    .unwrap();
}

#[test]
fn test_extract_correct_output_with_ref() {
    let out_fp =
        std::env::temp_dir().join("test_extract_correct_output_with_ref.tsv");

    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "-i",
        "25",
        "--force",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();

    check_mod_profiles_same(
        &out_fp,
        &Path::new(
            "tests/resources/bc_anchored_10_reads.sorted.methylprofile_ref.tsv",
        )
        .to_path_buf(),
    )
    .context("test_extract_correct_output_with_ref, failed")
    .unwrap();
}

#[test]
fn test_extract_duplex_correct_output() {
    let out_fp_sorted = std::env::temp_dir()
        .join("test_extract_duplex_correct_output_sorted.tsv");
    let out_fp =
        std::env::temp_dir().join("test_extract_duplex_correct_output.tsv");

    run_modkit(&[
        "extract",
        "tests/resources/duplex_modbam.sorted.bam",
        out_fp_sorted.to_str().unwrap(),
        "--region",
        "chr17",
        "--force",
    ])
    .unwrap();

    run_modkit(&[
        "extract",
        "tests/resources/duplex_modbam.bam",
        out_fp.to_str().unwrap(),
        "--region",
        "chr17",
        "--force",
    ])
    .unwrap();

    check_mod_profiles_same(&out_fp, &out_fp_sorted)
        .context("test_extract_duplex_correct_output, different outputs with and without index")
        .unwrap();

    check_mod_profiles_same(
        &out_fp,
        &Path::new("tests/resources/duplex_sorted.tsv").to_path_buf(),
    )
    .context("test_extract_correct_output, failed")
    .unwrap();
}

#[test]
fn test_extract_include_sites() {
    let out_fp = std::env::temp_dir().join("test_extract_include_sites.bed");
    let include_bed_fp = "tests/resources/CGI_ladder_3.6kb_ref_CG.bed";
    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "-i",
        "25",
        "--include-bed",
        include_bed_fp,
        "--force",
    ])
    .unwrap();

    let bed_positions =
        parse_bed_file(&Path::new(include_bed_fp).to_path_buf());
    let mod_profile = parse_mod_profile(&out_fp).unwrap();
    for (_read_id, data) in mod_profile {
        for item in data {
            let sites = bed_positions
                .get(&item.contig)
                .expect(&format!("expect to find {}", &item.contig));
            let x = (item.ref_pos, item.strand);
            assert!(sites.contains(&x), "{}", format!("should find {:?}", x))
        }
    }
}

#[test]
fn test_extract_include_sites_duplex_regression() {
    let out_fp =
        std::env::temp_dir().join("test_extract_include_sites_duplex.bed");
    let include_bed_fp = "tests/resources/hg38_chr17_CG0_snip.bed";
    run_modkit(&[
        "extract",
        "tests/resources/duplex_modbam.sorted.bam",
        out_fp.to_str().unwrap(),
        "--include-bed",
        include_bed_fp,
        "--force",
    ])
    .unwrap();
    check_against_expected_text_file(
        out_fp.to_str().unwrap(),
        "tests/resources/test_extract_include_sites_duplex_regression_expected.bed",
    );
}

#[test]
fn test_extract_exclude_sites() {
    let out_fp = std::env::temp_dir().join("test_extract_exclude_sites.bed");
    let exclude_bed_fp = "tests/resources/CGI_ladder_3.6kb_ref_CG_exclude.bed";
    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "-i",
        "25",
        "-v",
        exclude_bed_fp,
        "--force",
    ])
    .unwrap();

    let bed_positions =
        parse_bed_file(&Path::new(exclude_bed_fp).to_path_buf());
    let mod_profile = parse_mod_profile(&out_fp).unwrap();
    for (_read_id, data) in mod_profile {
        for item in data {
            let sites = bed_positions
                .get(&item.contig)
                .expect(&format!("expect to find {}", &item.contig));
            let x = (item.ref_pos, item.strand);
            assert!(
                !sites.contains(&x),
                "{}",
                format!("should not find {:?}", x)
            )
        }
    }
}

#[test]
fn test_extract_unmapped_bam_correct_output() {
    let out_fp = std::env::temp_dir()
        .join("test_extract_unmapped_bam_correct_output.tsv");
    let out_fp_unmapped = std::env::temp_dir()
        .join("test_extract_unmapped_bam_correct_output_unmapped.tsv");
    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.unmapped.bam",
        out_fp_unmapped.to_str().unwrap(),
        "-i",
        "25",
        "--force",
    ])
    .unwrap();

    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "-i",
        "25",
        "--force",
    ])
    .unwrap();
    let mod_profile_unmapped = parse_mod_profile(&out_fp_unmapped)
        .unwrap()
        .into_iter()
        .map(|(read_name, prof)| {
            let mut q_positions =
                prof.into_iter().map(|p| p.q_pos).collect::<Vec<usize>>();
            assert!(!q_positions.is_empty());
            q_positions.sort();
            (read_name, q_positions)
        })
        .collect::<HashMap<String, Vec<usize>>>();
    let mod_profile = parse_mod_profile(&out_fp)
        .unwrap()
        .into_iter()
        .map(|(read_name, prof)| {
            let mut q_positions =
                prof.into_iter().map(|p| p.q_pos).collect::<Vec<usize>>();
            assert!(!q_positions.is_empty());
            q_positions.sort();
            (read_name, q_positions)
        })
        .collect::<HashMap<String, Vec<usize>>>();
    assert_eq!(mod_profile, mod_profile_unmapped);
}

#[test]
fn test_extract_collapse_correct_output() {
    let out_fp =
        std::env::temp_dir().join("test_extract_collapse_correct_output.tsv");
    run_modkit(&[
        "extract",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        out_fp.to_str().unwrap(),
        "--ignore",
        "h",
        "-i",
        "25",
        "--force",
    ])
    .unwrap();

    check_mod_profiles_same(
        &out_fp,
        &Path::new(
            "tests/resources/bc_anchored_10_reads.sorted.methylprofile_ignoreh.tsv",
        )
            .to_path_buf(),
    )
        .context("test_extract_collapse_correct_output, output didn't match")
        .unwrap();
}

#[test]
fn test_extract_implicit_mod_calls() {
    let out_fp =
        std::env::temp_dir().join("test_extract_implicit_mod_calls.tsv");
    run_modkit(&[
        "extract",
        "tests/resources/implicit_mod_tags.bam",
        out_fp.to_str().unwrap(),
        "--force",
    ])
    .unwrap();

    check_mod_profiles_same(
        &out_fp,
        &Path::new("tests/resources/extract_with_implicit.tsv").to_path_buf(),
    )
    .context("test_extract_implicit_mod_calls, output didn't match")
    .unwrap();
}
