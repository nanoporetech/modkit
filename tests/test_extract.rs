use crate::common::check_against_expected_text_file;
use anyhow::{anyhow, Context};
use common::run_modkit;
use derive_new::new;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::{Path, PathBuf};

mod common;

#[derive(new, Eq, PartialEq, Debug)]
struct ModData {
    q_pos: usize,
    ref_pos: i64,
    mod_code: char,
    strand: char,
    contig: String,
    data: String,
}

impl PartialOrd for ModData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ModData {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.q_pos.cmp(&other.q_pos) {
            Ordering::Equal => match self.mod_code.cmp(&other.mod_code) {
                Ordering::Equal => self.strand.cmp(&other.strand),
                ord => ord,
            },
            ord => ord,
        }
    }
}

fn parse_mod_profile(
    fp: &PathBuf,
) -> anyhow::Result<HashMap<String, Vec<ModData>>> {
    let mut reader =
        BufReader::new(File::open(fp)?).lines().map(|l| l.unwrap());
    let mut agg = HashMap::new();
    let _ = reader.next(); // discard header
    while let Some(line) = reader.next() {
        let parts = line.split_ascii_whitespace().collect::<Vec<&str>>();
        let read_id = parts[0].to_owned();
        let q_pos = parts[1].parse::<usize>().unwrap();
        let ref_pos = parts[2].parse::<i64>().unwrap();
        let mod_code = parts[10].parse::<char>().unwrap();
        let strand = parts[5].parse::<char>().unwrap();
        let contig = parts[3].to_owned();
        agg.entry(read_id)
            .or_insert(Vec::new())
            .push(ModData::new(q_pos, ref_pos, mod_code, strand, contig, line));
    }
    for (_, dat) in agg.iter_mut() {
        dat.sort()
    }

    Ok(agg)
}

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
        ModData::new(0, 1, 'm', '+', "".to_string(), "".to_string());
    let mod_data2 =
        ModData::new(0, 1, 'h', '+', "".to_string(), "".to_string());
    let mod_data3 =
        ModData::new(1, 1, 'h', '+', "".to_string(), "".to_string());
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
                        "read {read}'s profile is different expected"
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
        "--include",
        include_bed_fp,
        "--force",
        "--mapped",
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
