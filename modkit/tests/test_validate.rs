use crate::common::run_modkit;
use anyhow::Context;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use tempfile::tempdir;

mod common;

fn run_validate(
    bam_bed_pairs: &[(&Path, &Path)],
    output_file: &Path,
) -> (String, HashMap<String, String>) {
    let mut args = vec!["validate".to_string()];
    for (bam, bed) in bam_bed_pairs {
        args.push("--bam-and-bed".to_string());
        args.push(bam.to_str().unwrap().to_string());
        args.push(bed.to_str().unwrap().to_string());
    }
    args.extend([
        "--out-filepath".to_string(),
        output_file.to_str().unwrap().to_string(),
        "--filter-threshold".to_string(),
        "0".to_string(),
        "--suppress-progress".to_string(),
        "--threads".to_string(),
        "1".to_string(),
    ]);
    let args = args.iter().map(String::as_str).collect::<Vec<_>>();
    run_modkit(&args).unwrap();

    let output = fs::read_to_string(output_file).unwrap();
    let fields = output
        .lines()
        .map(|line| {
            let (key, value) = line.split_once(": ").unwrap();
            (key.to_string(), value.to_string())
        })
        .collect();
    (output, fields)
}

#[test]
fn test_validate_help() {
    run_modkit(&["validate", "--help"])
        .context("validate --help failed")
        .unwrap();
}

#[test]
fn test_validate_expected() {
    let output_file = std::env::temp_dir().join("test_validate_output.tsv");
    run_modkit(&[
        "validate",
        "--bam-and-bed",
        "../tests/resources/input_5mC.bam",
        "../tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed",
        "--bam-and-bed",
        "../tests/resources/input_C.bam",
        "../tests/resources/CGI_ladder_3.6kb_ref_CG_C.bed",
        "--out-filepath",
        output_file.to_str().unwrap(),
    ])
    .context("should run validate")
    .unwrap();

    let reader = BufReader::new(File::open(output_file).unwrap());
    for line in reader.lines() {
        let line_content = line.unwrap();
        if line_content.starts_with("raw_accuracy") {
            let accuracy: f32 = line_content
                .split_whitespace()
                .skip(1)
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or_default();
            assert_eq!(accuracy, 85.305214);
        }
        if line_content.starts_with("filtered_accuracy") {
            let accuracy: f32 = line_content
                .split_whitespace()
                .skip(1)
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or_default();
            assert_eq!(accuracy, 89.00287);
        }
    }
}

#[test]
fn test_validate_reopens_bam_for_each_truth_bed() {
    let temp_dir = tempdir().unwrap();
    let bam = Path::new("../tests/resources/input_5mC.bam");
    let truth = Path::new("../tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed");
    let truth_lines = fs::read_to_string(truth)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect::<Vec<_>>();
    assert_eq!(truth_lines.len(), 48);

    let mut split_a_lines = Vec::new();
    let mut split_b_lines = Vec::new();
    for (idx, line) in truth_lines.iter().enumerate() {
        if idx % 2 == 0 {
            split_a_lines.push(line.clone());
        } else {
            split_b_lines.push(line.clone());
        }
    }
    let truth_set = truth_lines.iter().collect::<HashSet<_>>();
    let split_a_set = split_a_lines.iter().collect::<HashSet<_>>();
    let split_b_set = split_b_lines.iter().collect::<HashSet<_>>();
    assert_eq!(truth_set.len(), truth_lines.len());
    assert!(split_a_set.is_disjoint(&split_b_set));
    assert_eq!(
        split_a_set.union(&split_b_set).copied().collect::<HashSet<_>>(),
        truth_set
    );

    let split_a = temp_dir.path().join("truth-a.bed");
    let split_b = temp_dir.path().join("truth-b.bed");
    fs::write(&split_a, format!("{}\n", split_a_lines.join("\n"))).unwrap();
    fs::write(&split_b, format!("{}\n", split_b_lines.join("\n"))).unwrap();

    let (union_output, union) = run_validate(
        &[(bam, truth)],
        &temp_dir.path().join("union-output.tsv"),
    );
    assert_eq!(
        union["full_contingency_table"],
        "[[\"ground_truth_label\",\"C\",\"h\",\"m\",\"No Call\",\"A\",\"G\",\"T\",\"Deletion\"],[\"m\",16625,20432,120926,9069,4826,8337,859,4608]]"
    );
    assert_eq!(
        union["raw_contingency_table"],
        "[[\"ground_truth_label\",\"C\",\"h\",\"m\"],[\"m\",16625,20432,120926]]"
    );

    let (split_ab_output, split_ab) = run_validate(
        &[(bam, split_a.as_path()), (bam, split_b.as_path())],
        &temp_dir.path().join("split-ab-output.tsv"),
    );
    assert_eq!(
        split_ab["full_contingency_table"],
        union["full_contingency_table"]
    );
    assert_eq!(split_ab, union);
    assert_eq!(split_ab_output, union_output);

    let (split_ba_output, split_ba) = run_validate(
        &[(bam, split_b.as_path()), (bam, split_a.as_path())],
        &temp_dir.path().join("split-ba-output.tsv"),
    );
    assert_eq!(
        split_ba["full_contingency_table"],
        union["full_contingency_table"]
    );
    assert_eq!(split_ba, union);
    assert_eq!(split_ba_output, union_output);
}
