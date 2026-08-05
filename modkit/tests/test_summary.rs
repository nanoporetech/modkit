use crate::common::{
    run_modkit, run_simple_summary, run_simple_summary_with_collapse_method,
    run_simple_summary_with_edge_filter, run_summary_with_include_positions,
};
use anyhow::Context;
use mod_kit::mod_bam::{CollapseMethod, EdgeFilter};
use mod_kit::mod_base_code::{BaseState, DnaBase};
use rust_htslib::bam::{self, record::Aux, Read};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::process::Command;
use tempfile::tempdir;

mod common;

fn run_cpg_summary(bam_fp: &Path) -> HashMap<String, String> {
    let output = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "summary",
            bam_fp.to_str().unwrap(),
            "--reference",
            "../tests/resources/CGI_ladder_3.6kb_ref.fa",
            "--cpg",
            "--filter-threshold",
            "0",
            "--tsv",
            "--suppress-progress",
            "--threads",
            "1",
            "--io-threads",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "summary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| {
            let (key, value) = line.split_once('\t').unwrap();
            (key.to_string(), value.to_string())
        })
        .collect()
}

#[test]
fn test_summary_cpg_excludes_unmapped_calls() {
    let mapped_bam =
        Path::new("../tests/resources/bc_anchored_10_reads.sorted.bam");
    let temp_dir = tempdir().unwrap();
    let mixed_bam = temp_dir.path().join("mapped-and-unmapped.bam");
    let mut reader = bam::Reader::from_path(mapped_bam).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut writer =
        bam::Writer::from_path(&mixed_bam, &header, bam::Format::Bam).unwrap();
    for record in reader.records() {
        writer.write(&record.unwrap()).unwrap();
    }

    let mut unmapped = bam::Record::new();
    unmapped.set(b"unmapped", None, b"CCCC", &[255; 4]);
    unmapped.set_unmapped();
    unmapped.set_tid(-1);
    unmapped.set_pos(-1);
    unmapped.set_mtid(-1);
    unmapped.set_mpos(-1);
    unmapped.push_aux(b"MM", Aux::String("C+m?,0,0,0,0;")).unwrap();
    unmapped.push_aux(b"ML", Aux::ArrayU8((&[255; 4][..]).into())).unwrap();
    unmapped.push_aux(b"MN", Aux::I32(4)).unwrap();
    writer.write(&unmapped).unwrap();
    drop(writer);
    bam::index::build(&mixed_bam, None, bam::index::Type::Bai, 1).unwrap();

    let mixed_summary = run_cpg_summary(&mixed_bam);
    let mapped_summary = run_cpg_summary(mapped_bam);
    assert_eq!(mixed_summary["total_reads_used"], "10");
    assert_eq!(mixed_summary["count_reads_C"], "10");
    assert_eq!(mixed_summary["C_total_mod_calls"], "77");
    assert_eq!(mixed_summary, mapped_summary);
}

/// tests that the summary from a BAM is the same if run with the BAI
/// or without (just taking the first N reads). In this case we use all
/// of the reads in the BAM
#[test]
fn test_summary_with_regions() {
    let tf = std::env::temp_dir().join("test_summary_with_regions.bam");
    let bam_fp =
        Path::new("../tests/resources/bc_anchored_10_reads.sorted.bam");
    std::fs::copy(bam_fp, tf.clone()).unwrap();

    let summary_with_index =
        run_simple_summary(bam_fp.to_str().unwrap(), 25).unwrap();
    let summary_without_index =
        run_simple_summary(tf.to_str().unwrap(), 25).unwrap();

    assert_eq!(summary_with_index, summary_without_index);
}

#[test]
fn test_summary_ignore() {
    let bam_fp =
        Path::new("../tests/resources/bc_anchored_10_reads.sorted.bam");
    let summary_wo_collapse =
        run_simple_summary(bam_fp.to_str().unwrap(), 25).unwrap();
    let summary_w_collapse = run_simple_summary_with_collapse_method(
        bam_fp.to_str().unwrap(),
        25,
        &CollapseMethod::ReDistribute('h'.into()),
    )
    .unwrap();

    let mod_codes = summary_wo_collapse
        .mod_call_counts
        .values()
        .flat_map(|mod_code_counts| {
            mod_code_counts.keys().map(|mc| *mc).collect::<Vec<BaseState>>()
        })
        .collect::<HashSet<BaseState>>();
    let expected = HashSet::from([
        BaseState::Canonical(DnaBase::C),
        BaseState::Modified('m'.into()),
        BaseState::Modified('h'.into()),
    ]);

    assert_eq!(mod_codes, expected);
    let mod_codes = summary_w_collapse
        .mod_call_counts
        .values()
        .flat_map(|mod_code_counts| {
            mod_code_counts.keys().map(|mc| *mc).collect::<Vec<BaseState>>()
        })
        .collect::<HashSet<BaseState>>();
    let expected = HashSet::from([
        BaseState::Canonical(DnaBase::C),
        BaseState::Modified('m'.into()),
    ]);

    assert_eq!(mod_codes, expected);
}

#[test]
fn test_summary_edge_filter() {
    let bam_fp =
        Path::new("../tests/resources/bc_anchored_10_reads.sorted.bam");
    let summary_wo_edge_filter =
        run_simple_summary(bam_fp.to_str().unwrap(), 25)
            .context("test_summary_edge_filter failed to make control summary")
            .unwrap();

    let trim = 50;
    let edge_filter = EdgeFilter::new(trim, trim, false);
    let summary_w_edge_filter = run_simple_summary_with_edge_filter(
        bam_fp.to_str().unwrap(),
        25,
        &edge_filter,
    )
    .context("test_summary_edge_filter failed to make summary with edge filter")
    .unwrap();
    assert_eq!(
        summary_w_edge_filter.reads_with_mod_calls.get(&DnaBase::C).unwrap(),
        summary_wo_edge_filter.reads_with_mod_calls.get(&DnaBase::C).unwrap()
    );
    assert_eq!(
        summary_w_edge_filter.total_reads_used,
        summary_wo_edge_filter.total_reads_used
    );
    let total_mod_calls_wo_filter = summary_wo_edge_filter
        .mod_call_counts
        .get(&DnaBase::C)
        .unwrap()
        .values()
        .sum::<u64>();
    let total_mod_calls_w_filter = summary_w_edge_filter
        .mod_call_counts
        .get(&DnaBase::C)
        .unwrap()
        .values()
        .sum::<u64>();
    // weak assertion but better than nothing.
    assert!(total_mod_calls_wo_filter > total_mod_calls_w_filter);

    let adjusted_bam =
        std::env::temp_dir().join("test_summary_edge_filter_adjusted.bam");
    run_modkit(&[
        "adjust-mods",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        adjusted_bam.to_str().unwrap(),
        "--edge-filter",
        &format!("{}", trim),
    ])
    .context("test_summary_edge_filter failed to run adjust-mods")
    .unwrap();
    let summary_on_adjusted =
        run_simple_summary(adjusted_bam.to_str().unwrap(), 25)
            .context(
                "test_summary_edge_filter failed to make summary on adjusted \
                 bam",
            )
            .unwrap();
    assert_eq!(summary_w_edge_filter, summary_on_adjusted);
}

#[test]
fn test_summary_implicit_calls() {
    let bam_fp = Path::new("../tests/resources/single_read.bam").to_path_buf();
    let bedpositions =
        Path::new("../tests/resources/include_bed_summary_test.bed")
            .to_path_buf();

    let summary =
        run_summary_with_include_positions(&bam_fp, &bedpositions).unwrap();
    assert_eq!(
        summary
            .mod_call_counts
            .get(&DnaBase::A)
            .unwrap()
            .get(&BaseState::Canonical(DnaBase::A))
            .unwrap(),
        &8
    );
    assert_eq!(summary.reads_with_mod_calls.get(&DnaBase::A).unwrap(), &1);
    assert_eq!(summary.total_reads_used, 1);
    // let expected = ModSummary {
    //     reads_with_mod_calls: {
    //         A: 1,
    //     },
    //     mod_call_counts: {
    //         A: {
    //             A: 8,
    //         },
    //     },
    //     filtered_mod_call_counts: {
    //         A: {},
    //     },
    //     total_reads_used: 1,
    //     per_base_thresholds: {},
    //     region: None,
    // }
}
