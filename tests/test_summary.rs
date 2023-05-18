use crate::common::{
    run_simple_summary, run_simple_summary_with_collapse_method,
};
use mod_kit::mod_bam::CollapseMethod;
use mod_kit::mod_base_code::ModCode;
use std::collections::HashSet;
use std::path::Path;

mod common;

/// tests that the summary from a BAM is the same if run with the BAI
/// or without (just taking the first N reads). In this case we use all
/// of the reads in the BAM
#[test]
fn test_summary_with_regions() {
    let tf = std::env::temp_dir().join("test_summary_with_regions.bam");
    let bam_fp = Path::new("tests/resources/bc_anchored_10_reads.sorted.bam");
    std::fs::copy(bam_fp, tf.clone()).unwrap();

    let summary_with_index =
        run_simple_summary(bam_fp.to_str().unwrap(), 25).unwrap();
    let summary_without_index =
        run_simple_summary(tf.to_str().unwrap(), 25).unwrap();

    assert_eq!(summary_with_index, summary_without_index);
}

#[test]
fn test_summary_ignore() {
    let bam_fp = Path::new("tests/resources/bc_anchored_10_reads.sorted.bam");
    let summary_wo_collapse =
        run_simple_summary(bam_fp.to_str().unwrap(), 25).unwrap();
    let summary_w_collapse = run_simple_summary_with_collapse_method(
        bam_fp.to_str().unwrap(),
        25,
        &CollapseMethod::ReDistribute('h'),
    )
    .unwrap();

    let mod_codes = summary_wo_collapse
        .mod_call_counts
        .values()
        .flat_map(|mod_code_counts| {
            mod_code_counts
                .keys()
                .map(|mc| *mc)
                .collect::<Vec<ModCode>>()
        })
        .collect::<HashSet<ModCode>>();
    let expected = vec![ModCode::C, ModCode::m, ModCode::h]
        .into_iter()
        .collect::<HashSet<_>>();
    assert_eq!(mod_codes, expected);
    let mod_codes = summary_w_collapse
        .mod_call_counts
        .values()
        .flat_map(|mod_code_counts| {
            mod_code_counts
                .keys()
                .map(|mc| *mc)
                .collect::<Vec<ModCode>>()
        })
        .collect::<HashSet<ModCode>>();
    let expected = vec![ModCode::C, ModCode::m]
        .into_iter()
        .collect::<HashSet<_>>();
    assert_eq!(mod_codes, expected);
}
