use crate::common::run_simple_summary;
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
