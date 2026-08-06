use crate::common::{
    check_against_expected_text_file, check_legal_csv, run_modkit,
};
use std::fs;
use std::path::Path;

mod common;

#[test]
fn test_dmr_helps() {
    let _ = run_modkit(&["dmr", "pair", "--help"])
        .expect("failed to run modkit dmr pair help");
    let _ = run_modkit(&["dmr", "multi", "--help"])
        .expect("failed to run modkit dmr multi help");
}

#[test]
fn test_dmr_regression() {
    let out_bed = std::env::temp_dir().join("test_dmr_regression.bed");
    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "-f",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );

    let out_bed =
        std::env::temp_dir().join("foo").join("test_dmr_regression_2.bed");

    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "../tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "-f",
        "--header",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed).expect("should be a legal CSV");
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "../tests/resources/test_output_chr20-2.bed",
    );
}

fn run_segmented_dmr(
    output: &Path,
    segments: &Path,
    threads: &str,
    io_threads: &str,
) {
    run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "../tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        output.to_str().unwrap(),
        "--segment",
        segments.to_str().unwrap(),
        "--ref",
        "../tests/resources/GRCh38_chr20.fa",
        "--header",
        "--base",
        "C",
        "--max-coverages",
        "100",
        "100",
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ])
    .expect("segmented DMR run should succeed");
}

#[test]
fn segmentation_includes_last_site_and_is_thread_deterministic() {
    let temp_dir = tempfile::tempdir().unwrap();
    let sites_one = temp_dir.path().join("sites-threads-1.bed");
    let segments_one = temp_dir.path().join("segments-threads-1.bed");
    let sites_four = temp_dir.path().join("sites-threads-4.bed");
    let segments_four = temp_dir.path().join("segments-threads-4.bed");

    run_segmented_dmr(&sites_one, &segments_one, "1", "1");
    run_segmented_dmr(&sites_four, &segments_four, "4", "2");

    let sites_one = fs::read_to_string(sites_one).unwrap();
    let sites_four = fs::read_to_string(sites_four).unwrap();
    let segments_one = fs::read_to_string(segments_one).unwrap();
    let segments_four = fs::read_to_string(segments_four).unwrap();
    assert_eq!(sites_one.as_bytes(), sites_four.as_bytes());
    assert_eq!(segments_one.as_bytes(), segments_four.as_bytes());

    let site_rows = sites_one.lines().skip(1).collect::<Vec<_>>();
    let segment_rows = segments_one.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(site_rows.len(), 17_271);

    let segment_site_total = segment_rows
        .iter()
        .map(|row| row.split('\t').nth(5).unwrap().parse::<usize>().unwrap())
        .sum::<usize>();
    assert_eq!(segment_site_total, site_rows.len());

    let last_site_end = site_rows.last().unwrap().split('\t').nth(2).unwrap();
    let last_segment_end =
        segment_rows.last().unwrap().split('\t').nth(2).unwrap();
    assert_eq!(last_site_end, "10804378");
    assert_eq!(last_segment_end, last_site_end);
}

struct DmrOutput {
    sites: String,
    segments: String,
}

fn run_multi_contig_dmr(
    output_dir: &Path,
    label: &str,
    interval_size: &str,
    threads: &str,
    io_threads: &str,
) -> DmrOutput {
    let sites = output_dir.join(format!("sites-{label}.bed"));
    let segments = output_dir.join(format!("segments-{label}.bed"));

    run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "../tests/resources/dmr_contig_order_a.bed.gz",
        "-b",
        "../tests/resources/dmr_contig_order_b.bed.gz",
        "-o",
        sites.to_str().unwrap(),
        "--segment",
        segments.to_str().unwrap(),
        "--ref",
        "../tests/resources/dmr_contig_order.fa",
        "--header",
        "--base",
        "C",
        "--max-coverages",
        "100",
        "100",
        "--interval-size",
        interval_size,
        "--batch-size",
        "1",
        "--threads",
        threads,
        "--io-threads",
        io_threads,
        "--suppress-progress",
        "--force",
    ])
    .expect("multi-contig segmented DMR run should succeed");

    DmrOutput {
        sites: fs::read_to_string(sites).unwrap(),
        segments: fs::read_to_string(segments).unwrap(),
    }
}

fn parse_site_keys(output: &str) -> Vec<(String, u64, u64)> {
    output
        .lines()
        .skip(1)
        .map(|row| {
            let fields = row.split('\t').collect::<Vec<_>>();
            (
                fields[0].to_string(),
                fields[1].parse::<u64>().unwrap(),
                fields[2].parse::<u64>().unwrap(),
            )
        })
        .collect()
}

fn parse_segments(output: &str) -> Vec<(String, u64, u64, String, usize)> {
    output
        .lines()
        .skip(1)
        .map(|row| {
            let fields = row.split('\t').collect::<Vec<_>>();
            (
                fields[0].to_string(),
                fields[1].parse::<u64>().unwrap(),
                fields[2].parse::<u64>().unwrap(),
                fields[3].to_string(),
                fields[5].parse::<usize>().unwrap(),
            )
        })
        .collect()
}

#[test]
fn multi_contig_segmentation_is_invariant_to_batch_geometry_and_threads() {
    let temp_dir = tempfile::tempdir().unwrap();
    let interval_ten_one =
        run_multi_contig_dmr(temp_dir.path(), "i10-t1", "10", "1", "1");
    let interval_ten_four =
        run_multi_contig_dmr(temp_dir.path(), "i10-t4", "10", "4", "2");
    let interval_three_one =
        run_multi_contig_dmr(temp_dir.path(), "i3-t1", "3", "1", "1");
    let interval_three_four =
        run_multi_contig_dmr(temp_dir.path(), "i3-t4", "3", "4", "2");

    for (label, other) in [
        ("interval 10, four threads", &interval_ten_four),
        ("interval 3, one thread", &interval_three_one),
        ("interval 3, four threads", &interval_three_four),
    ] {
        assert!(
            interval_ten_one.sites.as_bytes() == other.sites.as_bytes(),
            "site output changed for {label}"
        );
        assert!(
            interval_ten_one.segments.as_bytes() == other.segments.as_bytes(),
            "segment output changed for {label}"
        );
    }

    let site_keys = parse_site_keys(&interval_ten_one.sites);
    let expected_site_keys = [("aa", 6u64), ("bb", 15u64), ("cc", 6u64)]
        .into_iter()
        .flat_map(|(chrom, size)| {
            (0..size).map(move |position| {
                (chrom.to_string(), position, position + 1)
            })
        })
        .collect::<Vec<_>>();
    assert_eq!(site_keys, expected_site_keys);

    let segments = parse_segments(&interval_ten_one.segments);
    assert_eq!(
        segments,
        vec![
            ("aa".to_string(), 0, 6, "different".to_string(), 6),
            ("bb".to_string(), 0, 15, "different".to_string(), 15),
            ("cc".to_string(), 0, 6, "different".to_string(), 6),
        ]
    );

    for (chrom, start, end) in site_keys {
        let covering_segments = segments
            .iter()
            .filter(|(segment_chrom, segment_start, segment_end, _, _)| {
                segment_chrom == &chrom
                    && *segment_start <= start
                    && end <= *segment_end
            })
            .count();
        assert_eq!(covering_segments, 1, "site {chrom}:{start}-{end}");
    }
}

// todo
//  test pair with explicit index
//  test multi
