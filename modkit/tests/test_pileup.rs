use anyhow::Context;
use itertools::Itertools;
use rust_htslib::bam;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use common::{check_against_expected_text_file, run_modkit};
use mod_kit::dmr::bedmethyl::BedMethylLine;
use mod_kit::mod_base_code::{ModCodeRepr, METHYL_CYTOSINE};

mod common;

fn make_ambiguous_query_pileup_bam(bam_path: &PathBuf, overrun: bool) {
    let mut header = bam::Header::new();
    let mut header_record = bam::header::HeaderRecord::new(b"HD");
    header_record.push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate");
    header.push_record(&header_record);
    let mut reference_record = bam::header::HeaderRecord::new(b"SQ");
    reference_record
        .push_tag(b"SN", "oligo_1512_adapters")
        .push_tag(b"LN", 156);
    header.push_record(&reference_record);

    let (sequence, cigar, mm_tag, ml_tag) = if overrun {
        (
            b"CNCTGTACTT".as_slice(),
            bam::record::CigarString(vec![
                bam::record::Cigar::SoftClip(1),
                bam::record::Cigar::Match(9),
            ]),
            "C+m?,0,0;",
            vec![255, 255],
        )
    } else {
        (
            b"NCTGTACTTN".as_slice(),
            bam::record::CigarString(vec![bam::record::Cigar::Match(10)]),
            "C+m?,0;",
            vec![255],
        )
    };
    let mut record = bam::Record::new();
    record.set(b"ambiguous-query", Some(&cigar), sequence, &[30; 10]);
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    record.push_aux(b"MM", bam::record::Aux::String(mm_tag)).unwrap();
    record
        .push_aux(b"ML", bam::record::Aux::ArrayU8((&ml_tag[..]).into()))
        .unwrap();

    let mut writer =
        bam::Writer::from_path(bam_path, &header, bam::Format::Bam).unwrap();
    writer.write(&record).unwrap();
    drop(writer);
    bam::index::build(bam_path, None, bam::index::Type::Bai, 1).unwrap();
}

type AmbiguousQueryPileupRow = (u64, ModCodeRepr, u64, u64);

fn run_ambiguous_query_pileup(
    input_bam: &PathBuf,
    output_bed: &PathBuf,
    optimized: bool,
) -> Result<Vec<AmbiguousQueryPileupRow>, String> {
    let mut args = vec![
        "pileup",
        input_bam.to_str().unwrap(),
        output_bed.to_str().unwrap(),
    ];
    if optimized {
        args.extend(["--modified-bases", "5mC"]);
    }
    args.extend([
        "--no-filtering",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--threads",
        "1",
    ]);
    run_modkit(&args).map_err(|error| error.to_string())?;

    Ok(BufReader::new(File::open(output_bed).unwrap())
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .map(|row| {
            (
                row.start(),
                row.raw_mod_code,
                row.count_methylated,
                row.valid_coverage,
            )
        })
        .collect())
}

#[test]
fn test_pileup_help() {
    let pileup_help_args = ["pileup", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

#[test]
fn test_pileup_no_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_nofilt.methyl.bed",
    );
}

#[test]
fn test_pileup_with_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_withfilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-f",
        "1.0",
        "-p",
        "0.25",
        "--seed",
        "42",
        "--include-unmapped",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt025.methyl.bed",
    );
}

#[test]
fn test_pileup_combine() {
    let test_adjusted_bam = std::env::temp_dir().join("test_combined.bed");
    let pileup_args = [
        "pileup",
        "--combine-mods",
        "--no-filtering",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        test_adjusted_bam.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_adjusted_bam.exists());

    check_against_expected_text_file(
        test_adjusted_bam.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_combined.methyl.bed",
    );
}

#[test]
#[ignore = "collapse no longer available in pileup since v0.6.0"]
fn test_pileup_collapse() {
    let test_collapsed_bam = std::env::temp_dir().join("test_collapsed.bam");
    let test_collapsed_bed = std::env::temp_dir().join("test_collapsed.bed");
    let test_restricted_bed = std::env::temp_dir().join("test_restricted.bed");

    let collapse_args = [
        "adjust-mods",
        "--ignore",
        "h",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        test_collapsed_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();
    assert!(test_collapsed_bam.exists());
    bam::index::build(
        test_collapsed_bam.clone(),
        None,
        bam::index::Type::Bai,
        1,
    )
    .unwrap();

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        test_collapsed_bam.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_collapsed_bed.exists());

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--ignore",
        "h",
        "--no-filtering",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        test_restricted_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_restricted_bed.exists());
    check_against_expected_text_file(
        test_restricted_bed.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_no_mod_calls() {
    let empty_bedfile =
        std::env::temp_dir().join("test_pileup_no_mod_calls_outbed.bed");
    let args = [
        "pileup",
        "--no-filtering",
        "../tests/resources/empty-tags.sorted.bam",
        empty_bedfile.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    let reader = BufReader::new(File::open(empty_bedfile).unwrap());
    let lines = reader.lines().collect::<Vec<Result<String, _>>>();
    assert_eq!(lines.len(), 0);
}

#[test]
fn test_pileup_skips_ambiguous_query_and_keeps_downstream_modification() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_bam = temp_dir.path().join("ambiguous-query.bam");
    let optimized_output = temp_dir.path().join("optimized.bed");
    let generic_output = temp_dir.path().join("generic.bed");
    make_ambiguous_query_pileup_bam(&input_bam, false);
    let observed = vec![
        (
            "optimized",
            run_ambiguous_query_pileup(&input_bam, &optimized_output, true),
        ),
        (
            "generic",
            run_ambiguous_query_pileup(&input_bam, &generic_output, false),
        ),
    ];
    let expected_row = (1, METHYL_CYTOSINE, 1, 1);

    assert_eq!(
        observed,
        vec![
            ("optimized", Ok(vec![expected_row])),
            ("generic", Ok(vec![expected_row])),
        ]
    );
}

#[test]
fn test_pileup_skips_ambiguous_query_after_unaligned_modification() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_bam = temp_dir.path().join("ambiguous-query-overrun.bam");
    let optimized_output = temp_dir.path().join("optimized.bed");
    let generic_output = temp_dir.path().join("generic.bed");
    make_ambiguous_query_pileup_bam(&input_bam, true);
    let observed = vec![
        (
            "optimized",
            run_ambiguous_query_pileup(&input_bam, &optimized_output, true),
        ),
        (
            "generic",
            run_ambiguous_query_pileup(&input_bam, &generic_output, false),
        ),
    ];
    let expected_row = (1, METHYL_CYTOSINE, 1, 1);

    assert_eq!(
        observed,
        vec![
            ("optimized", Ok(vec![expected_row])),
            ("generic", Ok(vec![expected_row])),
        ]
    );
}

#[test]
fn test_pileup_old_tags() {
    let updated_file =
        std::env::temp_dir().join("test_pileup_old_tags_updated.bam");
    run_modkit(&[
        "update-tags",
        "../tests/resources/HG002_small.ch20._other.sorted.bam",
        "--mode",
        "ambiguous",
        "--no-implicit-probs",
        updated_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(updated_file.exists());
    bam::index::build(updated_file.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    let out_file = std::env::temp_dir().join("test_pileup_old_tags.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        updated_file.to_str().unwrap(),
        out_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(out_file.exists());
    check_against_expected_text_file(
        out_file.to_str().unwrap(),
        "../tests/resources/pileup-old-tags-regressiontest.methyl.bed",
    );
}

#[test]
fn test_pileup_with_region() {
    let temp_file = std::env::temp_dir().join("test_pileup_with_region.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--region",
        "oligo_1512_adapters:0-50",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_nofilt_oligo_1512_adapters_10_50.\
         bed",
    );
}

#[test]
fn test_pileup_duplex_reads() {
    let temp_file = std::env::temp_dir().join("test_pileup_duplex_reads.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/duplex_modbam.sorted.bam",
        temp_file.to_str().unwrap(),
        "--region",
        "chr17",
        "--no-filtering",
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/duplex_modbam_pileup_nofilt.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering() {
    let temp_file = std::env::temp_dir().join("test_cpg_motif_filtering.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering_compressed_ref() {
    let temp_file = std::env::temp_dir()
        .join("test_pileup_cpg_motif_filtering_compressed_ref.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa.gz",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering_strand_combine() {
    let temp_file = std::env::temp_dir()
        .join("test_cpg_motif_filtering_strand_combine.bed");
    let interval_sizes =
        ["10", "88", "89", "90", "91", "92", "93", "94", "10000"];
    for interval_size in interval_sizes {
        run_modkit(&[
            "pileup",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
            temp_file.to_str().unwrap(),
            "--no-filtering",
            "-i",
            interval_size,
            "--cpg",
            "--combine-strands",
            "--modified-bases",
            "5mC",
            "5hmC",
            "--ref",
            "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        ])
        .unwrap();
        check_against_expected_text_file(
            temp_file.to_str().unwrap(),
            "../tests/resources/\
             bc_anchored_10_reads_nofilt_cg_motif_strand_combine.bed",
        );
    }
}

#[test]
#[ignore = "ignored bases not supported in puleup since v0.6.0"]
fn test_pileup_presets_traditional_same_as_options() {
    let preset_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options.bed");
    let options_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options2.bed");

    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        preset_temp_file.to_str().unwrap(),
        "--no-filtering",
        "--mixed-delim",
        "--preset",
        "traditional",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        options_temp_file.to_str().unwrap(),
        "--cpg",
        "--no-filtering",
        "--mixed-delim",
        "--ignore",
        "h",
        "--combine-strands",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        preset_temp_file.to_str().unwrap(),
        options_temp_file.to_str().unwrap(),
    );
}

#[test]
#[ignore = "duplicated reads no longer checked since v0.6.0"]
fn test_pileup_duplicated_reads_ignored() {
    let control_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_control.bed");
    let test_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_marked.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        control_fp.to_str().unwrap(),
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "../tests/resources/duplicated.marked.fixed.bam",
        test_fp.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        control_fp.to_str().unwrap(),
        test_fp.to_str().unwrap(),
    );
}

#[test]
#[ignore = "TODO: refactor this test"]
fn test_pileup_edge_filter_regression() {
    let adjusted_bam =
        std::env::temp_dir().join("test_pileup_edge_filter_adjusted.bam");
    let edge_filter_bed = std::env::temp_dir()
        .join("test_pileup_edge_filter_edge_filter_50.pileup.bed");
    let edge_filter_bed_2 = std::env::temp_dir()
        .join("test_pileup_edge_filter_edge_filter_50_2.pileup.bed");

    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--mixed-delim",
        "--edge-filter",
        "50",
    ])
    .context("test_pileup_edge_filter_regression failed to make bedMethyl")
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_edge_filter50.bed",
    );

    run_modkit(&[
        "adjust-mods",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        adjusted_bam.to_str().unwrap(),
        "--edge-filter",
        "50",
    ])
    .context("test_pileup_edge_filter_regression failed to run adjust-mods")
    .unwrap();

    bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    run_modkit(&[
        "pileup",
        adjusted_bam.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
        "--mixed-delim",
        "--no-filtering",
        "--edge-filter",
        "50",
    ])
    .context(
        "test_pileup_edge_filter_regression failed to make bedMethyl on \
         adjusted bam",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
    );
}

#[test]
#[ignore = "TODO: refactor this test"]
fn test_pileup_edge_filter_asymmetric_regression() {
    let adjusted_bam = std::env::temp_dir()
        .join("test_pileup_edge_filter_asymmetric_regression.bam");
    let edge_filter_bed = std::env::temp_dir().join(
        "test_pileup_edge_filter_asymmetric_regression_filter_50.pileup.bed",
    );
    let edge_filter_bed_2 = std::env::temp_dir()
        .join("test_pileup_edge_filter_asymmetric_regression_50_2.pileup.bed");

    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--edge-filter",
        "50,50",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make \
         bedMethyl with 50,50",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_edge_filter50.bed",
    );

    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--edge-filter",
        "50,0",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make \
         bedMethyl with 50,0",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_edge_filter50-0.bed",
    );

    run_modkit(&[
        "adjust-mods",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        adjusted_bam.to_str().unwrap(),
        "--edge-filter",
        "50,0",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make adjust \
         mods BAM",
    )
    .unwrap();

    bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    run_modkit(&[
        "pileup",
        adjusted_bam.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
        "--no-filtering",
        "--mixed-delim",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make pileup \
         on adjusted BAM",
    )
    .unwrap();

    check_against_expected_text_file(
        edge_filter_bed_2.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_edge_filter50-0.bed",
    );
}

#[test]
#[ignore = "partition tags removed in v0.6.0"]
fn test_pileup_partition_tags_partitioned() {
    let tmp_dir =
        std::env::temp_dir().join("test_pileup_partition_tags_partitioned");
    let control_file =
        std::env::temp_dir().join("test_pileup_partition_tags_control.bed");

    // control BED, all of the partitioned BED files should be the same as this
    // one
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        control_file.to_str().unwrap(),
        "--no-filtering",
    ])
    .context("failed to run modkit on control")
    .unwrap();

    // run partitioned on HP and RG tags. This test file has 2 HP tags {1, 2}
    // and 3 read groups {A, B, C}. So we expect 6 files, all the same as the
    // control
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        tmp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--no-filtering",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();

    let mut count = 0;
    for result in tmp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            control_file.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
#[ignore = "partition tags removed in v0.6.0"]
fn test_pileup_partition_tags_bedgraph() {
    let tmp_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_bedgraph_partitioned");
    let control_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_bedgraph_control");

    let collect_bedgraph_files =
        |dir_path: &PathBuf| -> std::io::Result<Vec<PathBuf>> {
            dir_path.read_dir().map(|read_dir| {
                read_dir
                    .filter_map(|dir| match dir {
                        Ok(dir) => {
                            if dir.path().extension().and_then(|fp| fp.to_str())
                                == Some("bedgraph")
                            {
                                Some(dir.path())
                            } else {
                                None
                            }
                        }
                        Err(_) => None,
                    })
                    .collect::<Vec<PathBuf>>()
            })
        };

    // control BED, all of the partitioned BED files should be the same as this
    // one
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        control_dir.to_str().unwrap(),
        "--no-filtering",
        "--bedgraph",
    ])
    .context("failed to run modkit on control bedgraph")
    .unwrap();

    let control_bedgraph_files = collect_bedgraph_files(&control_dir)
        .unwrap()
        .into_iter()
        .map(|fp| {
            let file_name = fp.file_name().unwrap().to_str().unwrap();
            match (file_name.starts_with("h"), file_name.contains("positive")) {
                (true, true) => (('h', "positive"), fp),
                (true, false) => (('h', "negative"), fp),
                (false, true) => (('m', "positive"), fp),
                (false, false) => (('m', "negative"), fp),
            }
        })
        .collect::<HashMap<(char, &str), PathBuf>>();

    // run partitioned on HP and RG tags. This test file has 2 HP tags {1, 2}
    // and 3 read groups {A, B, C}. So we expect 6 files, all the same as the
    // control
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        tmp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--no-filtering",
        "--bedgraph",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();

    let mut count = 0;
    for result in tmp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        if dir_entry.extension().and_then(|s| s.to_str()) != Some("bedgraph") {
            continue;
        }
        let file_name = dir_entry.file_name().unwrap().to_str().unwrap();
        let stripped = file_name.replace(".bedgraph", "");
        let parts = stripped.split('_').collect::<Vec<&str>>();
        let mod_code = parts[2].parse::<char>().unwrap();
        let strand = parts[3];
        let key = (mod_code, strand);
        let file_to_compare_to = control_bedgraph_files.get(&key).unwrap();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            file_to_compare_to.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 24);
}

#[test]
fn test_pileup_with_filt_position_filter() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_with_filt_position_filter.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--filter-threshold",
        "0.6894531",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "../tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt_positions_fixed_thresh.\
         methyl.bed",
    );
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-p",
        "0.25",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "../tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt_positions_025.methyl.bed",
    );
}

#[test]
#[ignore = "traditional preset removed in v0.6.0"]
fn test_pileup_with_filter_positions_and_traditional() {
    let temp_file = std::env::temp_dir()
        .join("test_pileup_with_filter_positions_and_traditional.bed");
    run_modkit(&[
        "pileup",
        "--mixed-delimiters",
        "-i",
        "25", // use small interval to make sure chunking works
        "-p",
        "0.25",
        "--preset",
        "traditional",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "../tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt_positions_025_traditional.\
         methyl.bed",
    );
}

#[test]
#[ignore = "partition tags removed in v0.6.0"]
fn test_pileup_partition_tags_combine_strands() {
    let exp_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_combine_strands_partitioned");
    let control_file = std::env::temp_dir()
        .join("test_pileup_partition_tags_combine_strands_control.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        control_file.to_str().unwrap(),
        "--combine-strands",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--cpg",
        "--no-filtering",
    ])
    .context("failed to run modkit on control")
    .unwrap();
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        exp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--combine-strands",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--cpg",
        "--no-filtering",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();
    let mut count = 0;
    for result in exp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            control_file.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
fn test_pileup_motifs_cg0_cgcg2() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_motifs_cg0_cgcg2.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test1.bed",
    );

    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads-2.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test2.bed",
    );
}

#[test]
#[ignore = "multiple motifs and combine strands not supported in v0.6.0"]
fn test_pileup_motifs_cg0_cgcg2_combined() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_motifs_cg0_cgcg2_combined.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--combine-strands",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test1_combine_strands.bed",
    );

    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads-2.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--combine-strands",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test2_combine_strands.bed",
    );
}

#[test]
fn test_pileup_chebi_code_same_output() {
    let ord_bm_line = |a: &BedMethylLine, b: &BedMethylLine| -> Ordering {
        match a.chrom.cmp(&b.chrom) {
            Ordering::Equal => match a.start().cmp(&b.start()) {
                Ordering::Equal => a.raw_mod_code.cmp(&b.raw_mod_code),
                o @ _ => o,
            },
            o @ _ => o,
        }
    };
    let adjusted_bam =
        std::env::temp_dir().join("test_chebi_code_same_output_hmc2chEBI.bam");
    let pileup =
        std::env::temp_dir().join("test_chebi_code_same_output_pileup.bed");
    let expected_fp = "../tests/resources/modbam.modpileup_nofilt.methyl.bed";
    let expected = BufReader::new(File::open(expected_fp).unwrap())
        .lines()
        .map(|l| BedMethylLine::parse(&l.unwrap()).unwrap())
        .sorted_by(|a, b| ord_bm_line(a, b))
        .collect::<Vec<BedMethylLine>>();
    for to_code in [ModCodeRepr::ChEbi(76792), ModCodeRepr::Code('c')] {
        run_modkit(&[
            "adjust-mods",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
            adjusted_bam.to_str().unwrap(),
            "--convert",
            "h",
            &to_code.to_string(),
        ])
        .with_context(|| format!("failed to change 5hmC to {to_code}"))
        .unwrap();

        bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
            .unwrap();

        run_modkit(&[
            "pileup",
            adjusted_bam.to_str().unwrap(),
            pileup.to_str().unwrap(),
            "-i",
            "25", // use small interval to make sure chunking works
            "--no-filtering",
        ])
        .context("failed to generate pileup")
        .unwrap();

        let observed = BufReader::new(File::open(pileup.clone()).unwrap())
            .lines()
            .map(|l| {
                let bm = BedMethylLine::parse(&l.unwrap()).unwrap();
                if bm.raw_mod_code != METHYL_CYTOSINE {
                    assert_eq!(bm.raw_mod_code, to_code);
                    BedMethylLine::new(
                        bm.chrom,
                        bm.interval,
                        ModCodeRepr::Code('h'),
                        bm.strand,
                        bm.count_methylated,
                        bm.valid_coverage,
                        bm.count_canonical,
                        bm.count_other,
                        bm.count_delete,
                        bm.count_fail,
                        bm.count_diff,
                        bm.count_nocall,
                    )
                } else {
                    bm
                }
            })
            .sorted_by(|a, b| ord_bm_line(a, b))
            .collect::<Vec<BedMethylLine>>();
        assert_eq!(expected, observed);
    }
}

#[test]
fn test_pileup_with_header() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--with-header",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/pileup_with_header.bed",
    );
}

#[test]
fn test_pilep_phased() {
    let test_tmp_dir = std::env::temp_dir();
    let phased_output_dir = test_tmp_dir.join("test_pilep_phased_pileup");
    let normal_output_pileup =
        test_tmp_dir.join("test_pilep_regular_pileup.bedmethyl");
    let args = [
        "pileup",
        "../tests/resources/test.sorted.phased.bam",
        phased_output_dir.to_str().unwrap(),
        "--modified-bases",
        "5mC",
        "--phased",
        "--ref",
        "../tests/resources/genome_for_phased_test.fasta",
        "--no-filtering",
    ];
    run_modkit(&args).unwrap();
    let args = [
        "pileup",
        "../tests/resources/test.sorted.phased.bam",
        normal_output_pileup.to_str().unwrap(),
        "--modified-bases",
        "5mC",
        "--ref",
        "../tests/resources/genome_for_phased_test.fasta",
        "--no-filtering",
    ];
    run_modkit(&args).unwrap();
    let combined_phased = phased_output_dir.join("combined.bedmethyl");
    assert!(combined_phased.exists() && combined_phased.is_file());
    assert!(normal_output_pileup.exists() && normal_output_pileup.is_file());
    check_against_expected_text_file(
        combined_phased.to_str().unwrap(),
        normal_output_pileup.to_str().unwrap(),
    );
}
