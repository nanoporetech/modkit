use anyhow::Context;
use itertools::Itertools;
use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{Format, Header, Record, Writer as BamWriter};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use common::{check_against_expected_text_file, run_modkit};
use mod_kit::dmr::bedmethyl::BedMethylLine;
use mod_kit::mod_base_code::{ModCodeRepr, METHYL_CYTOSINE};

mod common;

fn read_bed_sites(path: &str) -> HashSet<(String, u32, char)> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let fields = line.split('\t').collect::<Vec<_>>();
            (
                fields[0].to_string(),
                fields[1].parse::<u32>().unwrap(),
                fields[5].parse::<char>().unwrap(),
            )
        })
        .collect()
}

fn write_dynamic_slot_fixture_with_probs(
    root: &Path,
    name: &str,
    probabilities: &[[u8; 3]],
) -> (PathBuf, PathBuf) {
    let bam_path = root.join(format!("{name}.bam"));
    let fasta_path = root.join("reference.fa");

    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 3);
    header.push_record(&sq);

    let mut writer =
        BamWriter::from_path(&bam_path, &header, Format::Bam).unwrap();
    for (i, probabilities) in probabilities.iter().enumerate() {
        let cigar = CigarString(vec![Cigar::Match(3)]);
        let mut record = Record::new();
        record.set(
            format!("read-{i}").as_bytes(),
            Some(&cigar),
            b"ACT",
            &[30; 3],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.push_aux(b"MM", Aux::String("A+m?,0;C+m?,0;T+g?,0;")).unwrap();
        record
            .push_aux(b"ML", Aux::ArrayU8((&probabilities[..]).into()))
            .unwrap();
        record.push_aux(b"MN", Aux::U32(3)).unwrap();
        record.push_aux(b"NM", Aux::U32(0)).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(bam_path.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    File::create(&fasta_path).unwrap().write_all(b">chr1\nACT\n").unwrap();
    File::create(root.join("reference.fa.fai"))
        .unwrap()
        .write_all(b"chr1\t3\t6\t3\t4\n")
        .unwrap();

    (bam_path, fasta_path)
}

fn write_dynamic_slot_fixture(root: &Path) -> (PathBuf, PathBuf) {
    write_dynamic_slot_fixture_with_probs(
        root,
        "dynamic-slots",
        &[[255, 255, 255]],
    )
}

fn write_reverse_motif_record(
    writer: &mut BamWriter,
    name: &str,
    sequence: &[u8],
    cigar: CigarString,
    mm_tag: &str,
    probability: u8,
    edit_distance: u32,
) {
    let qualities = vec![30; sequence.len()];
    let probabilities = [probability];
    let mut record = Record::new();
    record.set(name.as_bytes(), Some(&cigar), sequence, &qualities);
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    record.set_flags(16);
    record.push_aux(b"MM", Aux::String(mm_tag)).unwrap();
    record.push_aux(b"ML", Aux::ArrayU8((&probabilities[..]).into())).unwrap();
    record.push_aux(b"MN", Aux::U32(sequence.len() as u32)).unwrap();
    record.push_aux(b"NM", Aux::U32(edit_distance)).unwrap();
    writer.write(&record).unwrap();
}

fn write_combined_non_cpg_anchor_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let bam_path = root.join("combined-non-cpg-anchor.bam");
    let fasta_path = root.join("combined-non-cpg-anchor.fa");

    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 4);
    header.push_record(&sq);

    let mut writer =
        BamWriter::from_path(&bam_path, &header, Format::Bam).unwrap();
    let full_match = || CigarString(vec![Cigar::Match(4)]);

    // GATC is reverse-complement palindromic. Focusing offset 1 makes the
    // forward A anchor position 1 and the reverse T anchor position 2.
    for (name, probability) in
        [("modified", 255), ("canonical", 0), ("filtered", 128)]
    {
        write_reverse_motif_record(
            &mut writer,
            name,
            b"GATC",
            full_match(),
            "A+m?,0;",
            probability,
            0,
        );
    }
    write_reverse_motif_record(
        &mut writer,
        "no-call",
        b"GATC",
        full_match(),
        "C+m?,0;",
        255,
        0,
    );
    write_reverse_motif_record(
        &mut writer,
        "mismatch",
        b"GACC",
        full_match(),
        "C+m?,0;",
        255,
        1,
    );
    write_reverse_motif_record(
        &mut writer,
        "deletion",
        b"GAC",
        CigarString(vec![Cigar::Match(2), Cigar::Del(1), Cigar::Match(1)]),
        "C+m?,0;",
        255,
        1,
    );
    drop(writer);
    bam::index::build(bam_path.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    File::create(&fasta_path).unwrap().write_all(b">chr1\nGATC\n").unwrap();
    File::create(root.join("combined-non-cpg-anchor.fa.fai"))
        .unwrap()
        .write_all(b"chr1\t4\t6\t4\t5\n")
        .unwrap();

    (bam_path, fasta_path)
}

fn run_combined_non_cpg_anchor_pileup(
    root: &Path,
    bam_path: &Path,
    fasta_path: &Path,
    high_depth: bool,
    threads: usize,
) -> Vec<u8> {
    let mode = if high_depth { "high" } else { "normal" };
    let output_path =
        root.join(format!("combined-non-cpg-{mode}-{threads}.bed"));
    let threads = threads.to_string();
    let mut args = vec![
        "pileup",
        bam_path.to_str().unwrap(),
        output_path.to_str().unwrap(),
        "--ref",
        fasta_path.to_str().unwrap(),
        "--motif",
        "GATC",
        "1",
        "--combine-strands",
        "--modified-bases",
        "A:m",
        "--filter-threshold",
        "0.9",
        "--interval-size",
        "4",
        "--threads",
        &threads,
        "--suppress-progress",
    ];
    if high_depth {
        args.extend_from_slice(&["--high-depth", "--max-depth", "10"]);
    }

    run_modkit(&args)
        .with_context(|| {
            format!(
                "combined non-CpG pileup failed for mode {mode}, threads \
                 {threads}"
            )
        })
        .unwrap();
    std::fs::read(output_path).unwrap()
}

fn run_combined_slot_pileup(
    root: &Path,
    bam_path: &Path,
    fasta_path: &Path,
    bases: &[&str],
    high_depth: bool,
    threads: usize,
) -> Vec<u8> {
    let mode = if high_depth { "high" } else { "normal" };
    let output_path =
        root.join(format!("combine-{}-{mode}-{threads}.bed", bases.join("")));
    let threads = threads.to_string();
    let mut args = vec![
        "pileup",
        bam_path.to_str().unwrap(),
        output_path.to_str().unwrap(),
        "--ref",
        fasta_path.to_str().unwrap(),
        "--modified-bases",
    ];
    args.extend_from_slice(bases);
    args.extend_from_slice(&[
        "--combine-mods",
        "--no-filtering",
        "--interval-size",
        "1",
        "--threads",
        &threads,
        "--suppress-progress",
    ]);
    if high_depth {
        args.extend_from_slice(&["--high-depth", "--max-depth", "10"]);
    }

    run_modkit(&args)
        .with_context(|| {
            format!(
                "combine pileup failed for bases {}, mode {mode}, threads {}",
                bases.join(","),
                threads
            )
        })
        .unwrap();
    std::fs::read(output_path).unwrap()
}

fn assert_combined_slot_parity(bases: &[&str], expected: &str) {
    let temp_dir = tempfile::tempdir().unwrap();
    let root = temp_dir.path();
    let (bam_path, fasta_path) = write_dynamic_slot_fixture(root);
    let baseline =
        run_combined_slot_pileup(root, &bam_path, &fasta_path, bases, false, 1);
    assert_eq!(baseline, expected.as_bytes());

    for threads in [1, 2, 3, 8] {
        for high_depth in [false, true] {
            let observed = run_combined_slot_pileup(
                root,
                &bam_path,
                &fasta_path,
                bases,
                high_depth,
                threads,
            );
            assert_eq!(
                observed,
                baseline,
                "bases {}, high_depth {high_depth}, threads {threads}",
                bases.join(",")
            );
        }
    }
}

fn run_explicit_pair_slot_pileup(
    root: &Path,
    bam_path: &Path,
    fasta_path: &Path,
    high_depth: bool,
    threads: usize,
) -> Vec<u8> {
    let mode = if high_depth { "high" } else { "normal" };
    let output_path = root.join(format!("explicit-{mode}-{threads}.bed"));
    let threads = threads.to_string();
    let mut args = vec![
        "pileup",
        bam_path.to_str().unwrap(),
        output_path.to_str().unwrap(),
        "--ref",
        fasta_path.to_str().unwrap(),
        "--modified-bases",
        "A:m",
        "C:m",
        "--no-filtering",
        "--interval-size",
        "1",
        "--threads",
        &threads,
        "--suppress-progress",
    ];
    if high_depth {
        args.extend_from_slice(&["--high-depth", "--max-depth", "10"]);
    }

    run_modkit(&args)
        .with_context(|| {
            format!(
                "explicit pair pileup failed for mode {mode}, threads \
                 {threads}"
            )
        })
        .unwrap();
    std::fs::read(output_path).unwrap()
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
fn test_pileup_cpg_combined_cytosine_debug_regression() {
    let temp_file = std::env::temp_dir()
        .join("test_pileup_cpg_combined_cytosine_debug_regression.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--cpg",
        "--modified-bases",
        "C",
        "--combine-mods",
        "--no-filtering",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--threads",
        "1",
    ])
    .expect("valid overlapping CpG and internal C masks should not panic");

    let observed = read_bed_sites(temp_file.to_str().unwrap());
    let expected = read_bed_sites(
        "../tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
    assert!(!observed.is_empty());
    assert_eq!(observed, expected);
}

#[test]
fn test_pileup_combined_c_compact_slot_matches_high_depth() {
    assert_combined_slot_parity(
        &["C"],
        "chr1\t1\t2\tC\t1\t+\t1\t2\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0\n",
    );
}

#[test]
fn test_pileup_combined_act_compact_slots_match_high_depth() {
    assert_combined_slot_parity(
        &["A", "C", "T"],
        concat!(
            "chr1\t0\t1\tA\t1\t+\t0\t1\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0\n",
            "chr1\t1\t2\tC\t1\t+\t1\t2\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0\n",
            "chr1\t2\t3\tT\t1\t+\t2\t3\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0\n",
        ),
    );
}

#[test]
fn test_pileup_explicit_same_code_slots_are_keyed_by_base() {
    let temp_dir = tempfile::tempdir().unwrap();
    let root = temp_dir.path();
    let (bam_path, fasta_path) = write_dynamic_slot_fixture_with_probs(
        root,
        "explicit-pairs",
        &[[255, 255, 255], [0, 0, 0]],
    );
    let expected = concat!(
        "chr1\t0\t1\tm\t2\t+\t0\t1\t255,0,0\t2\t50.00\t1\t1\t0\t0\t0\t0\t0\n",
        "chr1\t1\t2\tm\t2\t+\t1\t2\t255,0,0\t2\t50.00\t1\t1\t0\t0\t0\t0\t0\n",
    );

    for threads in [1, 2, 3, 8] {
        for high_depth in [false, true] {
            let observed = run_explicit_pair_slot_pileup(
                root,
                &bam_path,
                &fasta_path,
                high_depth,
                threads,
            );
            assert_eq!(
                observed,
                expected.as_bytes(),
                "high_depth {high_depth}, threads {threads}"
            );
        }
    }
}

#[test]
fn test_pileup_combined_non_cpg_reverse_statuses_share_anchor() {
    let temp_dir = tempfile::tempdir().unwrap();
    let root = temp_dir.path();
    let (bam_path, fasta_path) = write_combined_non_cpg_anchor_fixture(root);
    let expected =
        b"chr1\t1\t2\tm\t2\t.\t1\t2\t255,0,0\t2\t50.00\t1\t1\t0\t1\t1\t1\t1\n";

    for threads in [1, 2] {
        for high_depth in [false, true] {
            let observed = run_combined_non_cpg_anchor_pileup(
                root,
                &bam_path,
                &fasta_path,
                high_depth,
                threads,
            );
            assert_eq!(
                observed, expected,
                "high_depth {high_depth}, threads {threads}"
            );
        }
    }
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
