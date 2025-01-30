use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
};

use common::run_modkit;
use mod_kit::dmr::bedmethyl::BedMethylLine;

mod common;

#[test]
fn test_merge_bedmethyl_help() {
    run_modkit(&["bedmethyl", "--help"]).unwrap();
    run_modkit(&["bedmethyl", "merge", "--help"]).unwrap();
    run_modkit(&["bedmethyl", "tobigwig", "--help"]).unwrap();
}

#[test]
fn test_bedmethyl_merge() {
    let sizes = "chr20\t64444167\n";
    let sizes_fp = std::env::temp_dir().join("test_merge_bedmethyl_sizes.tsv");
    {
        let w = File::create(&sizes_fp).unwrap();
        let mut w = BufWriter::new(w);
        w.write(sizes.as_bytes()).unwrap();
        drop(w);
    };
    assert!(sizes_fp.exists());
    let bed_fp = "tests/resources/\
                  lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.\
                  gz";
    let out_bed = std::env::temp_dir().join("test_merge_bedmethyl.bed");

    run_modkit(&[
        "bedmethyl",
        "merge",
        bed_fp,
        bed_fp,
        "-g",
        sizes_fp.to_str().unwrap(),
        "-o",
        out_bed.to_str().unwrap(),
        "--force",
    ])
    .unwrap();
    assert!(out_bed.exists());
    // dbg!(&out_bed);

    let mut buff = vec![];
    let _ = rust_htslib::bgzf::Reader::from_path(bed_fp)
        .unwrap()
        .read_to_end(&mut buff);
    let input_records = buff
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .collect::<Vec<BedMethylLine>>();
    let merged_records = BufReader::new(File::open(out_bed).unwrap())
        .lines()
        .map(|line| BedMethylLine::parse(&line.unwrap()).unwrap())
        .collect::<Vec<BedMethylLine>>();
    assert_eq!(input_records.len(), merged_records.len());

    for (x, y) in input_records.into_iter().zip(merged_records) {
        assert_eq!(x.chrom, y.chrom);
        assert_eq!(x.start(), y.start());
        assert_eq!(x.raw_mod_code, y.raw_mod_code);
        assert_eq!(x.strand, y.strand);
        assert_eq!(x.count_methylated * 2, y.count_methylated);
        assert_eq!(x.valid_coverage * 2, y.valid_coverage);
        assert_eq!(x.count_canonical * 2, y.count_canonical);
        assert_eq!(x.count_other * 2, y.count_other);
        assert_eq!(x.count_delete * 2, y.count_delete);
        assert_eq!(x.count_fail * 2, y.count_fail);
        assert_eq!(x.count_diff * 2, y.count_diff);
        assert_eq!(x.count_nocall * 2, y.count_nocall);
    }
}
