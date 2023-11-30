use crate::dmr::util::DmrInterval;
use anyhow::{anyhow, Context};
use derive_new::new;
use itertools::Itertools;
use nom::character::complete::{multispace1, none_of};
use nom::combinator::map_res;
use nom::multi::many1;
use nom::IResult;
use noodles::bgzf;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::BufRead;
use std::path::PathBuf;

use crate::mod_base_code::{DnaBase, ModCodeRepr, SUPPORTED_CODES};
use crate::parsing_utils::{
    consume_char, consume_digit, consume_float, consume_string,
    consume_string_from_list,
};
use crate::position_filter::Iv;
use crate::util::StrandRule;

#[derive(new, Debug, PartialEq, Eq)]
pub struct BedMethylLine {
    pub chrom: String,
    pub interval: Iv,
    pub raw_mod_code: ModCodeRepr,
    pub strand: StrandRule,
    pub count_methylated: u64,
    pub valid_coverage: u64,
}

fn parse_bedmethyl_line(l: &str) -> IResult<&str, BedMethylLine> {
    let comma_parser = |x| consume_string_from_list(x, ",");
    let mut parse_modcode = map_res(comma_parser, |raw| {
        ModCodeRepr::parse(raw)
            .with_context(|| format!("{raw} invalid mod code representation"))
    });
    let mut parse_strand = map_res(consume_char, |x| StrandRule::try_from(x));

    let (rest, chrom) = consume_string(l)?;
    let (rest, start) = consume_digit(rest)?;
    let (rest, stop) = consume_digit(rest)?;
    let (rest, _) = multispace1(rest)?;
    // let (rest, raw_mod_code) = consume_char_from_list(rest, ",")?;
    let (rest, raw_mod_code) = parse_modcode(rest)?;
    let (rest, valid_coverage) = consume_digit(rest)?;
    let (rest, strand) = parse_strand(rest)?;
    let (rest, _discard) = many1(consume_digit)(rest)?;
    let (rest, _discard_too) = many1(none_of(" \t"))(rest)?;
    let (rest, _score_again) = consume_digit(rest)?;
    let (rest, _pct_methyl) = consume_float(rest)?;
    let (_rest, count_methylated) = consume_digit(rest)?;

    let interval = Iv {
        start,
        stop,
        val: (),
    };
    Ok((
        rest,
        BedMethylLine::new(
            chrom.to_string(),
            interval,
            raw_mod_code,
            strand,
            count_methylated,
            valid_coverage,
        ),
    ))
}

impl BedMethylLine {
    pub fn parse(line: &str) -> anyhow::Result<Self> {
        parse_bedmethyl_line(line)
            .map(|(_, this)| this)
            .map_err(|e| {
                anyhow!(
                    "failed to parse bedmethyl line {line}, {}",
                    e.to_string()
                )
            })
    }

    pub fn start(&self) -> u64 {
        self.interval.start
    }

    pub fn stop(&self) -> u64 {
        self.interval.stop
    }

    pub fn check_mod_code_supported(&self) -> bool {
        SUPPORTED_CODES.contains(&self.raw_mod_code)
    }

    pub fn check_base(
        &self,
        dna_base: DnaBase,
        additional_mappings: Option<&HashMap<ModCodeRepr, DnaBase>>,
    ) -> bool {
        if self.raw_mod_code.check_base(dna_base) {
            true
        } else {
            if let Some(mappings) = additional_mappings {
                mappings
                    .get(&self.raw_mod_code)
                    .map(|b| *b == dna_base)
                    .unwrap_or(false)
            } else {
                false
            }
        }
    }
}

pub(super) fn load_regions_from_bedmethyl(
    bedmethyl_fp: &PathBuf,
) -> anyhow::Result<Vec<DmrInterval>> {
    let reader = File::open(bedmethyl_fp).map(bgzf::Reader::new)?;
    let regions = reader
        .lines()
        .filter_map(|l| l.ok())
        .map(|l| {
            BedMethylLine::parse(&l).map(|bm| {
                let interval = bm.interval;
                let name = format!(
                    "{}:{}-{}",
                    &bm.chrom, interval.start, interval.stop
                );
                let chrom = bm.chrom;
                DmrInterval::new(interval, chrom, name)
            })
        })
        .collect::<anyhow::Result<BTreeSet<DmrInterval>>>()?
        .into_iter()
        .collect_vec();
    Ok(regions)
}

#[cfg(test)]
mod bedmethylline_tests {
    use crate::dmr::bedmethyl::BedMethylLine;
    use crate::position_filter::Iv;

    #[test]
    fn test_parse_bedmethyl_line() {
        for sep in ['\t', ' '] {
            let line = format!("chr20\t10034963\t10034964\tm,CG,0\t19\t-\t10034963\t10034964\t255,0,0\t19{sep}94.74{sep}18{sep}1{sep}0{sep}0{sep}1{sep}0{sep}2");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let start = 10034963;
            let stop = 10034964;
            let iv = Iv {
                start,
                stop,
                val: (),
            };
            let expected = BedMethylLine::new(
                "chr20".to_string(),
                iv,
                'm'.into(),
                '-'.try_into().unwrap(),
                18,
                19,
            );
            assert_eq!(bm_line, expected);
            let line = format!("chr20\t10034963\t10034964\tm\t19\t-\t10034963\t10034964\t255,0,0\t19{sep}94.74{sep}18{sep}1{sep}0{sep}0{sep}1{sep}0{sep}2");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            assert_eq!(bm_line, expected);

            let line = format!("oligo_1512_adapters\t9\t10\th\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0 ");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let expected = BedMethylLine::new(
                "oligo_1512_adapters".to_string(),
                Iv {
                    start: 9,
                    stop: 10,
                    val: (),
                },
                'h'.into(),
                '+'.try_into().unwrap(),
                2,
                4,
            );
            assert_eq!(bm_line, expected);
        }
    }

    #[test]
    fn test_parse_bedmethyl_line_chebi_code() {
        for sep in ['\t', ' '] {
            let line = format!("oligo_1512_adapters\t9\t10\t76792\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0 ");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let expected = BedMethylLine::new(
                "oligo_1512_adapters".to_string(),
                Iv {
                    start: 9,
                    stop: 10,
                    val: (),
                },
                76792u32.into(),
                '+'.try_into().unwrap(),
                2,
                4,
            );
            assert_eq!(bm_line, expected);
            let line = format!("oligo_1512_adapters\t9\t10\t76792,CG,0\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0 ");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            assert_eq!(bm_line, expected);
        }
    }
}
