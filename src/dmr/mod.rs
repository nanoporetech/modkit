use std::cmp::Ordering;
use std::fmt::Debug;

use anyhow::anyhow;
use derive_new::new;
use nom::character::complete::{multispace1, none_of, one_of};
use nom::multi::{many0, many1};
use nom::IResult;
use noodles::csi::index::{
    reference_sequence::bin::Chunk as IndexChunk, Index as CsiIndex,
};

use crate::parsing_utils::{
    consume_char, consume_char_from_list, consume_digit, consume_float,
    consume_string, consume_string_spaces,
};
use crate::position_filter::Iv;

mod model;
pub mod subcommand;

#[derive(new, Clone, Debug, Eq, PartialEq)]
pub(crate) struct DmrInterval {
    interval: Iv,
    chrom: String,
    name: String,
}

impl DmrInterval {
    fn parse_bed_line(line: &str) -> IResult<&str, Self> {
        let (rest, chrom) = consume_string(line)?;
        let (rest, start) = consume_digit(rest)?;
        let (rest, stop) = consume_digit(rest)?;

        let (rest, interval, name) = many0(one_of(" \t\r\n"))(rest)
            .and_then(|(rest, _)| consume_string_spaces(rest))
            .map(|(rest, name)| {
                let interval = Iv {
                    start,
                    stop,
                    val: (),
                };
                (rest, interval, name)
            })
            .unwrap_or_else(|_| {
                let interval = Iv {
                    start,
                    stop,
                    val: (),
                };
                let name = format!("{}:{}-{}", chrom, start, stop);
                (rest, interval, name)
            });

        Ok((
            rest,
            Self {
                interval,
                chrom,
                name,
            },
        ))
    }

    fn parse_str(line: &str) -> anyhow::Result<Self> {
        Self::parse_bed_line(line)
            .map(|(_, this)| this)
            .map_err(|e| anyhow!("{}", e.to_string()))
    }

    fn start(&self) -> u64 {
        self.interval.start
    }

    fn stop(&self) -> u64 {
        self.interval.stop
    }

    fn get_index_chunks(
        &self,
        index: &CsiIndex,
        chrom_id: usize,
    ) -> std::io::Result<Vec<IndexChunk>> {
        let start =
            noodles::core::Position::new((self.start() + 1) as usize).unwrap();
        let end =
            noodles::core::Position::new((self.stop() + 1) as usize).unwrap();
        let interval = noodles::core::region::Interval::from(start..=end);
        index.query(chrom_id, interval)
    }
}

impl PartialOrd for DmrInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl Ord for DmrInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.chrom.cmp(&other.chrom) {
            Ordering::Equal => self.interval.cmp(&other.interval),
            o @ _ => o,
        }
    }
}

#[derive(new, Debug, PartialEq, Eq)]
struct BedMethylLine {
    chrom: String,
    interval: Iv,
    raw_mod_code: char,
    // need to make this the object
    strand: char,
    count_methylated: u64,
    valid_coverage: u64,
}

fn parse_bedmethyl_line(l: &str) -> IResult<&str, BedMethylLine> {
    let (rest, chrom) = consume_string(l)?;
    let (rest, start) = consume_digit(rest)?;
    let (rest, stop) = consume_digit(rest)?;
    let (rest, _) = multispace1(rest)?;
    let (rest, raw_mod_code) = consume_char_from_list(rest, ",")?;
    let (rest, valid_coverage) = consume_digit(rest)?;
    let (rest, strand) = consume_char(rest)?;
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
    fn parse(line: &str) -> anyhow::Result<Self> {
        parse_bedmethyl_line(line)
            .map(|(_, this)| this)
            .map_err(|e| {
                anyhow!(
                    "failed to parse bedmethyl line {line}, {}",
                    e.to_string()
                )
            })
    }

    fn start(&self) -> u64 {
        self.interval.start
    }

    fn stop(&self) -> u64 {
        self.interval.stop
    }
}

#[cfg(test)]
mod dmr_mod_tests {
    use crate::dmr::{BedMethylLine, DmrInterval};
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
            let expected =
                BedMethylLine::new("chr20".to_string(), iv, 'm', '-', 18, 19);
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
                'h',
                '+',
                2,
                4,
            );
            assert_eq!(bm_line, expected);
        }
    }

    #[test]
    fn test_parse_rois() {
        let obs = DmrInterval::parse_str(
            "chr20\t279148\t279507\tCpG: 39 359\t39\t260\t21.7\t72.4\t0.83",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "CpG: 39 359".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507\tCpGby_any_other_name\t39\t260\t21.7\t72.4\t0.83").unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "CpGby_any_other_name".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507\t").unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "chr20:279148-279507".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507 ").unwrap();
        assert_eq!(obs, expected);
    }
}
