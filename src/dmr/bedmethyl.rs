use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::{Itertools, MinMaxResult};
use nom::character::complete::{multispace1, none_of};
use nom::combinator::map_res;
use nom::multi::many1;
use nom::IResult;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::llr_model::AggregatedCounts;
use crate::genome_positions::StrandedPosition;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::parsing_utils::{
    consume_char, consume_digit, consume_float, consume_string,
    consume_string_from_list,
};
use crate::position_filter::Iv;
use crate::util::{Strand, StrandRule};

#[derive(new, Debug, PartialEq, Eq)]
pub struct BedMethylLine {
    pub chrom: String,
    pub interval: Iv,
    pub raw_mod_code: ModCodeRepr,
    pub strand: StrandRule,
    pub count_methylated: u64,
    pub valid_coverage: u64,
    pub count_canonical: u64,
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
    let (rest, count_methylated) = consume_digit(rest)?;
    let (_rest, count_canonical) = consume_digit(rest)?;

    let interval = Iv { start, stop, val: () };
    Ok((
        rest,
        BedMethylLine::new(
            chrom.to_string(),
            interval,
            raw_mod_code,
            strand,
            count_methylated,
            valid_coverage,
            count_canonical,
        ),
    ))
}

impl BedMethylLine {
    pub fn parse(line: &str) -> anyhow::Result<Self> {
        parse_bedmethyl_line(line).map(|(_, this)| this).map_err(|e| {
            anyhow!("failed to parse bedmethyl line {line}, {}", e.to_string())
        })
    }

    pub fn start(&self) -> u64 {
        self.interval.start
    }

    pub fn stop(&self) -> u64 {
        self.interval.stop
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

    /// panics if mod code isn't in `code_lookup`
    pub(crate) fn get_stranded_position(
        &self,
        code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
    ) -> StrandedPosition<DnaBase> {
        let strand = match self.strand {
            StrandRule::Positive | StrandRule::Both => Strand::Positive,
            StrandRule::Negative => Strand::Negative,
        };
        let dna_base = if strand == Strand::Negative {
            code_lookup.get(&self.raw_mod_code).unwrap().complement()
        } else {
            *code_lookup.get(&self.raw_mod_code).unwrap()
        };

        StrandedPosition { position: self.start(), strand, value: dna_base }
    }

    pub(crate) fn frac_modified(&self) -> f32 {
        self.count_methylated as f32 / self.valid_coverage as f32
    }
}

pub(super) fn aggregate_counts2(
    bm_lines: &[BedMethylLine],
    code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
) -> anyhow::Result<AggregatedCounts> {
    let lines = bm_lines.iter().collect::<Vec<&BedMethylLine>>();
    aggregate_counts(&lines, code_lookup)
}

pub(super) fn aggregate_counts(
    bm_lines: &[&BedMethylLine],
    code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
) -> anyhow::Result<AggregatedCounts> {
    assert_eq!(
        bm_lines.iter().map(|l| &l.chrom).collect::<HashSet<_>>().len(),
        1
    );
    // group by position because multiple mods are distributed over more than
    // one bedmethyl record
    let grouped_by_position: FxHashMap<
        StrandedPosition<DnaBase>,
        Vec<&BedMethylLine>,
    > = bm_lines.iter().fold(FxHashMap::default(), |mut acc, bm_line| {
        acc.entry(bm_line.get_stranded_position(code_lookup))
            .or_insert(Vec::new())
            .push(*bm_line);
        acc
    });
    let (counts_per_code, total) = grouped_by_position.into_iter().try_fold(
        (HashMap::new(), 0),
        |(mut acc, mut total_so_far), (_pos, grouped)| {
            let valid_covs = grouped
                .iter()
                .map(|bml| bml.valid_coverage)
                .collect::<FxHashSet<u64>>();
            let canonical_counts = grouped
                .iter()
                .map(|bml| bml.count_canonical)
                .collect::<FxHashSet<u64>>();
            let valid_coverage = grouped[0].valid_coverage as usize;
            if valid_covs.len() != 1 || canonical_counts.len() != 1 {
                let mut message = format!(
                    "invalid data found, should not have more than 1 score or \
                     number of canonical calls per position for a base. "
                );
                match grouped.iter().minmax_by(|a, b| a.start().cmp(&b.start()))
                {
                    MinMaxResult::NoElements => {}
                    MinMaxResult::MinMax(s, t) => message.push_str(&format!(
                        "starting at {}, ending at {}",
                        s.start(),
                        t.stop()
                    )),
                    MinMaxResult::OneElement(s) => {
                        message.push_str(&format!("starting at {}", s.start()))
                    }
                }
                bail!(message)
            }

            // check that the sum of canonical counts and
            // modified counts is equal to the valid coverage
            let mut check = grouped[0].count_canonical as usize;
            for x in &grouped {
                *acc.entry(x.raw_mod_code).or_insert(0) +=
                    x.count_methylated as usize;
                check += x.count_methylated as usize;
            }
            if check != valid_coverage {
                let mut message = format!(
                    "invalid data, valid coverage ({valid_coverage}) is not \
                     equal to the sum of canonical and modified counts \
                     ({check}), "
                );
                match grouped.iter().minmax_by(|a, b| a.start().cmp(&b.start()))
                {
                    MinMaxResult::NoElements => {}
                    MinMaxResult::MinMax(s, t) => message.push_str(&format!(
                        "chrom {} starting at {}, ending at {}",
                        s.chrom,
                        s.start(),
                        t.stop()
                    )),
                    MinMaxResult::OneElement(s) => message.push_str(&format!(
                        "chrom: {} starting at {}",
                        s.chrom,
                        s.start()
                    )),
                }
                bail!(message)
            }
            total_so_far += valid_coverage;
            Ok((acc, total_so_far))
        },
    )?;

    // todo don't need this match
    match AggregatedCounts::try_new(counts_per_code, total) {
        Ok(x) => Ok(x),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod bedmethylline_tests {
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    use crate::dmr::bedmethyl::{aggregate_counts2, BedMethylLine};
    use crate::genome_positions::GenomePositions;
    use crate::mod_base_code::{DnaBase, ModCodeRepr, MOD_CODE_TO_DNA_BASE};
    use crate::position_filter::Iv;
    use crate::util::StrandRule;

    #[test]
    #[rustfmt::skip]
    fn test_parse_bedmethyl_line() {
        for sep in ['\t', ' '] {
            let line = format!("chr20\t10034963\t10034964\tm,CG,0\t19\t-\t10034963\t10034964\t255,0,0\t19{sep}94.74{sep}18{sep}1{sep}0{sep}0{sep}1{sep}0{sep}2");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let start = 10034963;
            let stop = 10034964;
            let iv = Iv { start, stop, val: () };
            let expected = BedMethylLine::new(
                "chr20".to_string(),
                iv,
                'm'.into(),
                '-'.try_into().unwrap(),
                18,
                19,
                1

            );
            assert_eq!(bm_line, expected);
            let line = format!("chr20\t10034963\t10034964\tm\t19\t-\t10034963\t10034964\t255,0,0\t19{sep}94.74{sep}18{sep}1{sep}0{sep}0{sep}1{sep}0{sep}2");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            assert_eq!(bm_line, expected);

            let line = format!("oligo_1512_adapters\t9\t10\th\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0 ");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let expected = BedMethylLine::new(
                "oligo_1512_adapters".to_string(),
                Iv { start: 9, stop: 10, val: () },
                'h'.into(),
                '+'.try_into().unwrap(),
                2,
                4,
                1
            );
            assert_eq!(bm_line, expected);
        }
    }

    #[test]
    #[rustfmt::skip]
    fn test_parse_bedmethyl_line_chebi_code() {
        for sep in ['\t', ' '] {
            let line = format!("oligo_1512_adapters\t9\t10\t76792\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            let expected = BedMethylLine::new(
                "oligo_1512_adapters".to_string(),
                Iv { start: 9, stop: 10, val: () },
                76792u32.into(),
                '+'.try_into().unwrap(),
                2,
                4,
                1
            );
            assert_eq!(bm_line, expected);
            let line = format!("oligo_1512_adapters\t9\t10\t76792,CG,0\t4\t+\t9\t10\t255,0,0\t4{sep}50.00{sep}2{sep}1{sep}1{sep}0{sep}0{sep}2{sep}0");
            let bm_line = BedMethylLine::parse(&line).unwrap();
            assert_eq!(bm_line, expected);
        }
    }

    #[test]
    fn test_aggregate_counts() {
        let fh =
            File::open("tests/resources/modbam.modpileup_nofilt.methyl.bed")
                .unwrap();
        let dna_bases = vec![DnaBase::C];
        let all_contigs = HashSet::from(["oligo_1512_adapters".to_string()]);
        let mp = indicatif::MultiProgress::new();
        mp.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        let genome_positions = GenomePositions::new_from_sequences(
            &dna_bases,
            &Path::new("tests/resources/CGI_ladder_3.6kb_ref.fa").to_path_buf(),
            false,
            &all_contigs,
            &mp,
        )
        .unwrap();
        let positions = genome_positions
            .get_positions("oligo_1512_adapters", &(72..73), StrandRule::Both)
            .unwrap()
            .into_iter()
            .collect::<HashSet<_>>();
        let bedmethyl_lines = BufReader::new(fh)
            .lines()
            .map(|l| BedMethylLine::parse(&l.unwrap()).unwrap())
            .filter(|l| {
                positions
                    .contains(&l.get_stranded_position(&MOD_CODE_TO_DNA_BASE))
            })
            .collect::<Vec<BedMethylLine>>();
        let counts =
            aggregate_counts2(&bedmethyl_lines, &MOD_CODE_TO_DNA_BASE).unwrap();
        assert_eq!(&counts.string_counts(), "h:2,m:4");
        assert_eq!(counts.total, 6);
        let filtered_bm_lines = bedmethyl_lines
            .into_iter()
            .filter(|l| l.raw_mod_code == ModCodeRepr::Code('m'))
            .collect::<Vec<BedMethylLine>>();
        assert!(aggregate_counts2(&filtered_bm_lines, &MOD_CODE_TO_DNA_BASE)
            .is_err());
    }

    #[test]
    fn test_frac_modified() {
        let fp = "tests/resources/head_21839.bed";
        let reader = BufReader::new(File::open(fp).unwrap());
        let expected = vec![0f32, 0.07692308f32, 0f32];
        for (i, line) in reader.lines().map(|r| r.unwrap()).enumerate() {
            let record = BedMethylLine::parse(&line).unwrap();
            assert_eq!(record.frac_modified(), expected[i]);
        }
    }
}
