use anyhow::bail;
use derive_new::new;
use log::{debug, error, info};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{self, Read};
use rustc_hash::{FxHashMap, FxHashSet};
use std::ops::AddAssign;

use crate::errs::{MkError, MkResult};
use crate::find_motifs::motif_bed::OverlappingRegex;
use crate::mod_bam::{
    format_mm_ml_tag, BaseModProbs, CollapseMethod, EdgeFilter, ModBaseInfo,
    SeqPosBaseModProbs,
};
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{format_errors_table, get_query_name_string, get_ticker};

#[derive(new)]
pub(crate) struct OverlappingRegexOffset(OverlappingRegex, usize);

impl OverlappingRegexOffset {
    pub(crate) fn as_str(&self) -> &str {
        self.0.as_str()
    }

    fn find_iter<'a>(
        &'a self,
        seq: &'a str,
    ) -> impl Iterator<Item = usize> + 'a {
        self.0.find_iter(seq).map(|m| m.start().saturating_add(self.1))
    }
}

struct SequenceMotifs<'a> {
    complex: Vec<&'a OverlappingRegexOffset>,
    simple: FxHashSet<u8>,
}

impl<'a> SequenceMotifs<'a> {
    fn new(motifs: &'a [OverlappingRegexOffset]) -> Self {
        let (complex, simple) = motifs.iter().fold(
            (Vec::new(), Vec::new()),
            |(mut complex, mut simple), next| {
                if let Ok(base) = next
                    .as_str()
                    .parse::<char>()
                    .map_err(|_e| MkError::InvalidDnaBase)
                    .and_then(|c| DnaBase::parse(c))
                {
                    simple.push(base.char() as u8);
                } else {
                    complex.push(next)
                };
                (complex, simple)
            },
        );
        Self { complex, simple: simple.into_iter().collect() }
    }

    fn find_positions(&self, record: &bam::Record) -> FxHashSet<usize> {
        let seq = if record.is_reverse() {
            bio::alphabets::dna::revcomp(record.seq().as_bytes())
        } else {
            record.seq().as_bytes()
        };
        let mut positions = self
            .simple
            .par_iter()
            .flat_map(|&needle| {
                memchr::memchr_iter(needle, &seq).collect::<Vec<usize>>()
            })
            .collect::<FxHashSet<usize>>();
        if self.complex.is_empty() {
            positions
        } else {
            let seq = seq.into_iter().map(|b| b as char).collect::<String>();
            let complex_positions = self
                .complex
                .par_iter()
                .flat_map(|motif| motif.find_iter(&seq).collect::<Vec<usize>>())
                .collect::<FxHashSet<usize>>();

            positions.op_mut(complex_positions);
            positions
        }
    }
}

impl SeqPosBaseModProbs {
    fn filter_motif_positions(
        self,
        positions: &FxHashSet<usize>,
        discard: bool,
    ) -> Self {
        let probs = self
            .pos_to_base_mod_probs
            .into_iter()
            .filter(|(pos, _base_mod_probs)| {
                if discard {
                    !positions.contains(pos)
                } else {
                    positions.contains(pos)
                }
            })
            .collect::<FxHashMap<usize, BaseModProbs>>();

        Self::new(crate::mod_bam::SkipMode::Explicit, probs)
    }
}

fn adjust_mod_probs<'a>(
    mut record: bam::Record,
    methods: &[CollapseMethod],
    caller: Option<&MultipleThresholdModCaller>,
    edge_filter: Option<&EdgeFilter>,
    filter_only: bool,
    sequence_motifs: &Option<SequenceMotifs<'a>>,
    discard_motifs: bool,
) -> MkResult<bam::Record> {
    let mod_base_info = ModBaseInfo::new_from_record(&record)?;
    let mm_style = mod_base_info.mm_style;
    let ml_style = mod_base_info.ml_style;

    let mut mm_agg = String::new();
    let mut ml_agg = Vec::new();

    let record_name = get_query_name_string(&record)?;
    let (converters, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();

    let positions =
        sequence_motifs.as_ref().map(|ms| ms.find_positions(&record));

    for (base, strand, seq_pos_mod_probs) in mod_prob_iter {
        let converter = converters.get(&base).unwrap();
        // edge filter
        let trimmed_seq_pos_base_mod_probs = if let Some(edge_filter) =
            edge_filter
        {
            // remove the positions at the ends, also update to Ambiguous mode
            match seq_pos_mod_probs
                .edge_filter_positions(edge_filter, record.seq_len())
            {
                Some(x) => Some(x),
                None => {
                    debug!(
                        "all base mod positions for record {record_name} and \
                         canonical base {base} were filtered out"
                    );
                    None
                }
            }
        } else {
            Some(seq_pos_mod_probs)
        };
        if let Some(mut seq_pos_mod_probs) = trimmed_seq_pos_base_mod_probs {
            // collapse/convert
            for method in methods {
                seq_pos_mod_probs = seq_pos_mod_probs.into_collapsed(method);
            }
            // call mods
            // todo refactor
            match caller {
                Some(caller) => {
                    if filter_only {
                        seq_pos_mod_probs = caller
                            .filter_seq_pos_mod_probs(&base, seq_pos_mod_probs);
                    } else {
                        seq_pos_mod_probs = caller
                            .call_seq_pos_mod_probs(&base, seq_pos_mod_probs);
                    }
                }
                _ => {}
            }
            // motif filter
            if let Some(positions) = positions.as_ref() {
                seq_pos_mod_probs = seq_pos_mod_probs
                    .filter_motif_positions(positions, discard_motifs)
            }

            let (mm, mut ml) = format_mm_ml_tag(
                seq_pos_mod_probs,
                base,
                &converter.cumulative_counts,
                strand,
            );
            mm_agg.push_str(&mm);
            ml_agg.extend_from_slice(&mut ml);
        } else {
            continue;
        }
    }

    record.remove_aux(mm_style.as_bytes())?;
    record.remove_aux(ml_style.as_bytes())?;

    let mm = Aux::String(&mm_agg);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml_agg;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record.push_aux(mm_style.as_bytes(), mm)?;
    record.push_aux(ml_style.as_bytes(), ml)?;

    Ok(record)
}

pub(crate) fn adjust_modbam(
    reader: &mut bam::Reader,
    writer: &mut bam::Writer,
    collapse_methods: &[CollapseMethod],
    threshold_caller: Option<&MultipleThresholdModCaller>,
    edge_filter: Option<&EdgeFilter>,
    fail_fast: bool,
    motifs: &Option<Vec<OverlappingRegexOffset>>,
    discard_motifs: bool,
    verb: &'static str,
    suppress_progress: bool,
    filter_only: bool,
) -> anyhow::Result<()> {
    let spinner = get_ticker();
    if suppress_progress {
        spinner.set_draw_target(indicatif::ProgressDrawTarget::hidden())
    }
    spinner.set_message(verb);
    let mut total = 0usize;
    let mut error_counts = FxHashMap::<String, usize>::default();
    let sequence_motifs = motifs.as_ref().map(|x| SequenceMotifs::new(x));
    for (i, result) in reader
        .records()
        .map(|r| r.map_err(|e| MkError::HtsLibError(e)))
        .enumerate()
    {
        match result {
            Ok(record) => {
                match adjust_mod_probs(
                    record,
                    &collapse_methods,
                    threshold_caller,
                    edge_filter,
                    filter_only,
                    &sequence_motifs,
                    discard_motifs,
                ) {
                    Err(mk_error) => {
                        if fail_fast {
                            spinner.set_draw_target(
                                indicatif::ProgressDrawTarget::hidden(),
                            );
                            error!("encountered error, failing fast");
                            bail!("{mk_error}")
                        } else {
                            error_counts
                                .entry(mk_error.to_string())
                                .or_insert(0usize)
                                .add_assign(1usize);
                        }
                    }
                    Ok(record) => {
                        if let Err(err) = writer
                            .write(&record)
                            .map_err(|e| MkError::HtsLibError(e))
                        {
                            if fail_fast {
                                spinner.set_draw_target(
                                    indicatif::ProgressDrawTarget::hidden(),
                                );
                                error!("encountered error, failing fast");
                                bail!("{err}")
                            } else {
                                error_counts
                                    .entry(err.to_string())
                                    .or_insert(0usize)
                                    .add_assign(1usize);
                            }
                        } else {
                            spinner.inc(1);
                            total = i + 1;
                        }
                    }
                }
            }
            Err(mk_error) => {
                if fail_fast {
                    spinner.set_draw_target(
                        indicatif::ProgressDrawTarget::hidden(),
                    );
                    error!("encountered error, failing fast");
                    bail!("{mk_error}")
                } else {
                    error_counts
                        .entry(mk_error.to_string())
                        .or_insert(0usize)
                        .add_assign(1usize);
                }
            }
        }
    }
    spinner.finish_and_clear();

    info!("done, {} records processed", total,);

    if !error_counts.is_empty() {
        info!("error/skip counts:");
        let error_table = format_errors_table(&error_counts);
        info!("\n{error_table}");
    }

    Ok(())
}

#[cfg(test)]
mod adjust_tests {
    use std::collections::HashMap;

    use rust_htslib::bam::{self, Read};

    use crate::{
        find_motifs::motif_bed::RegexMotif,
        read_ids_to_base_mod_probs::ReadBaseModProfile,
    };

    use super::{adjust_mod_probs, OverlappingRegexOffset, SequenceMotifs};

    #[test]
    fn test_motif_filtering() {
        let bam_fp = "tests/resources/testing_all_context_calls.bam";
        let mut reader = bam::Reader::from_path(bam_fp).unwrap();

        let regex_motifs = vec![
            RegexMotif::parse_string("CG", 0).unwrap(),
            RegexMotif::parse_string("DRACH", 2).unwrap(),
        ]
        .into_iter()
        .map(|motif| {
            let offset = motif.forward_offset();
            OverlappingRegexOffset::new(motif.forward_pattern, offset)
        })
        .collect::<Vec<OverlappingRegexOffset>>();
        let sequence_motifs = Some(SequenceMotifs::new(&regex_motifs));

        let checks = HashMap::from([
            (
                'A' as u8,
                (
                    RegexMotif::parse_string("DRACH", 2)
                        .unwrap()
                        .forward_pattern,
                    0usize,
                ),
            ),
            (
                'C' as u8,
                (
                    RegexMotif::parse_string("CG", 0).unwrap().forward_pattern,
                    2usize,
                ),
            ),
        ]);

        let iter = reader
            .records()
            .map(|r| r.unwrap())
            .map(|record| {
                adjust_mod_probs(
                    record,
                    &[],
                    None,
                    None,
                    false,
                    &sequence_motifs,
                    false,
                )
                .unwrap()
            })
            .map(|record| {
                ReadBaseModProfile::from_record(&record, None, None, 5).unwrap()
            });

        let mut tested = false;
        for read in iter {
            for profile in read.profile {
                let kmer = format!("{}", profile.query_kmer);
                let (motif_for_base, match_position) =
                    checks.get(&profile.query_kmer.get_nt(2).unwrap()).unwrap();
                let matches = motif_for_base
                    .find_iter(&kmer)
                    .filter(|pos| pos.start() == *match_position)
                    .count()
                    > 0;
                assert!(matches, "should match {kmer}");
                tested = true;
            }
        }
        assert!(tested);

        let bam_fp = "tests/resources/testing_all_context_calls.bam";
        let mut reader = bam::Reader::from_path(bam_fp).unwrap();
        let iter = reader
            .records()
            .map(|r| r.unwrap())
            .map(|record| {
                adjust_mod_probs(
                    record,
                    &[],
                    None,
                    None,
                    false,
                    &sequence_motifs,
                    true,
                )
                .unwrap()
            })
            .map(|record| {
                ReadBaseModProfile::from_record(&record, None, None, 5).unwrap()
            });

        tested = false;
        for read in iter {
            for profile in read.profile {
                let kmer = format!("{}", profile.query_kmer);
                let (motif_for_base, match_position) =
                    checks.get(&profile.query_kmer.get_nt(2).unwrap()).unwrap();
                let matches = motif_for_base
                    .find_iter(&kmer)
                    .filter(|pos| pos.start() == *match_position)
                    .count()
                    > 0;
                assert!(!matches, "should not match {kmer}");
                tested = true;
            }
        }
        assert!(tested);
    }
}
