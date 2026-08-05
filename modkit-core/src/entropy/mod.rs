use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::{Itertools, MinMaxResult};
use log::{debug, info};
use nom::character::complete::multispace1;
use nom::IResult;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;

use crate::entropy::methylation_entropy::{
    calc_me_entropy, EntropyPattern, EntropySymbol,
};
use crate::errs::{MkError, MkResult};
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::motifs::motif_bed::RegexMotif;
use crate::read_ids_to_base_mod_probs::{PositionModCalls, ReadBaseModProfile};
use crate::reads_sampler::sampling_schedule::ReferenceSequencesLookup;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::percentile_linear_interp;
use crate::util::{record_is_not_primary, ReferenceRecord, Strand};

mod methylation_entropy;
pub mod subcommand;
mod writers;

type BaseAndPosition = (DnaBase, u64);

fn half_open_interval<I>(positions: I) -> Range<u64>
where
    I: IntoIterator<Item = u64>,
{
    match positions.into_iter().minmax() {
        MinMaxResult::MinMax(start, end) => {
            start..end.checked_add(1).expect("reference interval overflow")
        }
        MinMaxResult::OneElement(position) => {
            position
                ..position
                    .checked_add(1)
                    .expect("reference interval overflow")
        }
        MinMaxResult::NoElements => {
            unreachable!("cannot build an interval without positions")
        }
    }
}

#[derive(Debug)]
pub(super) enum GenomeWindow {
    CombineStrands {
        interval: Range<u64>,
        neg_to_pos_positions: FxHashMap<BaseAndPosition, BaseAndPosition>,
        read_patterns: Vec<Vec<BaseModCall>>,
        position_valid_coverages: Vec<u32>,
    },
    Stranded {
        // todo instead of having pos/neg for everything, make one struct and
        // have  an optional for all of it
        pos_interval: Option<Range<u64>>,
        neg_interval: Option<Range<u64>>,
        pos_positions: Option<Vec<BaseAndPosition>>,
        neg_positions: Option<Vec<BaseAndPosition>>,
        pos_read_patterns: Vec<Vec<BaseModCall>>,
        neg_read_patterns: Vec<Vec<BaseModCall>>,
        pos_position_valid_coverages: Vec<u32>,
        neg_position_valid_coverages: Vec<u32>,
    },
}

impl GenomeWindow {
    fn new_combine_strands(
        interval: Range<u64>,
        num_positions: usize,
        neg_to_pos_positions: FxHashMap<BaseAndPosition, BaseAndPosition>,
    ) -> Self {
        let position_valid_coverages = vec![0u32; num_positions];
        Self::CombineStrands {
            interval,
            neg_to_pos_positions,
            read_patterns: Vec::new(),
            position_valid_coverages,
        }
    }

    fn new_stranded(
        pos_positions: Option<Vec<BaseAndPosition>>,
        neg_positions: Option<Vec<BaseAndPosition>>,
        num_positions: usize,
    ) -> Self {
        let pos_interval = pos_positions.as_ref().map(|positions| {
            half_open_interval(positions.iter().map(|(_, position)| *position))
        });
        let neg_interval = neg_positions.as_ref().map(|positions| {
            half_open_interval(positions.iter().map(|(_, position)| *position))
        });

        #[cfg(debug_assertions)]
        let check = |positions: Option<&Vec<BaseAndPosition>>| {
            if let Some(ps) = positions {
                ps.iter().skip(1).fold(ps[0].1, |last, (_, next)| {
                    assert!(last < *next, "needs to be sorted");
                    *next
                });
            }
        };

        #[cfg(debug_assertions)]
        check(pos_positions.as_ref());
        #[cfg(debug_assertions)]
        check(neg_positions.as_ref());

        let pos_position_valid_coverages = vec![0u32; num_positions];
        let neg_position_valid_coverages = vec![0u32; num_positions];
        // debug!(
        //     "interval {pos_interval:?}, {neg_interval:?} \n\t> pos: \
        //      {pos_positions:?} neg {neg_positions:?}"
        // );
        Self::Stranded {
            pos_interval,
            neg_interval,
            pos_positions,
            neg_positions,
            pos_read_patterns: Vec::new(),
            neg_read_patterns: Vec::new(),
            pos_position_valid_coverages,
            neg_position_valid_coverages,
        }
    }

    #[inline]
    fn inc_coverage(&mut self, pos: usize, strand: &Strand) {
        match self {
            Self::CombineStrands { position_valid_coverages, .. } => {
                assert!(
                    pos < position_valid_coverages.len(),
                    "pos is larger than the window size?"
                );
                position_valid_coverages[pos] += 1;
            }
            Self::Stranded {
                pos_position_valid_coverages,
                neg_position_valid_coverages,
                ..
            } => match strand {
                Strand::Positive => {
                    assert!(
                        pos < pos_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    pos_position_valid_coverages[pos] += 1;
                }
                Strand::Negative => {
                    assert!(
                        pos < neg_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    neg_position_valid_coverages[pos] += 1;
                }
            },
        };
    }

    fn add_pattern(&mut self, strand: &Strand, pattern: Vec<BaseModCall>) {
        match self {
            Self::Stranded { pos_read_patterns, neg_read_patterns, .. } => {
                match strand {
                    Strand::Positive => pos_read_patterns.push(pattern),
                    Strand::Negative => neg_read_patterns.push(pattern),
                }
            }
            Self::CombineStrands { read_patterns, .. } => {
                read_patterns.push(pattern);
            }
        }
    }

    fn leftmost(&self) -> u64 {
        match (self.start(&Strand::Positive), self.start(&Strand::Negative)) {
            (Some(x), Some(y)) => std::cmp::min(x, y),
            (Some(x), None) => x,
            (None, Some(x)) => x,
            _ => unreachable!(
                "should always have either a positive or negative interval!"
            ),
        }
    }

    fn exclusive_end(&self) -> u64 {
        match (self.end(&Strand::Positive), self.end(&Strand::Negative)) {
            (Some(x), Some(y)) => std::cmp::max(x, y),
            (Some(x), None) => x,
            (None, Some(x)) => x,
            _ => unreachable!(
                "should always have either a positive or negative interval!"
            ),
        }
    }

    fn start(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.start),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.start),
                Strand::Negative => neg_interval.as_ref().map(|x| x.start),
            },
        }
    }

    fn end(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.end),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.end),
                Strand::Negative => neg_interval.as_ref().map(|x| x.end),
            },
        }
    }

    fn add_read_to_patterns(
        &mut self,
        ref_pos_to_basemod_call: &FxHashMap<BaseAndPosition, BaseModCall>,
        reference_start: i64,
        reference_end: i64,
        strand: Strand,
        max_filtered_positions: usize,
    ) {
        // check that the read fully covers the interval
        let reference_start = if reference_start >= 0 {
            Some(reference_start as u64)
        } else {
            None
        };
        let reference_end = if reference_start
            .map(|x| reference_end > x as i64)
            .unwrap_or(false)
        {
            Some(reference_end as u64)
        } else {
            None
        };

        let overlaps = reference_start
            .and_then(|s| reference_end.map(|t| (s, t)))
            .map(|(s, t)| match (self.start(&strand), self.end(&strand)) {
                (Some(wind_start), Some(wind_end)) => {
                    s <= wind_start && t >= wind_end
                }
                _ => false,
            })
            // .map(|(s, t)| s <= self.start() && t >= self.end())
            .unwrap_or(false);
        if !overlaps {
            return;
        }

        let pattern: Vec<BaseModCall> = match strand {
            Strand::Positive => match &self {
                Self::Stranded { pos_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .values()
                        .map(|p| {
                            let call = ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (p, call)
                        })
                        .sorted_by(|((_, a), _), ((_, b), _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
            Strand::Negative => match &self {
                Self::Stranded { neg_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .iter()
                        .map(|(neg_position, positive_position)| {
                            let call = ref_pos_to_basemod_call
                                .get(neg_position)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (positive_position, call)
                        })
                        .sorted_by(|((_, a), _), ((_, b), _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
        };

        if pattern.iter().filter(|&bmc| bmc == &BaseModCall::Filtered).count()
            > max_filtered_positions
        {
            // skip when too many filtered positions
            return;
        }

        for (i, call) in pattern.iter().enumerate() {
            match call {
                BaseModCall::Filtered => {}
                _ => self.inc_coverage(i, &strand),
            }
        }
        self.add_pattern(&strand, pattern);
    }

    fn get_mod_code_lookup(&self) -> FxHashMap<ModCodeRepr, EntropySymbol> {
        // looks complicated, but it just iterates over either the positive and
        // negative read patterns or the positive-combined read patterns
        let read_patterns: Box<dyn Iterator<Item = &Vec<BaseModCall>>> =
            match self {
                Self::Stranded {
                    pos_read_patterns, neg_read_patterns, ..
                } => {
                    Box::new(pos_read_patterns.iter().chain(neg_read_patterns))
                }
                Self::CombineStrands { read_patterns, .. } => {
                    Box::new(read_patterns.iter())
                }
            };

        read_patterns
            .flat_map(|pattern| {
                pattern.iter().filter_map(|call| match call {
                    BaseModCall::Modified(_, code) => Some(*code),
                    _ => None,
                })
            })
            .collect::<BTreeSet<ModCodeRepr>>()
            .into_iter()
            .enumerate()
            .map(|(id, code)| {
                // save 0 for canonical
                let id = id
                    .checked_add(1)
                    .and_then(|id| u32::try_from(id).ok())
                    .filter(|id| *id < u32::MAX)
                    .expect(
                        "modification state count must fit below the reserved \
                         wildcard symbol",
                    );
                (code, EntropySymbol::called(id))
            })
            .collect::<FxHashMap<ModCodeRepr, EntropySymbol>>()
    }

    fn encode_patterns(
        &self,
        chrom_id: u32,
        strand: Strand,
        patterns: &Vec<Vec<BaseModCall>>,
        mod_code_lookup: &FxHashMap<ModCodeRepr, EntropySymbol>,
        position_valid_coverages: &[u32],
        min_coverage: u32,
    ) -> MkResult<Vec<EntropyPattern>> {
        // todo remove these checks after testing
        assert!(
            self.start(&strand).is_some(),
            "start should be Some when encoding pattern for strand \
             {strand:?}, {patterns:?}"
        );
        assert!(
            self.end(&strand).is_some(),
            "end should be Some when encoding pattern for strand {strand:?}, \
             {patterns:?}"
        );

        if position_valid_coverages.iter().all(|x| *x >= min_coverage) {
            let encoded = patterns
                .iter()
                .map(|pat| {
                    let pattern = pat
                        .iter()
                        .map(|call| match call {
                            BaseModCall::Canonical(_) => {
                                EntropySymbol::CANONICAL
                            }
                            BaseModCall::Modified(_, code) => {
                                *mod_code_lookup.get(code).unwrap()
                            }
                            BaseModCall::Filtered => EntropySymbol::FILTERED,
                        })
                        .collect::<Vec<EntropySymbol>>()
                        .into_boxed_slice();
                    // todo remove after testing
                    assert_eq!(
                        pattern.len(),
                        position_valid_coverages.len(),
                        "pattern {pattern:?} is the wrong size? \
                         {position_valid_coverages:?}"
                    );
                    pattern
                })
                .collect();
            Ok(encoded)
        } else {
            let zero_coverage =
                position_valid_coverages.iter().all(|&cov| cov == 0);
            if zero_coverage {
                return Err(MkError::EntropyZeroCoverage {
                    chrom_id,
                    start: self.start(&strand).unwrap(),
                    end: self.end(&strand).unwrap(),
                });
            } else {
                let err = MkError::EntropyInsufficientCoverage {
                    chrom_id,
                    start: self.start(&strand).unwrap(),
                    end: self.end(&strand).unwrap(),
                };
                return Err(err);
            }
        }
    }

    fn into_entropy(
        &self,
        chrom_id: u32,
        min_valid_coverage: u32,
    ) -> WindowEntropy {
        let window_size = self.size();
        let constant = 1f32 / window_size as f32; // todo make this configurable

        let mod_code_lookup = self.get_mod_code_lookup();
        let positive_encoded_patterns = match &self {
            Self::CombineStrands {
                read_patterns,
                position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Positive,
                read_patterns,
                &mod_code_lookup,
                position_valid_coverages,
                min_valid_coverage,
            )),
            Self::Stranded {
                pos_interval: Some(_),
                pos_read_patterns,
                pos_position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Positive,
                pos_read_patterns,
                &mod_code_lookup,
                &pos_position_valid_coverages,
                min_valid_coverage,
            )),
            _ => None,
        };
        let negative_patterns = match &self {
            Self::Stranded {
                neg_interval: Some(_),
                neg_read_patterns,
                neg_position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Negative,
                neg_read_patterns,
                &mod_code_lookup,
                neg_position_valid_coverages,
                min_valid_coverage,
            )),
            _ => None,
        };
        // left for debugging
        // debug!(
        //     "{}:{}-{} (+), {:?}",
        //     chrom,
        //     self.leftmost(),
        //     self.rightmost(),
        //     &positive_encoded_patterns
        // );
        // if let Some(nps) = negative_patterns.as_ref() {
        //     debug!(
        //         "{}:{}-{} (-), {:?}",
        //         chrom,
        //         self.leftmost(),
        //         self.rightmost(),
        //         &nps
        //     );
        // }

        // TODO: make sure there is a proper entropy test
        #[cfg(debug_assertions)]
        {
            if let Some(Ok(patterns)) = positive_encoded_patterns.as_ref() {
                debug_assert!(
                    patterns.iter().all(|x| x.len() == window_size),
                    "patterns are the wrong size {positive_encoded_patterns:?}"
                );
            }
            if let Some(Ok(neg_patterns)) = negative_patterns.as_ref() {
                debug_assert!(neg_patterns
                    .iter()
                    .all(|x| x.len() == window_size));
            }
        }

        let pos_me_entropy = positive_encoded_patterns.map(|maybe_patterns| {
            maybe_patterns.map(|patterns| {
                let me_entropy =
                    calc_me_entropy(&patterns, window_size, constant);
                let num_reads = patterns.len();
                let interval = self.start(&Strand::Positive).unwrap()
                    ..self.end(&Strand::Positive).unwrap();
                MethylationEntropy::new(me_entropy, num_reads, interval)
            })
        });

        let neg_me_entropy = negative_patterns.map(|maybe_patterns| {
            maybe_patterns.map(|patterns| {
                let me_entropy =
                    calc_me_entropy(&patterns, window_size, constant);
                let num_reads = patterns.len();
                let interval = self.start(&Strand::Negative).unwrap()
                    ..self.end(&Strand::Negative).unwrap();
                MethylationEntropy::new(me_entropy, num_reads, interval)
            })
        });

        WindowEntropy::new(chrom_id, pos_me_entropy, neg_me_entropy)
    }

    #[inline]
    fn size(&self) -> usize {
        match self {
            Self::Stranded { pos_position_valid_coverages, .. } => {
                pos_position_valid_coverages.len()
            }
            Self::CombineStrands { position_valid_coverages, .. } => {
                position_valid_coverages.len()
            }
        }
    }
}

pub(super) struct GenomeWindows {
    chrom_id: u32,
    entropy_windows: Vec<GenomeWindow>,
    region_name: Option<String>,
}

pub(super) enum EntropyCalculation {
    Windows(Vec<WindowEntropy>),
    Region(RegionEntropy),
}

impl GenomeWindows {
    fn new(
        chrom_id: u32,
        entropy_windows: Vec<GenomeWindow>,
        region_name: Option<String>,
    ) -> Self {
        assert!(!entropy_windows.is_empty());
        Self { chrom_id, entropy_windows, region_name }
    }

    fn get_range(&self) -> Range<u64> {
        let start = self
            .entropy_windows
            .iter()
            .map(GenomeWindow::leftmost)
            .min()
            .expect("self.entropy_windows should not be empty");
        let end = self
            .entropy_windows
            .iter()
            .map(GenomeWindow::exclusive_end)
            .max()
            .expect("self.entropy_windows should not be empty");
        start..end
    }

    fn get_fetch_definition(&self) -> FetchDefinition<'_> {
        let range = self.get_range();
        let start = range.start as i64;
        let end = range.end as i64;
        let chrom_id = self.chrom_id;
        FetchDefinition::Region(chrom_id as i32, start, end)
    }

    fn into_entropy_calculation(
        self,
        chrom_id: u32,
        min_coverage: u32,
    ) -> EntropyCalculation {
        // to appease the bC we have to get the interval
        // here, but it's only used if we're summarizing a region
        let interval = self.get_range();
        let window_entropies = self
            .entropy_windows
            .par_iter()
            .map(|ew| ew.into_entropy(chrom_id, min_coverage))
            .collect::<Vec<_>>();
        let chrom_id = self.chrom_id;
        if let Some(region_name) = self.region_name {
            let mut pos_entropies = Vec::with_capacity(window_entropies.len());
            let mut pos_num_reads = Vec::with_capacity(window_entropies.len());
            let mut pos_num_fails = 0usize;
            let mut neg_entropies = Vec::with_capacity(window_entropies.len());
            let mut neg_num_reads = Vec::with_capacity(window_entropies.len());
            let mut neg_num_fails = 0usize;

            for window_entropy in window_entropies.iter() {
                match window_entropy.pos_me_entropy.as_ref() {
                    Some(Ok(me_entropy)) => {
                        pos_entropies.push(me_entropy.me_entropy);
                        pos_num_reads.push(me_entropy.num_reads);
                    }
                    Some(Err(_e)) => {
                        pos_num_fails += 1;
                    }
                    None => {}
                }
                match window_entropy.neg_me_entropy.as_ref() {
                    Some(Ok(me_entropy)) => {
                        neg_entropies.push(me_entropy.me_entropy);
                        neg_num_reads.push(me_entropy.num_reads);
                    }
                    Some(Err(_e)) => {
                        neg_num_fails += 1;
                    }
                    // this means it was combine strands
                    None => {}
                }
            }

            // todo make sure the semantics here are what I want,
            //  should pos_entropy_stats be an Option?
            let pos_entropy_stats = DescriptiveStats::new(
                &pos_entropies,
                &pos_num_reads,
                pos_num_fails,
                chrom_id,
                &interval,
            );
            // if neg_entropies is empty and there are no fails, we never saw
            // any negative strand me entropies
            let neg_entropy_stats = if neg_entropies.is_empty()
                && neg_num_fails == 0
            {
                assert!(
                    neg_num_reads.is_empty(),
                    "neg num reads and window entropies should both be empty"
                );
                None
            } else {
                // this will fail correctly if there are neg_entropies is empty
                // but there are fails
                Some(DescriptiveStats::new(
                    &neg_entropies,
                    &neg_num_reads,
                    neg_num_fails,
                    chrom_id,
                    &interval,
                ))
            };

            let region_entropy = RegionEntropy::new(
                chrom_id,
                interval,
                pos_entropy_stats,
                neg_entropy_stats,
                region_name,
                window_entropies,
            );
            EntropyCalculation::Region(region_entropy)
        } else {
            EntropyCalculation::Windows(window_entropies)
        }
    }
}

#[derive(Debug, new)]
struct MotifHit {
    pos: u64,
    neg_position: Option<u64>,
    strand: Strand,
    base: DnaBase,
    motif_idx: usize,
}

#[derive(Debug)]
struct ReferenceSearchSpace {
    record: ReferenceRecord,
    sequence: Vec<char>,
    owner: Range<usize>,
}

#[derive(Debug)]
struct ScannedReference {
    record: ReferenceRecord,
    owner: Range<usize>,
    motif_hits: Vec<MotifHit>,
}

struct SlidingWindows {
    motifs: Vec<RegexMotif>,
    work_queue: VecDeque<ReferenceSearchSpace>,
    region_names: VecDeque<String>,
    window_size: usize,
    num_positions: usize,
    batch_size: usize,
    curr_position: usize,
    curr_hit_cursor: usize,
    curr_contig: ReferenceRecord,
    curr_owner: Range<usize>,
    curr_motif_hits: Vec<MotifHit>,
    curr_region_name: Option<String>,
    combine_strands: bool,
    done: bool,
}

impl SlidingWindows {
    const START_SEARCH_CHUNK_SIZE: usize = 10_000;

    fn motif_search_context(motifs: &[RegexMotif]) -> usize {
        motifs
            .iter()
            .map(|motif| motif.length().saturating_sub(1))
            .max()
            .unwrap_or(0)
    }

    fn find_motif_hits_in_owner(
        seq: &[char],
        motifs: &[RegexMotif],
        owner: Range<usize>,
        reference_start: u64,
        motif_search_context: usize,
    ) -> Vec<MotifHit> {
        assert!(owner.start <= owner.end);
        assert!(owner.end <= seq.len());
        if owner.is_empty() {
            return Vec::new();
        }

        let search_start = owner.start.saturating_sub(motif_search_context);
        let search_end = owner
            .end
            .saturating_add(motif_search_context)
            .min(seq.len());
        let subseq = seq[search_start..search_end].iter().collect::<String>();

        motifs
            .iter()
            .enumerate()
            .flat_map(|(motif_idx, motif)| {
                motif
                    .find_hits(&subseq)
                    .into_iter()
                    .filter_map(|(search_position, strand)| {
                        let local_position =
                            search_start.checked_add(search_position)?;
                        if !owner.contains(&local_position) {
                            return None;
                        }
                        let position = reference_start
                            .checked_add(local_position as u64)
                            .expect("reference position overflow");
                        let dna_base = DnaBase::parse(seq[local_position])
                            .expect("motif anchor must be a DNA base");
                        let base = if strand == Strand::Negative {
                            dna_base.complement()
                        } else {
                            dna_base
                        };
                        let neg_position = motif
                            .motif_info
                            .negative_strand_position(position as u32)
                            .map(|position| position as u64);
                        Some(MotifHit::new(
                            position,
                            neg_position,
                            strand,
                            base,
                            motif_idx,
                        ))
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn sort_and_dedup_motif_hits(
        mut motif_hits: Vec<MotifHit>,
        motifs: &[RegexMotif],
        combine_strands: bool,
        reference_name: &str,
    ) -> anyhow::Result<Vec<MotifHit>> {
        motif_hits.sort_by_key(|hit| {
            (hit.pos, hit.strand, hit.base, hit.motif_idx)
        });
        let mut deduped: Vec<MotifHit> = Vec::with_capacity(motif_hits.len());
        for hit in motif_hits {
            if let Some(previous) = deduped.last() {
                let same_anchor = previous.pos == hit.pos
                    && previous.strand == hit.strand
                    && previous.base == hit.base;
                if same_anchor {
                    if combine_strands
                        && hit.strand == Strand::Positive
                        && previous.neg_position != hit.neg_position
                    {
                        bail!(
                            "conflicting combined-strand motif partners at \
                             {reference_name}:{} for {} anchor: {:?} from \
                             motif {} versus {:?} from motif {}",
                            hit.pos,
                            hit.base,
                            previous.neg_position,
                            motifs[previous.motif_idx],
                            hit.neg_position,
                            motifs[hit.motif_idx],
                        )
                    }
                    continue;
                }
            }
            deduped.push(hit);
        }

        if combine_strands {
            let mut partner_to_positive = BTreeMap::new();
            for hit in deduped
                .iter()
                .filter(|hit| hit.strand == Strand::Positive)
            {
                let Some(negative_position) = hit.neg_position else {
                    continue;
                };
                let negative_partner = (hit.base, negative_position);
                if let Some(previous) =
                    partner_to_positive.insert(negative_partner, hit)
                {
                    if previous.pos != hit.pos {
                        bail!(
                            "conflicting combined-strand motif anchors at \
                             {reference_name}: negative {} partner {} maps \
                             to positive anchors {} from motif {} and {} \
                             from motif {}",
                            hit.base,
                            negative_position,
                            previous.pos,
                            motifs[previous.motif_idx],
                            hit.pos,
                            motifs[hit.motif_idx],
                        )
                    }
                }
            }
        }
        Ok(deduped)
    }

    fn scan_motif_hits_with_chunk_size(
        seq: &[char],
        motifs: &[RegexMotif],
        owner: Range<usize>,
        reference_start: u64,
        reference_name: &str,
        combine_strands: bool,
        chunk_size: usize,
    ) -> anyhow::Result<Vec<MotifHit>> {
        assert!(chunk_size > 0, "motif search chunk size must be positive");
        assert!(owner.start <= owner.end);
        assert!(owner.end <= seq.len());
        let motif_search_context = Self::motif_search_context(motifs);
        let owner_chunks = (owner.start..owner.end)
            .step_by(chunk_size)
            .map(|start| start..start.saturating_add(chunk_size).min(owner.end))
            .collect::<Vec<_>>();
        let motif_hits = owner_chunks
            .into_par_iter()
            .flat_map(|owner_chunk| {
                Self::find_motif_hits_in_owner(
                    seq,
                    motifs,
                    owner_chunk,
                    reference_start,
                    motif_search_context,
                )
                .into_par_iter()
            })
            .collect::<Vec<_>>();
        Self::sort_and_dedup_motif_hits(
            motif_hits,
            motifs,
            combine_strands,
            reference_name,
        )
    }

    fn scan_reference(
        search_space: ReferenceSearchSpace,
        motifs: &[RegexMotif],
        combine_strands: bool,
    ) -> anyhow::Result<ScannedReference> {
        let ReferenceSearchSpace { record, sequence, owner } = search_space;
        let motif_hits = Self::scan_motif_hits_with_chunk_size(
            &sequence,
            motifs,
            owner.clone(),
            record.start as u64,
            &record.name,
            combine_strands,
            Self::START_SEARCH_CHUNK_SIZE,
        )?;
        Ok(ScannedReference { record, owner, motif_hits })
    }

    fn preflight_combined_partner_conflicts(
        work_queue: &VecDeque<ReferenceSearchSpace>,
        motifs: &[RegexMotif],
        combine_strands: bool,
    ) -> anyhow::Result<()> {
        if !(combine_strands && motifs.len() > 1) {
            return Ok(());
        }

        // Conflicting partner mappings must fail before the output writer is
        // created. Multi-motif combined mode is uncommon, so preflight each
        // owner with bounded memory and discard its hits. Ordinary and
        // single-motif modes retain the one-scan-per-active-owner path.
        for search_space in work_queue {
            Self::scan_motif_hits_with_chunk_size(
                &search_space.sequence,
                motifs,
                search_space.owner.clone(),
                search_space.record.start as u64,
                &search_space.record.name,
                combine_strands,
                Self::START_SEARCH_CHUNK_SIZE,
            )?;
        }
        Ok(())
    }

    fn new_with_regions(
        reference_sequences_lookup: ReferenceSequencesLookup,
        regions_bed_fp: &PathBuf,
        motifs: Vec<RegexMotif>,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let motif_search_context = Self::motif_search_context(&motifs);
        let regions_iter =
            BufReader::new(File::open(regions_bed_fp).with_context(|| {
                format!("failed to load regions at {regions_bed_fp:?}")
            })?)
            .lines()
            // skip comments/headers
            .skip_while(|r| {
                r.as_ref().map(|l| l.starts_with('#')).unwrap_or(false)
            })
            // change the lines into Errors
            .map(|r| r.map_err(|e| anyhow!("failed to read line, {e}")))
            // Parse the lines
            .map(|r| r.and_then(|l| BedRegion::parse_str(&l)))
            // grab the subsequences, also collect up the errors for invalid BED
            // lines
            .map(|r| {
                r.and_then(|bed_region| {
                    let chrom = bed_region.chrom.as_str();
                    let sequence_length = reference_sequences_lookup
                        .sequence_length_by_name(chrom)?;
                    if bed_region.interval.end > sequence_length {
                        bail!(
                            "region {}:{}-{} exceeds reference length {}",
                            bed_region.chrom,
                            bed_region.interval.start,
                            bed_region.interval.end,
                            sequence_length
                        )
                    }
                    let search_start = bed_region
                        .interval
                        .start
                        .saturating_sub(motif_search_context);
                    let search_end = bed_region
                        .interval
                        .end
                        .saturating_add(motif_search_context)
                        .min(sequence_length);
                    let owner = bed_region.interval.start - search_start
                        ..bed_region.interval.end - search_start;
                    let seq = reference_sequences_lookup
                        .get_subsequence_by_name(chrom, search_start..search_end)?;
                    let tid = reference_sequences_lookup
                        .name_to_chrom_id(chrom)
                        .ok_or_else(|| anyhow!("missing reference ID for {chrom}"))?;
                    let reference_record = ReferenceRecord::new(
                        tid,
                        search_start as u32,
                        seq.len() as u32,
                        bed_region.chrom.clone(),
                    );
                    Ok((reference_record, bed_region.name, seq, owner))
                })
            });

        // accumulators for the above iterator, could have done this all in a
        // fold, but with 3 accumulators this is easier to look at and
        // ends up being the same thing
        let mut work_queue = VecDeque::new();
        let mut region_queue = VecDeque::new();
        let mut failures = HashMap::new();

        let mut add_failure = |cause: String| {
            *failures.entry(cause).or_insert(0) += 1;
        };

        for res in regions_iter {
            match res {
                Ok((reference_record, region_name, subseq, owner)) => {
                    work_queue.push_back(ReferenceSearchSpace {
                        record: reference_record,
                        sequence: subseq,
                        owner,
                    });
                    region_queue.push_back(region_name);
                }
                Err(e) => {
                    add_failure(e.to_string());
                }
            }
        }

        if !failures.is_empty() {
            debug!("failure reasons while parsing regions BED file");
            for (cause, count) in
                failures.iter().sorted_by(|(_, a), (_, b)| a.cmp(b))
            {
                debug!("\t {cause}: {count}")
            }
        }

        if work_queue.is_empty() {
            bail!("no valid regions parsed");
        }

        assert_eq!(region_queue.len(), work_queue.len());
        Self::preflight_combined_partner_conflicts(
            &work_queue,
            &motifs,
            combine_strands,
        )?;
        let (
            curr_contig,
            curr_owner,
            curr_motif_hits,
            curr_position,
            curr_region_name,
        ) = loop {
            let (search_space, region_name) =
                match (work_queue.pop_front(), region_queue.pop_front()) {
                    (Some(search_space), Some(region_name)) => {
                        anyhow::Ok((search_space, region_name))
                    }
                    _ => bail!(
                        "didn't find at least 1 sequence with valid start \
                         position"
                    ),
                }?;
            let scanned = Self::scan_reference(
                search_space,
                &motifs,
                combine_strands,
            )
            .with_context(|| {
                format!("failed to scan entropy region {region_name}")
            })?;
            if let Some(first_hit) = scanned.motif_hits.first() {
                let start_position = first_hit
                    .pos
                    .checked_sub(scanned.record.start as u64)
                    .and_then(|position| usize::try_from(position).ok())
                    .expect("motif hit must be local to its reference slice");
                info!(
                    "starting with region {region_name} at 0-based position \
                     {} on contig {}",
                    first_hit.pos,
                    &scanned.record.name
                );
                break (
                    scanned.record,
                    scanned.owner,
                    scanned.motif_hits,
                    start_position,
                    region_name,
                );
            } else {
                info!("region {region_name} has no valid positions, skipping");
                continue;
            }
        };
        debug!(
            "parsed {} regions, starting with {} on contig {}",
            region_queue.len() + 1usize,
            &curr_region_name,
            curr_contig.name
        );

        Ok(Self {
            motifs,
            work_queue,
            region_names: region_queue,
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_hit_cursor: 0,
            curr_contig,
            curr_owner,
            curr_motif_hits,
            curr_region_name: Some(curr_region_name),
            combine_strands,
            done: false,
        })
    }

    fn new(
        reference_sequence_lookup: ReferenceSequencesLookup,
        motifs: Vec<RegexMotif>,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let mut work_queue = reference_sequence_lookup
            .into_reference_sequences()
            .into_iter()
            .map(|(record, sequence)| {
                let owner = 0..sequence.len();
                ReferenceSearchSpace { record, sequence, owner }
            })
            .collect::<VecDeque<_>>();
        Self::preflight_combined_partner_conflicts(
            &work_queue,
            &motifs,
            combine_strands,
        )?;

        let (
            curr_contig,
            curr_owner,
            curr_motif_hits,
            curr_position,
        ) = loop {
            let search_space = work_queue.pop_front().ok_or_else(|| {
                anyhow!(
                    "didn't find at least 1 sequence with a valid start \
                     position"
                )
            })?;
            let scanned =
                Self::scan_reference(search_space, &motifs, combine_strands)?;
            if let Some(first_hit) = scanned.motif_hits.first() {
                let pos = first_hit
                    .pos
                    .checked_sub(scanned.record.start as u64)
                    .and_then(|position| usize::try_from(position).ok())
                    .expect("motif hit must be local to its reference slice");
                info!(
                    "starting with contig {} at 0-based position {pos}",
                    &scanned.record.name
                );
                break (
                    scanned.record,
                    scanned.owner,
                    scanned.motif_hits,
                    pos,
                );
            } else {
                info!(
                    "contig {} had no valid motif positions, skipping..",
                    scanned.record.name
                );
            }
        };

        Ok(Self {
            motifs,
            work_queue,
            region_names: VecDeque::new(),
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_hit_cursor: 0,
            curr_contig,
            curr_owner,
            curr_motif_hits,
            curr_region_name: None,
            combine_strands,
            done: false,
        })
    }

    #[inline]
    fn take_hits_if_enough(
        &self,
        motif_hits: &[&MotifHit],
    ) -> Option<Vec<BaseAndPosition>> {
        let positions = motif_hits
            .iter()
            .take(self.num_positions)
            .map(|mh| (mh.base, mh.pos))
            .sorted_by(|(_, a), (_, b)| a.cmp(b))
            .collect::<Vec<BaseAndPosition>>();
        if positions.len() == self.num_positions {
            Some(positions)
        } else {
            None
        }
    }

    #[inline]
    fn enough_hits_for_window(
        &self,
        pos_hits: &[&MotifHit],
        neg_hits: &[&MotifHit],
    ) -> Option<(GenomeWindow, u64)> {
        if self.combine_strands {
            let next_reference_position = pos_hits
                .first()?
                .pos
                .checked_add(1)
                .expect("reference position overflow");
            let neg_to_pos = pos_hits
                .iter()
                .take(self.num_positions)
                .filter_map(|motif_hit| {
                    assert_eq!(
                        motif_hit.strand,
                        Strand::Positive,
                        "logic error!"
                    );
                    motif_hit.neg_position.map(|np| {
                        ((motif_hit.base, np), (motif_hit.base, motif_hit.pos))
                    })
                })
                .collect::<FxHashMap<BaseAndPosition, BaseAndPosition>>();
            if neg_to_pos.len() < self.num_positions {
                None
            } else {
                let interval = half_open_interval(
                    neg_to_pos
                        .keys()
                        .chain(neg_to_pos.values())
                        .map(|(_, position)| *position),
                );
                Some((
                    GenomeWindow::new_combine_strands(
                        interval,
                        self.num_positions,
                        neg_to_pos,
                    ),
                    next_reference_position,
                ))
            }
        } else {
            if pos_hits.len() >= self.num_positions
                || neg_hits.len() >= self.num_positions
            {
                let pos_positions = self.take_hits_if_enough(pos_hits);
                let neg_positions = self.take_hits_if_enough(neg_hits);
                match (pos_positions, neg_positions) {
                    (Some(p), Some(n)) => {
                        assert_eq!(p.len(), self.num_positions);
                        assert!(!p.is_empty());
                        assert_eq!(n.len(), self.num_positions);
                        assert!(!n.is_empty());
                        let leftmost_positive_ref_pos = p
                            .iter()
                            .min_by(|(_, a), (_, b)| a.cmp(b))
                            .map(|(_, p)| *p)
                            .unwrap();
                        let leftmost_negative_ref_pos = n
                            .iter()
                            .min_by(|(_, a), (_, b)| a.cmp(b))
                            .map(|(_, p)| *p)
                            .unwrap();
                        if leftmost_positive_ref_pos < leftmost_negative_ref_pos
                        {
                            // debug!("(+) is lefter, using {p:?}");
                            Some((
                                GenomeWindow::new_stranded(
                                    Some(p),
                                    None,
                                    self.num_positions,
                                ),
                                leftmost_positive_ref_pos
                                    .checked_add(1)
                                    .expect("reference position overflow"),
                            ))
                        } else if leftmost_negative_ref_pos
                            < leftmost_positive_ref_pos
                        {
                            // debug!("(-) is lefter, using {n:?}");
                            Some((
                                GenomeWindow::new_stranded(
                                    None,
                                    Some(n),
                                    self.num_positions,
                                ),
                                leftmost_negative_ref_pos
                                    .checked_add(1)
                                    .expect("reference position overflow"),
                            ))
                        } else {
                            assert_eq!(
                                leftmost_positive_ref_pos,
                                leftmost_negative_ref_pos
                            );
                            // debug!("they are the same, using {p:?} and
                            // {n:?}");
                            Some((
                                GenomeWindow::new_stranded(
                                    Some(p),
                                    Some(n),
                                    self.num_positions,
                                ),
                                leftmost_positive_ref_pos
                                    .checked_add(1)
                                    .expect("reference position overflow"),
                            ))
                        }
                    }
                    (Some(p), None) => {
                        // debug!("(+) only, using {p:?}");
                        let next_reference_position = p[0]
                            .1
                            .checked_add(1)
                            .expect("reference position overflow");
                        Some((
                            GenomeWindow::new_stranded(
                                Some(p),
                                None,
                                self.num_positions,
                            ),
                            next_reference_position,
                        ))
                    }
                    (None, Some(n)) => {
                        // debug!("(-) only, using {n:?}");
                        let next_reference_position = n[0]
                            .1
                            .checked_add(1)
                            .expect("reference position overflow");
                        Some((
                            GenomeWindow::new_stranded(
                                None,
                                Some(n),
                                self.num_positions,
                            ),
                            next_reference_position,
                        ))
                    }
                    _ => None,
                }
            } else {
                None
            }
        }
    }

    fn advance_cursor(&mut self, next_reference_position: u64) {
        let next_local_position = next_reference_position
            .checked_sub(self.curr_contig.start as u64)
            .and_then(|position| usize::try_from(position).ok())
            .expect("next motif cursor must be local to its reference slice");
        assert!(
            next_local_position > self.curr_position,
            "motif cursor must advance: {} -> {}",
            self.curr_position,
            next_local_position
        );
        assert!(
            next_local_position <= self.curr_owner.end,
            "motif cursor must remain inside its owner"
        );
        self.curr_position = next_local_position;
        while self.curr_hit_cursor < self.curr_motif_hits.len()
            && self.curr_motif_hits[self.curr_hit_cursor].pos
                < next_reference_position
        {
            self.curr_hit_cursor += 1;
        }
    }

    fn next_window(&mut self) -> Option<GenomeWindow> {
        while !self.at_end_of_contig() {
            let current_reference_position = self
                .curr_contig
                .start
                .checked_add(self.curr_position as u32)
                .expect("reference position overflow")
                as u64;
            let owner_end = self
                .curr_contig
                .start
                .checked_add(self.curr_owner.end as u32)
                .expect("reference position overflow")
                as u64;
            let window_end = current_reference_position
                .saturating_add(self.window_size as u64)
                .min(owner_end);
            let relative_end = self.curr_motif_hits[self.curr_hit_cursor..]
                .partition_point(|hit| hit.pos < window_end);
            let window_hits = &self.curr_motif_hits
                [self.curr_hit_cursor..self.curr_hit_cursor + relative_end];
            let pos_hits = window_hits
                .iter()
                .filter(|hit| hit.strand == Strand::Positive)
                .collect::<Vec<_>>();
            let neg_hits = window_hits
                .iter()
                .filter(|hit| hit.strand == Strand::Negative)
                .collect::<Vec<_>>();

            if let Some((entropy_window, next_reference_position)) =
                self.enough_hits_for_window(&pos_hits, &neg_hits)
            {
                self.advance_cursor(next_reference_position);
                return Some(entropy_window);
            }

            let next_reference_position = self.curr_motif_hits
                [self.curr_hit_cursor..]
                .iter()
                .find(|hit| hit.pos > current_reference_position)
                .map(|hit| hit.pos)
                .unwrap_or(owner_end);
            self.advance_cursor(next_reference_position);
        }
        None
    }

    #[cfg(test)]
    fn find_start_position(
        seq: &[char],
        motifs: &[RegexMotif],
    ) -> Option<usize> {
        Self::scan_motif_hits_with_chunk_size(
            seq,
            motifs,
            0..seq.len(),
            0,
            "reference",
            false,
            Self::START_SEARCH_CHUNK_SIZE,
        )
        .expect("non-combined motif scan cannot have partner conflicts")
        .first()
        .and_then(|hit| usize::try_from(hit.pos).ok())
    }

    #[inline]
    fn at_end_of_contig(&self) -> bool {
        self.curr_position >= self.curr_owner.end
            || self.curr_hit_cursor >= self.curr_motif_hits.len()
    }

    fn update_current_contig(&mut self) {
        'search: loop {
            let Some(search_space) = self.work_queue.pop_front() else {
                assert!(self.region_names.is_empty());
                self.done = true;
                break 'search;
            };
            assert!(
                self.region_names.is_empty()
                    || self.region_names.len() == self.work_queue.len() + 1,
                "region names must remain aligned with search spaces"
            );
            let region_name = if self.region_names.is_empty() {
                None
            } else {
                self.region_names.pop_front()
            };
            let scanned = Self::scan_reference(
                search_space,
                &self.motifs,
                self.combine_strands,
            )
            .expect(
                "combined motif partner conflicts must be rejected during \
                 construction",
            );
            if let Some(first_hit) = scanned.motif_hits.first() {
                let start_position = first_hit
                    .pos
                    .checked_sub(scanned.record.start as u64)
                    .and_then(|position| usize::try_from(position).ok())
                    .expect("motif hit must be local to its reference slice");
                self.curr_contig = scanned.record;
                self.curr_position = start_position;
                self.curr_hit_cursor = 0;
                self.curr_owner = scanned.owner;
                self.curr_motif_hits = scanned.motif_hits;
                self.curr_region_name = region_name;
                break 'search;
            }

            if let Some(region_name) = region_name {
                debug!(
                    "skipping region {region_name}, no valid positions for \
                     motifs {:?}",
                    &self.motifs
                );
            } else {
                debug!(
                    "skipping {}, no valid positions for motifs {:?}",
                    &scanned.record.name, &self.motifs
                );
            }
        }
    }

    pub(super) fn total_length(&self) -> usize {
        self.work_queue
            .iter()
            .map(|search_space| {
                search_space.owner.end - search_space.owner.start
            })
            .sum::<usize>()
            + self.curr_owner.end
            - self.curr_owner.start
    }
}

impl Iterator for SlidingWindows {
    type Item = Vec<GenomeWindows>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::with_capacity(self.batch_size);
        let mut windows = Vec::new();
        loop {
            // stopping conditions
            if self.done || batch.len() >= self.batch_size {
                break;
            }

            // grab the next window
            if let Some(entropy_window) = self.next_window() {
                windows.push(entropy_window);
            }

            // update conditions
            if self.at_end_of_contig() {
                // need to rotate the windows since we're moving on to another
                // contig
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                let finished_region =
                    std::mem::replace(&mut self.curr_region_name, None);
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        finished_region,
                    );
                    batch.push(entropy_windows);
                }
                self.update_current_contig();
                continue;
            }

            // N.B. semantics, if the current region name is None, we're just
            // doing sliding windows over the genome and we can cut
            // this batch once the window size is the batch size.
            // otoh, if current_region is Some, we never cut the batch until
            // we've finished the contig so that an entire region ends up
            // in a single batch
            if self.curr_region_name.is_none()
                && windows.len() > self.batch_size
            {
                assert!(
                    self.region_names.is_empty(),
                    "region names should be empty here!"
                );
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        None,
                    );
                    batch.push(entropy_windows);
                }
            }
        }

        if !windows.is_empty() {
            assert!(
                self.region_names.is_empty(),
                "region names should be empty here also!"
            );
            let entropy_windows =
                GenomeWindows::new(self.curr_contig.tid, windows, None);
            batch.push(entropy_windows)
        }

        if batch.is_empty() {
            None
        } else {
            Some(batch)
        }
    }
}

#[derive(new, Debug)]
pub(super) struct MethylationEntropy {
    me_entropy: f32,
    num_reads: usize,
    interval: Range<u64>,
}

// todo make this an enum, one for regions
#[derive(new, Debug)]
pub(super) struct WindowEntropy {
    chrom_id: u32,
    pos_me_entropy: Option<MkResult<MethylationEntropy>>,
    neg_me_entropy: Option<MkResult<MethylationEntropy>>,
}

struct DescriptiveStats {
    mean_entropy: f32,
    median_entropy: f32,
    max_entropy: f32,
    min_entropy: f32,
    mean_num_reads: f32,
    max_num_reads: usize,
    min_num_reads: usize,
    failed_count: usize,
    successful_count: usize,
}

impl DescriptiveStats {
    fn mean(xs: &[f32]) -> f32 {
        xs.iter().sum::<f32>() / (xs.len() as f32)
    }

    fn new(
        measurements: &[f32],
        n_reads: &[usize],
        n_fails: usize,
        chrom_id: u32,
        interval: &Range<u64>,
    ) -> MkResult<Self> {
        if measurements.is_empty() {
            debug_assert!(
                n_reads.is_empty(),
                "measurements and reads should be empty together"
            );
            Err(MkError::EntropyZeroCoverage {
                chrom_id,
                start: interval.start,
                end: interval.end,
            })
        } else {
            debug_assert_eq!(
                measurements.len(),
                n_reads.len(),
                "measurements and n_reads should be the same length"
            );
            let mean_entropy = Self::mean(measurements);
            let median_entropy =
                percentile_linear_interp(measurements, 0.5f32)?;
            // safe because of above check
            let (min_entropy, max_entropy) = match measurements.iter().minmax()
            {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let mean_num_reads = Self::mean(
                &n_reads.iter().map(|&x| x as f32).collect::<Vec<_>>(),
            );
            let (min_num_reads, max_num_reads) = match n_reads.iter().minmax() {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let success_count = measurements.len();

            Ok(Self {
                mean_entropy,
                median_entropy,
                max_entropy,
                min_entropy,
                mean_num_reads,
                max_num_reads,
                min_num_reads,
                successful_count: success_count,
                failed_count: n_fails,
            })
        }
    }

    pub(super) fn to_row(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        strand: Strand,
        region_name: &str,
    ) -> String {
        use crate::util::TAB;

        format!(
            "\
            {chrom}{TAB}\
            {start}{TAB}\
            {end}{TAB}\
            {region_name}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}\n",
            self.mean_entropy,
            strand.to_char(),
            self.median_entropy,
            self.min_entropy,
            self.max_entropy,
            self.mean_num_reads,
            self.min_num_reads,
            self.max_num_reads,
            self.successful_count,
            self.failed_count
        )
    }
}

#[derive(new)]
pub(super) struct RegionEntropy {
    chrom_id: u32,
    interval: Range<u64>,
    pos_entropy_stats: MkResult<DescriptiveStats>,
    neg_entropy_stats: Option<MkResult<DescriptiveStats>>,
    region_name: String,
    window_entropies: Vec<WindowEntropy>,
}

#[derive(new)]
struct Message {
    mod_calls: FxHashMap<BaseAndPosition, BaseModCall>,
    reference_start: i64,
    reference_end: i64,
    strand: Strand,
    // _name: String,
}

fn process_bam_fp(
    bam_fp: &PathBuf,
    fetch_definition: FetchDefinition,
    caller: Arc<MultipleThresholdModCaller>,
    io_threads: usize,
) -> anyhow::Result<Vec<Message>> {
    let mut reader = bam::IndexedReader::from_path(bam_fp)?;
    reader.set_threads(io_threads)?;
    reader.fetch(fetch_definition)?;

    let record_iter = reader
        .records()
        .filter_map(|r| r.ok())
        .filter(|record| {
            !record.is_unmapped()
                && !(record_is_not_primary(&record) || record.seq_len() == 0)
        })
        .filter_map(|record| {
            String::from_utf8(record.qname().to_vec())
                .ok()
                .map(|name| (record, name))
        })
        .filter_map(|(record, name)| {
            match ModBaseInfo::new_from_record(&record) {
                Ok(modbase_info) => Some((modbase_info, record, name)),
                Err(run_error) => {
                    debug!(
                        "read {name}, failed to parse modbase info, \
                         {run_error}"
                    );
                    None
                }
            }
        });

    let mut messages = Vec::new();
    for (modbase_info, record, name) in record_iter {
        match ReadBaseModProfile::process_record(
            &record,
            &name,
            modbase_info,
            None,
            None,
            1,
        ) {
            Ok(profile) => {
                let position_calls = PositionModCalls::from_profile(&profile);
                let strands = position_calls
                    .iter()
                    .map(|p| p.mod_strand)
                    .collect::<HashSet<Strand>>();
                if strands.len() > 1 {
                    debug!("duplex not yet supported");
                } else {
                    let strand = if record.is_reverse() {
                        Strand::Negative
                    } else {
                        Strand::Positive
                    };
                    let mod_calls = position_calls
                        .into_iter()
                        .filter_map(|p| {
                            match (p.ref_position, p.alignment_strand) {
                                (Some(ref_pos), Some(aln_strand)) => {
                                    Some((p, ref_pos, aln_strand))
                                }
                                _ => None,
                            }
                        })
                        .map(|(p, ref_pos, _alignment_strand)| {
                            let mod_base_call = caller
                                .call(&p.canonical_base, &p.base_mod_probs);
                            ((p.canonical_base, ref_pos as u64), mod_base_call)
                        })
                        .collect::<FxHashMap<BaseAndPosition, BaseModCall>>();
                    let msg = Message::new(
                        mod_calls,
                        record.reference_start(),
                        record.reference_end(),
                        strand,
                    );
                    messages.push(msg);
                }
            }
            Err(e) => {
                debug!("read {name} failed to extract modbase info, {e}");
            }
        };
    }
    Ok(messages)
}

pub(super) fn process_entropy_window(
    mut entropy_windows: GenomeWindows,
    min_coverage: u32,
    max_filtered_positions: usize,
    io_threads: usize,
    caller: Arc<MultipleThresholdModCaller>,
    bam_fps: &[PathBuf],
) -> anyhow::Result<EntropyCalculation> {
    let bam_fp = &bam_fps[0];
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let chrom_id = entropy_windows.chrom_id;
    drop(reader);

    let results = bam_fps
        .into_par_iter()
        .map(|fp| {
            process_bam_fp(
                fp,
                entropy_windows.get_fetch_definition(),
                caller.clone(),
                io_threads,
            )
        })
        .collect::<Vec<anyhow::Result<Vec<Message>>>>();

    for message_result in results {
        match message_result {
            Ok(messages) => {
                for message in messages {
                    entropy_windows.entropy_windows.par_iter_mut().for_each(
                        |window| {
                            window.add_read_to_patterns(
                                &message.mod_calls,
                                message.reference_start,
                                message.reference_end,
                                message.strand,
                                max_filtered_positions,
                            )
                        },
                    );
                }
            }
            Err(e) => {
                debug!("failed to run bam {e}");
            }
        }
    }

    Ok(entropy_windows.into_entropy_calculation(chrom_id, min_coverage))
}

#[derive(new, Debug)]
struct BedRegion {
    chrom: String,
    interval: Range<usize>,
    name: String,
}

impl BedRegion {
    fn parser(raw: &str) -> IResult<&str, Self> {
        let n_parts = raw.split('\t').count();
        let (rest, chrom) = crate::parsing_utils::consume_string(raw)?;
        let (rest, start) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, stop) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, name) = if n_parts == 3 {
            (rest, format!("{chrom}:{start}-{stop}"))
        } else {
            let (rest, _leading_tab) = multispace1(rest)?;
            crate::parsing_utils::consume_string_spaces(rest)?
        };

        let interval = (start as usize)..(stop as usize);
        let this = Self { chrom, interval, name };
        Ok((rest, this))
    }

    fn parse_str(raw: &str) -> anyhow::Result<Self> {
        Self::parser(raw)
            .map_err(|e| anyhow!("failed to parse {raw} into BED3 line, {e}"))
            .and_then(|(_, this)| {
                if this.interval.end > this.interval.start {
                    Ok(this)
                } else {
                    bail!("end must be after start")
                }
            })
    }
}

#[cfg(test)]
mod entropy_mod_tests {
    use crate::entropy::methylation_entropy::EntropySymbol;
    use crate::entropy::{
        BedRegion, GenomeWindow, GenomeWindows, ReferenceSearchSpace,
        SlidingWindows,
    };
    use crate::mod_bam::BaseModCall;
    use crate::mod_base_code::{DnaBase, ModCodeRepr};
    use crate::motifs::motif_bed::RegexMotif;
    use crate::util::{ReferenceRecord, Strand};
    use rayon::prelude::*;
    use rayon::ThreadPoolBuilder;
    use rustc_hash::FxHashMap;
    use std::collections::{BTreeSet, HashSet, VecDeque};

    fn combined_window_with_code_count(code_count: usize) -> GenomeWindow {
        let read_patterns = (0..code_count)
            .map(|idx| {
                vec![BaseModCall::Modified(1.0, ModCodeRepr::ChEbi(idx as u32))]
            })
            .collect();
        GenomeWindow::CombineStrands {
            interval: 0..1,
            neg_to_pos_positions: FxHashMap::default(),
            read_patterns,
            position_valid_coverages: vec![code_count as u32],
        }
    }

    fn mixed_codes() -> Vec<ModCodeRepr> {
        vec![
            ModCodeRepr::Code('z'),
            ModCodeRepr::ChEbi(900),
            ModCodeRepr::Code('a'),
            ModCodeRepr::ChEbi(1),
            ModCodeRepr::Code('m'),
            ModCodeRepr::ChEbi(42),
            ModCodeRepr::Code('b'),
            ModCodeRepr::ChEbi(7),
            ModCodeRepr::Code('q'),
            ModCodeRepr::ChEbi(3),
            ModCodeRepr::Code('x'),
            ModCodeRepr::ChEbi(100),
            ModCodeRepr::Code('c'),
            ModCodeRepr::ChEbi(2),
            ModCodeRepr::Code('d'),
            ModCodeRepr::ChEbi(10),
            ModCodeRepr::Code('e'),
        ]
    }

    fn stranded_lookup_window(
        pos_codes: &[ModCodeRepr],
        neg_codes: &[ModCodeRepr],
    ) -> GenomeWindow {
        let make_patterns = |codes: &[ModCodeRepr]| {
            codes
                .iter()
                .map(|code| vec![BaseModCall::Modified(1.0, *code)])
                .collect::<Vec<_>>()
        };
        GenomeWindow::Stranded {
            pos_interval: None,
            neg_interval: None,
            pos_positions: None,
            neg_positions: None,
            pos_read_patterns: make_patterns(pos_codes),
            neg_read_patterns: make_patterns(neg_codes),
            pos_position_valid_coverages: Vec::new(),
            neg_position_valid_coverages: Vec::new(),
        }
    }

    fn p_q_w_calls(codes: &[ModCodeRepr]) -> Vec<Vec<BaseModCall>> {
        let p = codes
            .iter()
            .map(|code| BaseModCall::Modified(1.0, *code))
            .collect::<Vec<_>>();
        let mut q = p.clone();
        q[0] = BaseModCall::Canonical(1.0);
        let mut wildcard = p.clone();
        wildcard[0] = BaseModCall::Filtered;
        vec![p, q, wildcard]
    }

    fn oracle_window(
        codes: &[ModCodeRepr],
        pos_order: &[usize; 3],
        neg_order: &[usize; 3],
    ) -> GenomeWindow {
        let patterns = p_q_w_calls(codes);
        let reorder = |order: &[usize; 3]| {
            order.iter().map(|idx| patterns[*idx].clone()).collect::<Vec<_>>()
        };
        let mut coverages = vec![3; codes.len()];
        coverages[0] = 2;
        let interval_end = codes.len().saturating_sub(1) as u64;
        GenomeWindow::Stranded {
            pos_interval: Some(0..interval_end),
            neg_interval: Some(100..100 + interval_end),
            pos_positions: None,
            neg_positions: None,
            pos_read_patterns: reorder(pos_order),
            neg_read_patterns: reorder(neg_order),
            pos_position_valid_coverages: coverages.clone(),
            neg_position_valid_coverages: coverages,
        }
    }

    fn entropy_snapshot(window: &GenomeWindow) -> [(f32, usize); 2] {
        let entropy = window.into_entropy(7, 2);
        let pos = entropy.pos_me_entropy.unwrap().unwrap();
        let neg = entropy.neg_me_entropy.unwrap().unwrap();
        [(pos.me_entropy, pos.num_reads), (neg.me_entropy, neg.num_reads)]
    }

    fn motif_sequence(length: usize, start: usize, motif: &str) -> Vec<char> {
        let mut sequence = vec!['N'; length];
        for (idx, base) in motif.chars().enumerate() {
            sequence[start + idx] = base;
        }
        sequence
    }

    fn sliding_windows_for_test(
        sequence: &str,
        motifs: Vec<RegexMotif>,
        combine_strands: bool,
        num_positions: usize,
    ) -> SlidingWindows {
        let sequence = sequence.chars().collect::<Vec<_>>();
        let sequence_length = sequence.len();
        let scanned = SlidingWindows::scan_reference(
            ReferenceSearchSpace {
                record: ReferenceRecord::new(
                    0,
                    0,
                    sequence_length as u32,
                    "chr1".to_string(),
                ),
                sequence,
                owner: 0..sequence_length,
            },
            &motifs,
            combine_strands,
        )
        .unwrap();
        let curr_position = scanned.motif_hits.first().unwrap().pos as usize;
        SlidingWindows {
            motifs,
            work_queue: VecDeque::new(),
            region_names: VecDeque::new(),
            window_size: sequence_length,
            num_positions,
            batch_size: 1,
            curr_position,
            curr_hit_cursor: 0,
            curr_contig: scanned.record,
            curr_owner: scanned.owner,
            curr_motif_hits: scanned.motif_hits,
            curr_region_name: None,
            combine_strands,
            done: false,
        }
    }

    #[test]
    fn one_position_window_is_one_base_half_open() {
        let mut window =
            GenomeWindow::new_stranded(Some(vec![(DnaBase::C, 9)]), None, 1);
        window.inc_coverage(0, &Strand::Positive);
        window
            .add_pattern(&Strand::Positive, vec![BaseModCall::Canonical(1.0)]);

        let entropy = window.into_entropy(0, 1);
        assert_eq!(entropy.pos_me_entropy.unwrap().unwrap().interval, 9..10);
    }

    #[test]
    fn reverse_anchor_at_final_base_is_one_base_half_open() {
        let mut window =
            GenomeWindow::new_stranded(None, Some(vec![(DnaBase::C, 3)]), 1);
        window.inc_coverage(0, &Strand::Negative);
        window
            .add_pattern(&Strand::Negative, vec![BaseModCall::Canonical(1.0)]);

        let entropy = window.into_entropy(0, 1);
        assert_eq!(entropy.neg_me_entropy.unwrap().unwrap().interval, 3..4);
    }

    #[test]
    fn entropy_start_position_is_global_after_first_search_chunk() {
        let sequence = motif_sequence(20_010, 15_000, "GATC");
        let motif = RegexMotif::parse_string("GATC", 1).unwrap();

        assert_eq!(
            super::SlidingWindows::find_start_position(&sequence, &[motif]),
            Some(15_001)
        );
    }

    #[test]
    fn entropy_start_position_finds_motif_across_search_chunk_seam() {
        let sequence = motif_sequence(10_010, 9_998, "GATC");
        let motif = RegexMotif::parse_string("GATC", 1).unwrap();

        assert_eq!(
            super::SlidingWindows::find_start_position(&sequence, &[motif]),
            Some(9_999)
        );
    }

    #[test]
    fn fetch_range_uses_global_min_start_and_max_exclusive_end() {
        let first = GenomeWindow::new_combine_strands(
            10..21,
            0,
            FxHashMap::default(),
        );
        let second = GenomeWindow::new_combine_strands(
            15..17,
            0,
            FxHashMap::default(),
        );
        let windows = GenomeWindows::new(0, vec![first, second], None);

        assert_eq!(windows.get_range(), 10..21);
    }

    #[test]
    fn combined_window_advances_from_owned_anchor_not_left_partner() {
        let motifs = vec![RegexMotif::parse_string("GATC", 3).unwrap()];
        let mut windows =
            sliding_windows_for_test("GATC", motifs, true, 1);

        assert!(windows.next_window().is_some());
        assert_eq!(windows.curr_position, 4);
    }

    #[test]
    fn overlapping_motifs_do_not_duplicate_non_combined_anchor() {
        let motifs = vec![
            RegexMotif::parse_string("CG", 0).unwrap(),
            RegexMotif::parse_string("CGN", 0).unwrap(),
        ];
        let mut windows =
            sliding_windows_for_test("CGA", motifs, false, 2);

        assert!(windows.next_window().is_none());
    }

    #[test]
    fn combined_motifs_reject_two_anchors_for_one_negative_partner() {
        let motifs = vec![
            RegexMotif::parse_string("GATC", 1).unwrap(),
            RegexMotif::parse_string("GATC", 2).unwrap(),
        ];
        let hits = vec![
            super::MotifHit::new(
                3,
                Some(10),
                Strand::Positive,
                DnaBase::C,
                0,
            ),
            super::MotifHit::new(
                5,
                Some(10),
                Strand::Positive,
                DnaBase::C,
                1,
            ),
        ];

        let error = SlidingWindows::sort_and_dedup_motif_hits(
            hits,
            &motifs,
            true,
            "chr1",
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains("negative C partner 10"), "{error}");
        assert!(error.contains("positive anchors 3"), "{error}");
        assert!(error.contains("and 5"), "{error}");
        assert!(error.contains("motif GATC,1"), "{error}");
        assert!(error.contains("motif GATC,2"), "{error}");
    }

    #[test]
    fn combined_motifs_deduplicate_identical_anchor_partner_pairs() {
        let motifs = vec![
            RegexMotif::parse_string("CG", 0).unwrap(),
            RegexMotif::parse_string("CGN", 0).unwrap(),
        ];
        let hits = vec![
            super::MotifHit::new(
                4,
                Some(5),
                Strand::Positive,
                DnaBase::C,
                0,
            ),
            super::MotifHit::new(
                4,
                Some(5),
                Strand::Positive,
                DnaBase::C,
                1,
            ),
        ];

        let deduped = SlidingWindows::sort_and_dedup_motif_hits(
            hits,
            &motifs,
            true,
            "chr1",
        )
        .unwrap();
        assert_eq!(deduped.len(), 1);
    }

    #[test]
    fn nine_distinct_modification_codes_have_atomic_states() {
        let window = combined_window_with_code_count(9);
        let lookup = window.get_mod_code_lookup();
        assert_eq!(lookup.len(), 9);
        for id in 1..=9 {
            assert_eq!(
                lookup.get(&ModCodeRepr::ChEbi(id - 1)),
                Some(&EntropySymbol::called(id))
            );
        }
    }

    #[test]
    fn ten_distinct_modification_codes_have_atomic_states() {
        let window = combined_window_with_code_count(10);
        assert_eq!(window.get_mod_code_lookup().len(), 10);
    }

    #[test]
    fn mixed_code_strand_union_preserves_current_btree_ranking() {
        let codes = mixed_codes();
        let first = stranded_lookup_window(&codes[..8], &codes[8..]);
        let mut reversed_pos = codes[8..].to_vec();
        reversed_pos.reverse();
        let mut reversed_neg = codes[..8].to_vec();
        reversed_neg.reverse();
        let permuted = stranded_lookup_window(&reversed_pos, &reversed_neg);
        let first_lookup = first.get_mod_code_lookup();
        let permuted_lookup = permuted.get_mod_code_lookup();
        assert_eq!(first_lookup.len(), 17);
        assert_eq!(first_lookup, permuted_lookup);
        assert_eq!(
            first_lookup.keys().copied().collect::<HashSet<_>>(),
            codes.iter().copied().collect::<HashSet<_>>()
        );
        let distinct_symbols =
            first_lookup.values().copied().collect::<BTreeSet<_>>();
        let expected_symbols =
            (1..=17).map(EntropySymbol::called).collect::<BTreeSet<_>>();
        assert_eq!(distinct_symbols, expected_symbols);
        assert!(!distinct_symbols.contains(&EntropySymbol::CANONICAL));
        assert!(!distinct_symbols.contains(&EntropySymbol::FILTERED));

        let mut code_symbols = codes
            .iter()
            .filter_map(|code| match code {
                ModCodeRepr::Code(c) => Some((*c, first_lookup[code])),
                ModCodeRepr::ChEbi(_) => None,
            })
            .collect::<Vec<_>>();
        code_symbols.sort_by_key(|(code, _)| *code);
        assert!(code_symbols.windows(2).all(|pair| pair[0].1 < pair[1].1));
        let mut chebi_symbols = codes
            .iter()
            .filter_map(|code| match code {
                ModCodeRepr::ChEbi(id) => Some((*id, first_lookup[code])),
                ModCodeRepr::Code(_) => None,
            })
            .collect::<Vec<_>>();
        chebi_symbols.sort_by_key(|(code, _)| *code);
        assert!(chebi_symbols.windows(2).all(|pair| pair[0].1 < pair[1].1));

        let expected_entropy = [(1.0 / 17.0, 3), (1.0 / 17.0, 3)];
        assert_eq!(
            entropy_snapshot(&oracle_window(&codes, &[0, 1, 2], &[2, 0, 1])),
            expected_entropy
        );
        let mut reversed_codes = codes.clone();
        reversed_codes.reverse();
        assert_eq!(
            entropy_snapshot(&oracle_window(
                &reversed_codes,
                &[1, 2, 0],
                &[0, 2, 1],
            )),
            expected_entropy
        );
    }

    #[test]
    fn seventeen_codes_preserve_oracle_across_threads_and_permutations() {
        let codes = (1..=17).map(ModCodeRepr::ChEbi).collect::<Vec<_>>();
        let mut reversed_codes = codes.clone();
        reversed_codes.reverse();
        let expected = [(1.0 / 17.0, 3), (1.0 / 17.0, 3)];

        for threads in [1, 2, 4] {
            let windows = vec![
                oracle_window(&codes, &[0, 1, 2], &[2, 0, 1]),
                oracle_window(&reversed_codes, &[1, 2, 0], &[0, 2, 1]),
            ];
            let pool =
                ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let observed = pool.install(|| {
                windows.par_iter().map(entropy_snapshot).collect::<Vec<_>>()
            });
            assert_eq!(observed, vec![expected, expected]);
        }
    }

    #[test]
    fn test_bed_region_parsing() {
        let raw = "chr1\t100\t101\tfoo\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");
        let raw = "chr1\t100\t101\tfoo\t400\t.\tmorestuff\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");

        let raw = "chr20\t279148\t279507\tCpG: 39";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr20");
        assert_eq!(bed_region.interval, 279148usize..279507);
        assert_eq!(&bed_region.name, "CpG: 39");
    }
}
