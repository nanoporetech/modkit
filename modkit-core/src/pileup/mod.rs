use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use derive_new::new;
use indexmap::IndexSet;
use itertools::Itertools;
use log::{debug, error};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::{FetchDefinition, Read};
use rustc_hash::FxHashMap;
use sortedlist_rs::SortedList;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::errs::MkError;
use crate::interval_chunks::{FocusPositions, MultiChromCoordinates};
use crate::mod_bam::{BaseModCall, BaseModProbs, CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::motifs::motif_bed::MotifInfo;
use crate::read_cache::ReadCache;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::{percentile_linear_interp, Percenileable};
use crate::util::{
    get_query_name_string, get_stringable_aux, record_is_not_primary, SamTag,
    Strand, StrandRule,
};

pub(crate) mod duplex;
pub mod subcommand;

pub(crate) enum ThresholdingOptions {
    PerPosition { percentile: f32, min_coverage_per_position: usize },
    Global,
}

#[derive(Default)]
struct PositionModCalls {
    num_filtered: u32,
    num_canonical: u32,
    mod_calls: FxHashMap<ModCodeRepr, u32>,
}

impl PositionModCalls {
    fn add_canonical(&mut self) {
        self.num_canonical = self.num_canonical.saturating_add(1u32)
    }
    fn add_filtered(&mut self) {
        self.num_filtered = self.num_filtered.saturating_add(1u32);
    }
    fn incr_mod_code(&mut self, mod_code_repr: ModCodeRepr) {
        self.mod_calls
            .entry(mod_code_repr)
            .and_modify(|x| *x = x.saturating_add(1u32))
            .or_insert(1u32);
    }
}

#[derive(Debug, Clone)]
enum Feature {
    Delete,
    // Filtered,
    NoCall(DnaBase),
    ModCall(BaseModProbs, DnaBase),
}

impl Feature {
    fn from_base_mod_probs(
        base_mod_probs: BaseModProbs,
        read_base: DnaBase,
    ) -> Self {
        Self::ModCall(base_mod_probs, read_base)
    }
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct PositionStats {
    pub(crate) min_prob: f32,
    pub(crate) median_prob: f32,
    pub(crate) mean_prob: f32,
    pub(crate) std_prob: f32,
    pub(crate) max_prob: f32,
}

impl PositionStats {
    fn new_empty() -> Self {
        Self {
            min_prob: f32::NAN,
            median_prob: f32::NAN,
            max_prob: f32::NAN,
            mean_prob: f32::NAN,
            std_prob: f32::NAN,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub(crate) enum PositionThreshold {
    Single { thresh: f32, stats: PositionStats },
    CombineStrands { pos: f32, neg: f32 },
}

impl PositionThreshold {
    pub(crate) fn print_position_threshold(&self) -> String {
        match self {
            PositionThreshold::Single { thresh, .. } => thresh.to_string(),
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }

    pub(crate) fn min_prob(&self) -> String {
        match self {
            PositionThreshold::Single { stats, .. } => {
                stats.min_prob.to_string()
            }
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }

    pub(crate) fn median_prob(&self) -> String {
        match self {
            PositionThreshold::Single { stats, .. } => {
                stats.median_prob.to_string()
            }
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }

    pub(crate) fn mean_prob(&self) -> String {
        match self {
            PositionThreshold::Single { stats, .. } => {
                stats.mean_prob.to_string()
            }
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }

    pub(crate) fn std_prob(&self) -> String {
        match self {
            PositionThreshold::Single { stats, .. } => {
                stats.std_prob.to_string()
            }
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }

    pub(crate) fn max_prob(&self) -> String {
        match self {
            PositionThreshold::Single { stats, .. } => {
                stats.max_prob.to_string()
            }
            PositionThreshold::CombineStrands { .. } => todo!(),
        }
    }
}

#[derive(Debug, Copy, Clone, new)]
pub(crate) struct PileupFeatureCounts {
    pub raw_strand: char,
    pub filtered_coverage: u32,
    pub raw_mod_code: ModCodeRepr,
    pub fraction_modified: f32,
    pub n_canonical: u32,
    pub n_modified: u32,
    pub n_other_modified: u32,
    pub n_delete: u32,
    pub n_filtered: u32,
    pub n_diff: u32,
    pub n_nocall: u32,
    pub motif_idx: Option<usize>,
    pub position_threshold: Option<PositionThreshold>,
}

impl PileupFeatureCounts {
    fn new_empty(
        raw_strand: char,
        raw_mod_code: ModCodeRepr,
        motif_index: Option<usize>,
    ) -> Self {
        Self {
            raw_strand,
            filtered_coverage: 0,
            raw_mod_code,
            motif_idx: motif_index,
            fraction_modified: 0f32,
            n_canonical: 0,
            n_modified: 0,
            n_other_modified: 0,
            n_delete: 0,
            n_filtered: 0,
            n_diff: 0,
            n_nocall: 0,
            position_threshold: None,
        }
    }

    // could make this moniod
    fn combine_counts_ignore_strand(self, other: Self) -> Self {
        if self.raw_mod_code != other.raw_mod_code {
            error!(
                "shouldn't be combining counts with different mod codes!{} vs \
                 {}",
                self.raw_mod_code, other.raw_mod_code
            );
        }
        if (self.motif_idx.is_some() && other.motif_idx.is_some())
            && self.motif_idx != other.motif_idx
        {
            error!(
                "shouldn't be combining counts with different motif indices \
                 {:?} vs {:?}",
                self.motif_idx, other.motif_idx
            );
        }
        let n_modified = self.n_modified + other.n_modified;
        let n_canonical = self.n_canonical + other.n_canonical;
        let n_other_modified = self.n_other_modified + other.n_other_modified;
        let filtered_coverage =
            self.filtered_coverage + other.filtered_coverage;
        let n_delete = self.n_delete + other.n_delete;
        let n_filtered = self.n_filtered + other.n_filtered;
        let n_diff = self.n_diff + other.n_diff;
        let n_nocall = self.n_nocall + other.n_nocall;

        let fraction_modified = n_modified as f32 / filtered_coverage as f32;

        let motif_idx = self.motif_idx;
        Self::new(
            self.raw_strand,
            filtered_coverage,
            self.raw_mod_code,
            fraction_modified,
            n_canonical,
            n_modified,
            n_other_modified,
            n_delete,
            n_filtered,
            n_diff,
            n_nocall,
            motif_idx,
            None,
        )
    }

    fn strand(&self) -> Option<Strand> {
        match &self.raw_strand {
            '+' => Some(Strand::Positive),
            '-' => Some(Strand::Negative),
            _ => None,
        }
    }
}

impl From<BedMethylLine> for PileupFeatureCounts {
    fn from(item: BedMethylLine) -> PileupFeatureCounts {
        PileupFeatureCounts::new(
            item.strand.into(),
            item.valid_coverage.try_into().unwrap_or(0),
            item.raw_mod_code,
            item.frac_modified(),
            item.count_canonical.try_into().unwrap_or(0),
            item.count_methylated.try_into().unwrap_or(0),
            item.count_other.try_into().unwrap_or(0),
            item.count_delete.try_into().unwrap_or(0),
            item.count_fail.try_into().unwrap_or(0),
            item.count_diff.try_into().unwrap_or(0),
            item.count_nocall.try_into().unwrap_or(0),
            None,
            None,
        )
    }
}

#[allow(non_snake_case)]
#[derive(Default)]
struct Tally {
    n_delete: u32,
    // n_filtered: u32,
    basecall_counts: FxHashMap<DnaBase, u32>,
    // make this a Map<DnaBase, Vec<BaseModCall>>
    base_to_mod_probs: FxHashMap<DnaBase, SortedList<BaseModProbs>>,
}

impl Tally {
    fn add_feature(&mut self, feature: Feature) {
        match feature {
            // Feature::Filtered => unreachable!(),
            Feature::Delete => self.n_delete += 1,
            Feature::ModCall(probs, primary_base) => self
                .base_to_mod_probs
                .entry(primary_base)
                .or_insert(SortedList::new())
                .insert(probs),
            Feature::NoCall(dna_base) => {
                *self.basecall_counts.entry(dna_base).or_insert(0) += 1;
            }
        }
    }

    // all of the counts of calls (canonical and mod) that aren't
    // for the primary base of this mode code
    #[inline]
    fn diff_calls_count(
        &self,
        primary_base: &DnaBase,
        n_other_modcall: u32,
    ) -> u32 {
        let n_other_basecall = self
            .basecall_counts
            .iter()
            .filter_map(|(dna_base, count)| {
                if dna_base == primary_base {
                    None
                } else {
                    Some(*count)
                }
            })
            .sum::<u32>();
        // let n_other_modcall = self
        //     .modcall_counts
        //     .iter()
        //     .filter_map(|(dna_base, mod_counts)| {
        //         if dna_base == primary_base {
        //             None
        //         } else {
        //             // need counts here
        //             Some(mod_counts.values().map(|&x| x).sum::<u32>())
        //         }
        //     })
        //     .sum::<u32>();

        n_other_basecall + n_other_modcall
    }

    fn get_position_stats(
        &self,
        thresholding_options: &ThresholdingOptions,
    ) -> Option<FxHashMap<DnaBase, PositionStats>> {
        match thresholding_options {
            ThresholdingOptions::PerPosition { .. } => {
                let per_base_stats = self
                    .base_to_mod_probs
                    .iter()
                    .map(|(base, probs)| match probs.len() {
                        0 => (*base, PositionStats::new_empty()),
                        1..=2 => {
                            let min_prob = probs
                                .first()
                                .map(|x| x.argmax_base_mod_call().prob())
                                .unwrap();
                            let max_prob = probs
                                .last()
                                .map(|x| x.argmax_base_mod_call().prob())
                                .unwrap();
                            let median_prob = (min_prob + max_prob) / 2f32;
                            let mean_prob = (min_prob / max_prob) / 2f32;
                            let std_prob = {
                                let var = [min_prob, max_prob]
                                    .map(|x| (x - mean_prob).powi(2))
                                    .iter()
                                    .sum::<f32>()
                                    / 2f32;
                                var.sqrt()
                            };

                            let stats = PositionStats {
                                min_prob,
                                median_prob,
                                mean_prob,
                                std_prob,
                                max_prob,
                            };
                            (*base, stats)
                        }
                        _ => {
                            let min_prob = probs
                                .first()
                                .map(|x| x.argmax_base_mod_call().prob())
                                .unwrap();
                            let max_prob = probs
                                .last()
                                .map(|x| x.argmax_base_mod_call().prob())
                                .unwrap();
                            let median_prob =
                                percentile_linear_interp(&probs, 0.5f32)
                                    .unwrap();
                            let mean_prob = probs.mean();
                            let std_prob = probs.std(mean_prob);
                            let stats = PositionStats {
                                min_prob,
                                median_prob,
                                mean_prob,
                                std_prob,
                                max_prob,
                            };
                            (*base, stats)
                        }
                    })
                    .collect::<FxHashMap<DnaBase, PositionStats>>();
                Some(per_base_stats)
            }
            ThresholdingOptions::Global => None,
        }
    }

    fn get_position_caller(
        &self,
        thresholding_options: &ThresholdingOptions,
    ) -> Option<MultipleThresholdModCaller> {
        match thresholding_options {
            ThresholdingOptions::PerPosition {
                percentile,
                min_coverage_per_position,
            } => {
                let per_base_thresholds = self
                    .base_to_mod_probs
                    .iter()
                    .filter_map(|(dna_base, probs)| {
                        if probs.len() >= *min_coverage_per_position {
                            let thresh = match percentile_linear_interp(
                                &probs,
                                *percentile,
                            ) {
                                Ok(t) => t,
                                Err(
                                    MkError::PercentileNotEnoughDatapoints(_),
                                ) => 0.0f32,
                                _ => unreachable!(),
                            };

                            Some((*dna_base, thresh))
                        } else {
                            None
                        }
                    })
                    .collect::<HashMap<DnaBase, f32>>();

                Some(MultipleThresholdModCaller::new(
                    per_base_thresholds,
                    HashMap::new(),
                    0.0,
                ))
            }
            ThresholdingOptions::Global => None,
        }
    }

    fn get_mod_calls(
        &self,
        caller: &MultipleThresholdModCaller,
        position_caller: Option<&MultipleThresholdModCaller>,
        thresholding_options: &ThresholdingOptions,
    ) -> FxHashMap<DnaBase, PositionModCalls> {
        // todo rayon
        let iter: Box<dyn Iterator<Item = (DnaBase, BaseModCall)>> =
            match thresholding_options {
                ThresholdingOptions::PerPosition { .. } => Box::new(
                    self.base_to_mod_probs.iter().flat_map(|(base, probs)| {
                        let caller = match position_caller {
                            Some(x)
                                if x.per_base_thresholds.contains_key(base) =>
                            {
                                x
                            }
                            _ => caller,
                        };
                        probs.iter().map(|prob| {
                            let base_mod_call = caller.call(base, prob);
                            (*base, base_mod_call)
                        })
                    }),
                ),
                ThresholdingOptions::Global => Box::new(
                    self.base_to_mod_probs.iter().flat_map(|(base, probs)| {
                        probs.iter().map(|prob| {
                            let base_mod_call = caller.call(base, prob);
                            (*base, base_mod_call)
                        })
                    }),
                ),
            };
        iter.fold(
            FxHashMap::default(),
            |mut base_to_pos_calls, (base, call)| {
                let pos_calls = base_to_pos_calls
                    .entry(base)
                    .or_insert_with(|| PositionModCalls::default());
                match call {
                    BaseModCall::Canonical(_prob) => {
                        pos_calls.add_canonical();
                    }
                    BaseModCall::Modified(_prob, code) => {
                        pos_calls.incr_mod_code(code);
                    }
                    BaseModCall::Filtered => pos_calls.add_filtered(),
                }
                base_to_pos_calls
            },
        )
    }
}

struct FeatureVector<'a> {
    pos_tally: Tally,
    neg_tally: Tally,
    caller: &'a MultipleThresholdModCaller,
    thresholding_options: &'a ThresholdingOptions,
}

impl<'a> FeatureVector<'a> {
    pub(crate) fn new(
        caller: &'a MultipleThresholdModCaller,
        thresholding_options: &'a ThresholdingOptions,
    ) -> Self {
        Self {
            pos_tally: Tally::default(),
            neg_tally: Tally::default(),
            caller,
            thresholding_options,
        }
    }

    /// Add counts to the tally.
    pub(crate) fn add_feature(
        &mut self,
        alignment_strand: Strand,
        feature: Feature,
        read_strand: Strand,
        strand_rule: &StrandRule,
    ) {
        match strand_rule {
            StrandRule::Both => match (alignment_strand, read_strand) {
                (Strand::Positive, Strand::Positive) => {
                    self.pos_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Positive) => {
                    self.neg_tally.add_feature(feature)
                }

                (Strand::Positive, Strand::Negative) => {
                    self.neg_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Negative) => {
                    self.pos_tally.add_feature(feature)
                }
            },
            StrandRule::Positive => match (alignment_strand, read_strand) {
                (Strand::Positive, Strand::Positive) => {
                    self.pos_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Negative) => {
                    self.pos_tally.add_feature(feature)
                }
                _ => {}
            },
            StrandRule::Negative => match (alignment_strand, read_strand) {
                (Strand::Negative, Strand::Positive) => {
                    self.neg_tally.add_feature(feature)
                }

                (Strand::Positive, Strand::Negative) => {
                    self.neg_tally.add_feature(feature)
                }
                _ => {}
            },
        }
    }

    fn add_tally_to_counts(
        &self,
        counts: &mut Vec<PileupFeatureCounts>,
        tally: &Tally,
        strand: Strand,
        observed_mods: &FxHashMap<DnaBase, HashSet<ModCodeRepr>>,
        pileup_options: &PileupNumericOptions,
        motif_idxs: Option<&Vec<usize>>,
    ) {
        let position_caller =
            tally.get_position_caller(self.thresholding_options);
        let position_stats =
            tally.get_position_stats(self.thresholding_options);
        let base_to_mod_calls = tally.get_mod_calls(
            self.caller,
            position_caller.as_ref(),
            self.thresholding_options,
        );

        let iter = base_to_mod_calls.iter().map(
            |(primary_base, position_mod_calls)| {
                (
                    primary_base,
                    position_mod_calls,
                    tally.basecall_counts.get(primary_base).unwrap_or(&0),
                    position_caller
                        .as_ref()
                        .and_then(|caller| {
                            caller.per_base_thresholds.get(primary_base)
                        })
                        .and_then(|thresh| {
                            position_stats
                                .as_ref()
                                .and_then(|x| {
                                    x.get(primary_base)
                                        .map(|stats| stats)
                                        .copied()
                                })
                                .map(|stats| PositionThreshold::Single {
                                    thresh: *thresh,
                                    stats,
                                })
                        }),
                )
            },
        );
        for (primary_base, position_mod_calls, &n_nocall, position_threshold) in
            iter
        {
            let n_canonical = position_mod_calls.num_canonical;
            let n_filtered = position_mod_calls.num_filtered;
            let mod_code_counts = &position_mod_calls.mod_calls;

            // todo do this arithmatic at object construction
            let n_other_modcall = base_to_mod_calls
                .iter()
                .filter_map(|(dna_base, mod_calls)| {
                    if dna_base == primary_base {
                        None
                    } else {
                        Some(mod_calls.mod_calls.values().copied().sum::<u32>())
                    }
                })
                .sum::<u32>();

            let total_num_modified = mod_code_counts.values().sum::<u32>();
            let filtered_coverage = total_num_modified + n_canonical;
            if filtered_coverage == 0 {
                continue;
            }
            match pileup_options {
                PileupNumericOptions::Passthrough
                | PileupNumericOptions::Collapse(_) => {
                    for (&mod_code, &n_mod) in observed_mods
                        .get(primary_base)
                        .unwrap_or(&HashSet::new())
                        .iter()
                        .map(|mod_code| {
                            (
                                mod_code,
                                mod_code_counts.get(mod_code).unwrap_or(&0),
                            )
                        })
                    {
                        let n_diff = tally
                            .diff_calls_count(primary_base, n_other_modcall);
                        let n_other_mod =
                            total_num_modified.checked_sub(n_mod).unwrap_or(0);
                        let percent_modified =
                            n_mod as f32 / filtered_coverage as f32;

                        if let Some(idxs) = motif_idxs {
                            for &idx in idxs.iter() {
                                counts.push(PileupFeatureCounts {
                                    raw_strand: strand.to_char(),
                                    filtered_coverage,
                                    raw_mod_code: mod_code,
                                    fraction_modified: percent_modified,
                                    n_canonical,
                                    n_modified: n_mod,
                                    n_other_modified: n_other_mod,
                                    n_delete: tally.n_delete,
                                    n_filtered,
                                    n_diff,
                                    n_nocall,
                                    motif_idx: Some(idx),
                                    position_threshold,
                                });
                            }
                        } else {
                            counts.push(PileupFeatureCounts {
                                raw_strand: strand.to_char(),
                                filtered_coverage,
                                raw_mod_code: mod_code,
                                fraction_modified: percent_modified,
                                n_canonical,
                                n_modified: n_mod,
                                n_other_modified: n_other_mod,
                                n_delete: tally.n_delete,
                                n_filtered,
                                n_diff,
                                n_nocall,
                                motif_idx: None,
                                position_threshold,
                            });
                        }
                    }
                }
                PileupNumericOptions::Combine => {
                    let percent_modified =
                        total_num_modified as f32 / filtered_coverage as f32;
                    let n_diff =
                        tally.diff_calls_count(&primary_base, n_other_modcall);
                    if let Some(idxs) = motif_idxs.as_ref() {
                        for &idx in idxs.iter() {
                            counts.push(PileupFeatureCounts {
                                raw_strand: strand.to_char(),
                                filtered_coverage,
                                raw_mod_code: ModCodeRepr::any_mod_code(
                                    &primary_base,
                                ),
                                fraction_modified: percent_modified,
                                n_canonical,
                                n_modified: total_num_modified,
                                n_other_modified: 0,
                                n_delete: tally.n_delete,
                                n_filtered,
                                n_diff,
                                n_nocall,
                                motif_idx: Some(idx),
                                position_threshold,
                            })
                        }
                    } else {
                        counts.push(PileupFeatureCounts {
                            raw_strand: strand.to_char(),
                            filtered_coverage,
                            raw_mod_code: ModCodeRepr::any_mod_code(
                                &primary_base,
                            ),
                            fraction_modified: percent_modified,
                            n_canonical,
                            n_modified: total_num_modified,
                            n_other_modified: 0,
                            n_delete: tally.n_delete,
                            n_filtered,
                            n_diff,
                            n_nocall,
                            motif_idx: None,
                            position_threshold,
                        })
                    }
                }
            }
        }
    }

    pub fn decode(
        self,
        pos_observed_mods: &FxHashMap<DnaBase, HashSet<ModCodeRepr>>,
        neg_observed_mods: &FxHashMap<DnaBase, HashSet<ModCodeRepr>>,
        pileup_options: &PileupNumericOptions,
        positive_motif_idxs: Option<&Vec<usize>>,
        negative_motif_idxs: Option<&Vec<usize>>,
    ) -> Vec<PileupFeatureCounts> {
        let mut counts = Vec::new();
        // dbg!(&self.pos_tally, &pos_observed_mods);

        self.add_tally_to_counts(
            &mut counts,
            &self.pos_tally,
            Strand::Positive,
            pos_observed_mods,
            pileup_options,
            positive_motif_idxs,
        );
        self.add_tally_to_counts(
            &mut counts,
            &self.neg_tally,
            Strand::Negative,
            neg_observed_mods,
            pileup_options,
            negative_motif_idxs,
        );

        counts.sort_by(|a, b| match a.raw_strand.cmp(&b.raw_strand) {
            Ordering::Equal => a.raw_mod_code.cmp(&b.raw_mod_code),
            o @ _ => o,
        });
        counts
    }
}

#[inline]
fn select_pileup_feature_counts(
    mappings: Option<&FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>>,
    partition_key: PartitionKey,
    strand: Strand,
    motif_idx: usize,
) -> Vec<PileupFeatureCounts> {
    mappings
        .and_then(|hm| hm.get(&partition_key))
        .unwrap_or(&Vec::new())
        .iter()
        .filter(|pileup_feature_counts| {
            let strand_match = pileup_feature_counts.strand() == Some(strand);
            let motif_match =
                pileup_feature_counts.motif_idx == Some(motif_idx);
            strand_match && motif_match
        })
        .copied()
        .collect()
}

fn combine_strand_features(
    positive_motif_idxs_lut: &BTreeMap<u32, Vec<(MotifInfo, usize)>>,
    position_feature_counts: BTreeMap<
        u32,
        FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>,
    >,
) -> BTreeMap<u32, FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>> {
    let mut result = BTreeMap::new();
    for (positive_strand_pos, motifs_at_position) in positive_motif_idxs_lut {
        let positive_feature_mappings =
            position_feature_counts.get(&positive_strand_pos);

        // start summing up the motif counts
        'motif: for (motif, idx) in motifs_at_position {
            // this is the position on the negative strand corresponding to this
            // positive strand motif,
            // e.g. for CCGG, 0
            //   v
            // + CCGG
            // - GGCC
            // _____^ <- this position
            let negative_strand_pos =
                motif.negative_strand_position(*positive_strand_pos);
            if negative_strand_pos.is_none() {
                continue 'motif;
            }
            let negative_strand_pos = negative_strand_pos.unwrap();
            let negative_feature_mappings =
                position_feature_counts.get(&negative_strand_pos);
            let partition_keys = positive_feature_mappings
                .unwrap_or(&FxHashMap::default())
                .keys()
                .chain(
                    negative_feature_mappings
                        .unwrap_or(&FxHashMap::default())
                        .keys(),
                )
                .copied()
                .collect::<HashSet<PartitionKey>>();
            for partition_key in partition_keys {
                // gather the positive and negative strands PileupFeatureCounts
                // that will be combined together
                let positive_strand_features = select_pileup_feature_counts(
                    positive_feature_mappings,
                    partition_key,
                    Strand::Positive,
                    *idx,
                );
                let negative_strand_features = select_pileup_feature_counts(
                    negative_feature_mappings,
                    partition_key,
                    Strand::Negative,
                    *idx,
                );
                // group them by mod code, use BTreeMap here so that the mod
                // codes are in a consistent order
                let grouped_by_mod_code = positive_strand_features
                    .into_iter()
                    .chain(negative_strand_features)
                    .fold(BTreeMap::new(), |mut agg, feats| {
                        agg.entry(feats.raw_mod_code)
                            .or_insert(Vec::new())
                            .push(feats);
                        agg
                    });
                // could technically make this one giant chained call but..
                let mut combined = grouped_by_mod_code
                    .into_iter()
                    .map(|(mod_code, feature_counts)| {
                        assert_eq!(feature_counts.len(), 2);
                        feature_counts.into_iter().fold(
                            // use unknown/ambiguous strand because we're
                            // combining
                            PileupFeatureCounts::new_empty(
                                '.',
                                mod_code,
                                Some(*idx),
                            ),
                            |acc, next| {
                                acc.combine_counts_ignore_strand(next) // use moniod
                            },
                        )
                    })
                    .collect::<Vec<PileupFeatureCounts>>();
                result
                    .entry(*positive_strand_pos)
                    .or_insert(FxHashMap::default())
                    .entry(partition_key)
                    .or_insert(Vec::new())
                    .append(&mut combined)
            }
        }
    }

    result
}

#[derive(new)]
struct StrandPileup {
    pub(crate) bam_pileup: bam::pileup::Pileup,
    strand_rule: StrandRule,
}

#[derive(new)]
struct PileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    start_pos: u32,
    end_pos: u32,
    focus_positions: &'a FocusPositions,
}

impl<'a> Iterator for PileupIter<'a> {
    type Item = StrandPileup;

    fn next(&mut self) -> Option<Self::Item> {
        let mut pileup: Option<Self::Item> = None;
        while let Some(Ok(plp)) = self.pileups.next() {
            let off_end = plp.pos() >= self.end_pos;
            if off_end {
                // we're done
                return None;
            } else if plp.pos() < self.start_pos {
                // advance into region we're looking at
                continue;
            } else {
                let pos = plp.pos();
                if let Some(strand_rule) =
                    self.focus_positions.check_position(&pos)
                {
                    pileup = Some(StrandPileup::new(plp, strand_rule));
                    break;
                } else {
                    continue;
                }
            }
        }
        pileup
    }
}

#[derive(Debug, Hash, Eq, PartialEq, Copy, Clone, Ord, PartialOrd)]
pub enum PartitionKey {
    NoKey,
    Key(usize),
}

fn get_forward_read_base(
    alignment: &bam::pileup::Alignment,
    record: &bam::Record,
) -> Option<DnaBase> {
    alignment.qpos().and_then(|pos| {
        if pos >= record.seq_len() {
            debug!("Record position is not included in sequence?");
            None
        } else {
            DnaBase::parse(record.seq()[pos] as char).ok()
        }
    })
}

fn parse_tags_from_record(
    record: &bam::Record,
    tags: &[SamTag],
) -> Option<String> {
    let values = tags
        .iter()
        .map(|tag| get_stringable_aux(&record, tag))
        .collect::<Vec<Option<String>>>();
    let got_match = values.iter().any(|b| b.is_some());
    if !got_match {
        return None;
    }
    let key = values
        .into_iter()
        .map(|v| v.unwrap_or("missing".to_string()))
        .join("_");
    Some(key)
}

pub struct ModBasePileup {
    pub chrom_name: String,
    position_feature_counts:
        BTreeMap<u32, FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>>,
    pub(crate) skipped_records: usize,
    pub(crate) processed_records: usize,
    pub(crate) partition_keys: IndexSet<String>,
}

impl ModBasePileup {
    pub fn num_results(&self) -> usize {
        self.position_feature_counts.len()
    }

    // this might be slowing us down, use a BTreeMap instead
    #[inline]
    pub(crate) fn iter_counts_sorted(
        &self,
    ) -> impl Iterator<
        Item = (&u32, &FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>),
    > {
        // self.position_feature_counts.iter().sorted_by(|(x, _), (y, _)|
        // x.cmp(y))
        self.position_feature_counts.iter()
    }
}

pub enum PileupNumericOptions {
    Passthrough,
    Combine,
    Collapse(CollapseMethod),
}

impl PileupNumericOptions {
    fn get_collapse_method(&self) -> Option<&CollapseMethod> {
        match self {
            Self::Collapse(method) => Some(method),
            _ => None,
        }
    }
}

// todo make this function generic so it can be used for duplex
//  as well.
pub(crate) fn process_region_batch<T: AsRef<Path> + Copy + Sync>(
    chromosome_coordintes: &MultiChromCoordinates,
    bam_fp: T,
    caller: &MultipleThresholdModCaller,
    thresholding_options: &ThresholdingOptions,
    pileup_numeric_options: &PileupNumericOptions,
    force_allow: bool,
    combine_strands: bool,
    max_depth: u32,
    edge_filter: Option<&EdgeFilter>,
    partition_tags: Option<&Vec<SamTag>>,
) -> Vec<Result<ModBasePileup, String>> {
    // todo make this anyhow::Result
    chromosome_coordintes
        .0
        .par_iter()
        .map(|chrom_coords| {
            process_region(
                bam_fp,
                chrom_coords.chrom_tid,
                chrom_coords.start_pos,
                chrom_coords.end_pos,
                caller,
                thresholding_options,
                pileup_numeric_options,
                force_allow,
                combine_strands,
                max_depth,
                &chrom_coords.focus_positions,
                edge_filter,
                partition_tags,
            )
        })
        .collect()
}

fn process_region<T: AsRef<Path>>(
    bam_fp: T,
    chrom_tid: u32,
    start_pos: u32,
    end_pos: u32,
    caller: &MultipleThresholdModCaller,
    thresholding_options: &ThresholdingOptions,
    pileup_numeric_options: &PileupNumericOptions,
    force_allow: bool,
    combine_strands: bool,
    max_depth: u32,
    focus_positions: &FocusPositions,
    edge_filter: Option<&EdgeFilter>,
    partition_tags: Option<&Vec<SamTag>>,
    // multi_progress: MultiProgress,
) -> Result<ModBasePileup, String> {
    let mut bam_reader =
        bam::IndexedReader::from_path(bam_fp).map_err(|e| e.to_string())?;
    let chrom_name =
        String::from_utf8_lossy(bam_reader.header().tid2name(chrom_tid))
            .to_string();
    bam_reader
        .fetch(FetchDefinition::Region(
            chrom_tid as i32,
            start_pos as i64,
            end_pos as i64,
        ))
        .map_err(|e| e.to_string())?;

    let mut read_cache = ReadCache::new(
        pileup_numeric_options.get_collapse_method(),
        caller,
        edge_filter,
        force_allow,
    );
    let mut position_feature_counts = BTreeMap::<
        u32,
        FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>,
    >::new();
    // collection of all partition keys encountered, ordered so
    // we can use their index
    let mut partition_keys = IndexSet::new();
    let hts_pileup = {
        let mut tmp_pileup = bam_reader.pileup();
        tmp_pileup.set_max_depth(max_depth);
        tmp_pileup
    };
    let pileup_iter =
        PileupIter::new(hts_pileup, start_pos, end_pos, focus_positions);
    let mut dupe_reads = HashMap::new(); // optimize
    for pileup in pileup_iter {
        let pos = pileup.bam_pileup.pos();

        // make a mapping of partition keys to feature vectors for this position
        let mut feature_vectors = HashMap::new();

        // Also make mappings of the observed mod codes per partition key
        let mut pos_strand_observed_mod_codes = FxHashMap::<
            PartitionKey,
            FxHashMap<DnaBase, HashSet<ModCodeRepr>>,
        >::default();
        let mut neg_strand_observed_mod_codes = FxHashMap::<
            PartitionKey,
            FxHashMap<DnaBase, HashSet<ModCodeRepr>>,
        >::default();

        // used for warning about dupes, could make this a bloom filter for
        // better perf?
        let mut observed_read_ids_to_pos = HashMap::new(); // optimize

        let alignment_iter =
            pileup.bam_pileup.alignments().filter(|alignment| {
                if alignment.is_refskip() {
                    false
                } else {
                    let record = alignment.record();
                    !(record_is_not_primary(&record) || record.seq_len() == 0)
                }
            });
        for alignment in alignment_iter {
            assert!(!alignment.is_refskip());
            let record = alignment.record();
            let partition_key = if let Some(tags) = partition_tags {
                match parse_tags_from_record(&record, tags) {
                    Some(s) => {
                        if let Some(idx) = partition_keys.get_index_of(&s) {
                            PartitionKey::Key(idx)
                        } else {
                            let inserted = partition_keys.insert(s);
                            debug_assert!(inserted);
                            debug_assert!(partition_keys.len() > 0);
                            PartitionKey::Key(
                                partition_keys
                                    .len()
                                    .checked_sub(1)
                                    .unwrap_or(0),
                            )
                        }
                    }
                    None => PartitionKey::NoKey,
                }
            } else {
                PartitionKey::NoKey
            };

            // data structures we update per alignment/read
            let mut pos_strand_mod_codes_for_key =
                pos_strand_observed_mod_codes
                    .entry(partition_key)
                    .or_insert(FxHashMap::default());
            let mut neg_strand_mod_codes_for_key =
                neg_strand_observed_mod_codes
                    .entry(partition_key)
                    .or_insert(FxHashMap::default());
            let feature_vector = feature_vectors
                .entry(partition_key)
                .or_insert(FeatureVector::new(caller, thresholding_options));

            read_cache.add_mod_codes_for_record(
                &record,
                &mut pos_strand_mod_codes_for_key,
                &mut neg_strand_mod_codes_for_key,
            );

            // optimize, could use a smarter string implementation here
            if let Ok(read_name) = get_query_name_string(&record) {
                (*observed_read_ids_to_pos
                    .entry(read_name)
                    .or_insert(0usize)) += 1
            }

            // alignment stand is the strand the read is aligned to
            let alignment_strand = if record.is_reverse() {
                Strand::Negative
            } else {
                Strand::Positive
            };

            if alignment.is_del() {
                feature_vector.add_feature(
                    alignment_strand,
                    Feature::Delete,
                    Strand::Positive,
                    &pileup.strand_rule,
                );
                continue;
            }

            // not delete or skip, add base
            let read_base = get_forward_read_base(&alignment, &record);

            let read_base = if let Some(base) = read_base {
                if record.is_reverse() {
                    base.complement()
                } else {
                    base
                }
            } else {
                // skip because read base failed, should this read be added to
                // the skip list?
                continue;
            };

            match read_cache.get_mod_call(&record, pos, read_base) {
                // a read can report on the read-positive or read-negative
                // strand (see the docs for .get_mod_call above) so the
                // pos_call and neg_call below are _read oriented_, the
                // `read_strand` in add_feature (see the docs there too)
                // is meant to pass along the information regarding which
                // strand of a read the feature belongs to. In almost all
                // cases this is Positive, because we sequence single
                // stranded DNA. However, for duplex and other double-
                // stranded techs, you could have a read with a mod call on
                // the negative strand. You must pass along the
                // `alignment_strand` so that everything can be oriented to
                // the positive strand of the reference.
                (Some(pos_call), Some(neg_call)) => {
                    let pos_feature =
                        Feature::from_base_mod_probs(pos_call, read_base);
                    let neg_feature = Feature::from_base_mod_probs(
                        neg_call,
                        read_base.complement(),
                    );
                    feature_vector.add_feature(
                        alignment_strand,
                        pos_feature,
                        Strand::Positive,
                        &pileup.strand_rule,
                    );
                    feature_vector.add_feature(
                        alignment_strand,
                        neg_feature,
                        Strand::Negative,
                        &pileup.strand_rule,
                    );
                }
                (Some(pos_call), None) => {
                    let pos_feature =
                        Feature::from_base_mod_probs(pos_call, read_base);
                    feature_vector.add_feature(
                        alignment_strand,
                        pos_feature,
                        Strand::Positive,
                        &pileup.strand_rule,
                    );
                }
                (None, Some(neg_call)) => {
                    let neg_feature = Feature::from_base_mod_probs(
                        neg_call,
                        read_base.complement(),
                    );

                    feature_vector.add_feature(
                        alignment_strand,
                        neg_feature,
                        Strand::Negative,
                        &pileup.strand_rule,
                    );
                }
                (None, None) => feature_vector.add_feature(
                    alignment_strand,
                    Feature::NoCall(read_base),
                    Strand::Positive,
                    &pileup.strand_rule,
                ),
            }
        } // alignment loop
        let pileup_feature_counts = feature_vectors
            .into_iter()
            .map(|(partition_key, fv)| {
                let pos_strand_observed_mod_codes_for_key =
                    pos_strand_observed_mod_codes.get(&partition_key);
                let neg_strand_observed_mod_codes_for_key =
                    neg_strand_observed_mod_codes.get(&partition_key);

                let positive_motif_idxs =
                    focus_positions.get_positive_strand_motif_ids(&pos);
                let negative_motif_idxs =
                    focus_positions.get_negative_strand_motif_ids(&pos);
                (
                    partition_key,
                    fv.decode(
                        pos_strand_observed_mod_codes_for_key
                            .unwrap_or(&FxHashMap::default()),
                        neg_strand_observed_mod_codes_for_key
                            .unwrap_or(&FxHashMap::default()),
                        &pileup_numeric_options,
                        positive_motif_idxs.as_ref(),
                        negative_motif_idxs.as_ref(),
                    ),
                )
            })
            .collect::<FxHashMap<PartitionKey, Vec<PileupFeatureCounts>>>();

        position_feature_counts.insert(pos, pileup_feature_counts);
        observed_read_ids_to_pos
            .into_iter()
            .filter(|(_, count)| *count > 1usize)
            .for_each(|(read_id, count)| {
                dupe_reads.entry(read_id).or_insert(Vec::new()).push(count);
            })
    } // position loop

    let position_feature_counts = if combine_strands {
        match focus_positions {
            FocusPositions::MotifCombineStrands { positive_motifs, .. } => {
                combine_strand_features(
                    positive_motifs,
                    position_feature_counts,
                )
            }
            _ => {
                // todo multiprogress suspend
                error!(
                    "asked to combine strand information without any motifs"
                );
                position_feature_counts
            }
        }
    } else {
        position_feature_counts
    };

    let (processed_records, skipped_records) =
        read_cache.get_records_used_and_skipped();

    let should_warn = !dupe_reads.is_empty();
    for (read_id, counts) in dupe_reads {
        let avg_times =
            counts.iter().map(|c| *c as f32).sum::<f32>() / counts.len() as f32;
        debug!(
            "read {read_id} was observed multiple times, avg {avg_times} at \
             {} positions on contig {chrom_name}, between {start_pos} and \
             {end_pos}",
            counts.len()
        );
    }
    if should_warn {
        debug!("consider marking duplicate alignments");
    }

    Ok(ModBasePileup {
        chrom_name,
        position_feature_counts,
        processed_records,
        skipped_records,
        partition_keys,
    })
}

#[cfg(test)]
mod mod_pileup_tests {
    use std::collections::HashSet;

    use rust_htslib::bam::{self, Read};
    use rustc_hash::FxHashMap;

    use crate::mod_bam::BaseModProbs;
    use crate::mod_base_code::{HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE};
    use crate::pileup::{
        parse_tags_from_record, DnaBase, Feature, FeatureVector,
        PileupNumericOptions, StrandRule, ThresholdingOptions,
    };
    use crate::threshold_mod_caller::MultipleThresholdModCaller;
    use crate::util::{SamTag, Strand};

    fn canonical_probs_fact() -> BaseModProbs {
        let mut base_mod_probs =
            BaseModProbs::new_init(HYDROXY_METHYL_CYTOSINE, 0.0f32);
        base_mod_probs.add_base_mod_prob(METHYL_CYTOSINE, 0.0).unwrap();
        base_mod_probs
    }

    fn hydroxy_probs_fact() -> BaseModProbs {
        let mut base_mod_probs =
            BaseModProbs::new_init(HYDROXY_METHYL_CYTOSINE, 0.9941406f32);
        base_mod_probs.add_base_mod_prob(METHYL_CYTOSINE, 0.0f32).unwrap();
        base_mod_probs
    }

    fn methyl_probs_fact() -> BaseModProbs {
        let mut base_mod_probs =
            BaseModProbs::new_init(METHYL_CYTOSINE, 0.9941406f32);
        base_mod_probs
            .add_base_mod_prob(HYDROXY_METHYL_CYTOSINE, 0.0f32)
            .unwrap();
        base_mod_probs
    }

    #[test]
    fn test_feature_vector_basic() {
        let hmc = HYDROXY_METHYL_CYTOSINE;
        let mc = METHYL_CYTOSINE;
        let pos_observed_mods = FxHashMap::from_iter(
            [(DnaBase::C, HashSet::from([mc, hmc]))].into_iter(),
        );
        let neg_observed_mods = FxHashMap::default();
        let caller = MultipleThresholdModCaller::new_passthrough();
        let threshold_options = ThresholdingOptions::Global;

        let mut fv = FeatureVector::new(&caller, &threshold_options);
        fv.add_feature(
            Strand::Positive,
            Feature::NoCall(DnaBase::A),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(canonical_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(methyl_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(methyl_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::NoCall(DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &neg_observed_mods,
            &PileupNumericOptions::Passthrough,
            None,
            None,
        );
        // dbg!(&counts);
        assert_eq!(counts.len(), 2); // h and m, negative strand should not be there
        for pileup_counts in counts {
            assert_eq!(pileup_counts.filtered_coverage, 3);
            assert_eq!(pileup_counts.n_nocall, 1);
            assert_eq!(pileup_counts.n_diff, 1);
            assert_eq!(pileup_counts.raw_strand, Strand::Positive.to_char());
        }
        let mut fv = FeatureVector::new(&caller, &threshold_options);
        let neg_observed_mods = FxHashMap::from_iter(
            [(DnaBase::C, HashSet::from([mc, hmc]))].into_iter(),
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(canonical_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::ModCall(methyl_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &neg_observed_mods,
            &PileupNumericOptions::Passthrough,
            None,
            None,
        );
        assert_eq!(counts.len(), 4);
        counts
            .iter()
            .filter(|c| c.raw_strand == Strand::Negative.to_char())
            .for_each(|c| assert_eq!(c.n_diff, 2));
    }

    #[test]
    fn test_feature_vector_with_strand_rules() {
        let caller = MultipleThresholdModCaller::new_passthrough();
        let threshold_options = ThresholdingOptions::Global;
        let mut fv = FeatureVector::new(&caller, &threshold_options);

        let mc = METHYL_CYTOSINE;
        let pos_observed_mods =
            FxHashMap::from_iter([(DnaBase::C, HashSet::from([mc]))]);
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(methyl_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Positive,
        );
        // this feature should be ignored because it's on the wrong
        // strand
        fv.add_feature(
            Strand::Negative,
            Feature::ModCall(methyl_probs_fact(), DnaBase::C),
            Strand::Positive,
            &StrandRule::Positive,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &FxHashMap::default(),
            &PileupNumericOptions::Passthrough,
            None,
            None,
        );
        assert_eq!(counts.len(), 1);
        let count = &counts[0];
        // change alignment strand to Positive and this will be 2
        assert_eq!(count.n_modified, 1);
    }

    #[test]
    fn test_parse_tags_from_record() {
        let mut reader = bam::Reader::from_path(
            "../tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        )
        .unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let tags = [SamTag::parse(['H', 'P']), SamTag::parse(['R', 'G'])];
        let key = parse_tags_from_record(&record, &tags);
        assert_eq!(key, Some("1_A".to_string()));
        let tags = [SamTag::parse(['R', 'G']), SamTag::parse(['H', 'P'])];
        let key = parse_tags_from_record(&record, &tags);
        assert_eq!(key, Some("A_1".to_string()));
    }
}
