use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use anyhow::bail;
use common_macros::hash_map;
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator};

use log::{debug, error, info};
use rayon::prelude::*;

use crate::mod_bam::{prob_to_qual, BaseModCall, CollapseMethod, EdgeFilter};
use crate::mod_base_code::{BaseState, DnaBase, ModCodeRepr};
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs;
use crate::record_processor::WithRecords;
use crate::sample_probs::{ModHist, QualHist};
use crate::threshold_mod_caller::MultipleThresholdModCaller;

use crate::thresholds::calc_thresholds_per_base;
use crate::util::{get_master_progress_bar, Region};

/// Count statistics from a modBAM.
#[derive(Debug, new, PartialEq)]
pub struct ModSummary<'a> {
    /// For each canonical base, how many reads had
    /// base modification calls for this base.
    pub reads_with_mod_calls: HashMap<DnaBase, u64>,
    /// For each canonical base, how many of each base modification
    /// code were observed and not filtered out.
    pub mod_call_counts: HashMap<DnaBase, HashMap<BaseState, u64>>,
    /// For each canonical base, how many of each base modification
    /// code were observed but filtered out.
    pub filtered_mod_call_counts: HashMap<DnaBase, HashMap<BaseState, u64>>,
    /// Total number of reads used in the summary. Usually a summary is
    /// computed on a sub-sample of the reads in a modBAM (or a
    /// sub-region).
    pub total_reads_used: usize,
    /// Mapping of canonical base to the estimated base modification confidence
    /// threshold for base modification calls at that base.
    pub per_base_thresholds: HashMap<DnaBase, f32>,
    /// If a region is provided, this is a reference to that region.
    pub region: Option<&'a Region>,
    /// Mapping of which modcodes were observed for each base
    pub per_base_mod_codes: HashMap<DnaBase, HashSet<ModCodeRepr>>,
}

impl<'a> ModSummary<'a> {
    pub(crate) fn mod_bases(&self) -> String {
        self.mod_call_counts
            .keys()
            .map(|d| d.char().to_string())
            .collect::<Vec<String>>()
            .join(",")
    }

    fn unpack_counts(
        mod_hists: &[ModHist],
        base_totals: [u64; 4],
    ) -> HashMap<DnaBase, HashMap<BaseState, u64>> {
        let mut mod_call_counts =
            mod_hists.iter().fold(HashMap::new(), |mut agg, mod_hist| {
                let e =
                    agg.entry(mod_hist.dna_base).or_insert_with(HashMap::new);
                e.insert(
                    BaseState::Modified(mod_hist.mod_code),
                    mod_hist.total,
                );
                agg
            });
        for primary_base in (0..4).map(|x: usize| DnaBase::try_from(x).unwrap())
        {
            if base_totals[primary_base as usize] > 0 {
                let e = mod_call_counts
                    .entry(primary_base)
                    .or_insert_with(HashMap::new);
                e.insert(
                    BaseState::Canonical(primary_base),
                    base_totals[primary_base as usize],
                );
            }
        }

        mod_call_counts
    }

    fn threshold_qual_hist(
        qual_hist: &QualHist,
        per_base_thresholds: [f32; 4],
        mod_thresholds: &HashMap<ModCodeRepr, f32>,
        pass: bool,
    ) -> QualHist {
        let mut this = QualHist::default();
        for b in 0..4usize {
            let t = per_base_thresholds[b];
            let q = prob_to_qual(t);
            let counts = qual_hist.hist[b];
            let mut total = 0u64;
            let range = if pass {
                (q as usize)..256usize
            } else {
                0usize..(q as usize)
            };
            for idx in range {
                assert!(idx < counts.len());
                let count = counts[idx];
                total = total.saturating_add(count);
            }
            this.base_totals[b] = total;
        }
        for mod_hist in qual_hist.mods_hists.iter() {
            let t = mod_thresholds
                .get(&mod_hist.mod_code)
                .copied()
                .unwrap_or(per_base_thresholds[mod_hist.dna_base as usize]);
            let q = prob_to_qual(t);
            let mut total = 0u64;
            let range = if pass {
                (q as usize)..256usize
            } else {
                0usize..(q as usize)
            };
            for idx in range {
                assert!(idx < mod_hist.hist.len());
                let count = mod_hist.hist[idx];
                total = total.saturating_add(count);
            }
            this.mods_hists.push(ModHist {
                total,
                mod_code: mod_hist.mod_code,
                dna_base: mod_hist.dna_base,
                hist: [0u64; 256],
            });
        }
        this
    }

    pub(crate) fn from_qual_hist(
        qual_hist: QualHist,
        filter_percentile: f32,
        max_thresholds_per_base: Option<[f32; 4]>,
        region: Option<&'a Region>,
        base_thresholds: Option<[f32; 4]>,
        mod_thresholds: Option<&HashMap<ModCodeRepr, f32>>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Self> {
        let reads_with_mod_calls = hash_map!{
            DnaBase::A => qual_hist.num_records_with_base_mods[DnaBase::A as usize] as u64,
            DnaBase::C => qual_hist.num_records_with_base_mods[DnaBase::C as usize] as u64,
            DnaBase::G => qual_hist.num_records_with_base_mods[DnaBase::G as usize] as u64,
            DnaBase::T => qual_hist.num_records_with_base_mods[DnaBase::T as usize] as u64,
        }.into_iter().filter(|(_, c)| *c > 0).collect::<HashMap<DnaBase, u64>>();
        let base_thresholds = if let Some(x) = base_thresholds {
            x
        } else {
            if !qual_hist.has_probability_observations() {
                bail!(
                    "cannot calculate automatic thresholds because no \
                     modification probabilities were sampled"
                )
            }
            let mut agg = [0f32; 4];
            let base_level_threshs = qual_hist.base_level_percentiles(
                &[filter_percentile],
                multi_progress,
                false,
            );
            for b in (0usize..4).map(|x| DnaBase::try_from(x).unwrap()) {
                if let Some(res) = base_level_threshs.get(&b) {
                    let t = res[0].threshold;
                    if let Some(max_ts) = &max_thresholds_per_base {
                        let max_t = max_ts[b as usize];
                        agg[b as usize] = if t > max_t { max_t } else { t }
                    } else {
                        agg[b as usize] = t;
                    }
                }
            }
            agg
        };

        let mut base_totals = qual_hist.base_totals;
        for mod_hist in qual_hist.mods_hists.iter() {
            // debug!("{} has total={}", mod_hist.mod_code, mod_hist.total);
            let c = base_totals[mod_hist.dna_base as usize]
                .saturating_add(mod_hist.total);
            base_totals[mod_hist.dna_base as usize] = c;
        }

        let pass_qual_hist = Self::threshold_qual_hist(
            &qual_hist,
            base_thresholds,
            mod_thresholds.unwrap_or(&HashMap::new()),
            true,
        );
        let fail_qual_hist = Self::threshold_qual_hist(
            &qual_hist,
            base_thresholds,
            mod_thresholds.unwrap_or(&HashMap::new()),
            false,
        );

        let mod_call_counts = Self::unpack_counts(
            &pass_qual_hist.mods_hists,
            pass_qual_hist.base_totals,
        );
        let filtered_mod_call_counts = Self::unpack_counts(
            &fail_qual_hist.mods_hists,
            fail_qual_hist.base_totals,
        );
        let total_reads_used = qual_hist.ok_records;
        let per_base_thresholds = hash_map! {
             DnaBase::A => base_thresholds[DnaBase::A as usize],
             DnaBase::C => base_thresholds[DnaBase::C as usize],
             DnaBase::G => base_thresholds[DnaBase::G as usize],
             DnaBase::T => base_thresholds[DnaBase::T as usize],
        }
        .into_iter()
        .filter(|(b, _)| base_totals[*b as usize] > 0)
        .collect::<HashMap<DnaBase, f32>>();

        let per_base_mod_codes = qual_hist.mods_hists.iter().fold(
            HashMap::new(),
            |mut agg, mod_hist| {
                let e =
                    agg.entry(mod_hist.dna_base).or_insert_with(HashSet::new);
                e.insert(mod_hist.mod_code);
                agg
            },
        );

        Ok(Self {
            reads_with_mod_calls,
            mod_call_counts,
            filtered_mod_call_counts,
            total_reads_used,
            per_base_thresholds,
            region,
            per_base_mod_codes,
        })
    }
}

/// Compute summary statistics from the reads in a modBAM. See `ModSummary`
/// for more details.
pub fn summarize_modbam<'a>(
    bam_fp: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&'a Region>,
    filter_percentile: f32,
    filter_thresholds: Option<MultipleThresholdModCaller>,
    per_mod_thresholds: Option<HashMap<ModCodeRepr, f32>>,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    suppress_progress: bool,
) -> anyhow::Result<ModSummary<'a>> {
    let read_ids_to_base_mod_calls =
        get_sampled_read_ids_to_base_mod_probs::<ReadIdsToBaseModProbs>(
            bam_fp,
            threads,
            interval_size,
            sample_frac,
            num_reads,
            seed,
            region,
            collapse_method,
            edge_filter,
            position_filter,
            only_mapped,
            suppress_progress,
        )?;

    let threshold_caller = if let Some(ft) = filter_thresholds {
        // filter thresholds provided, use those
        ft
    } else {
        // calculate the filter thresholds at the requested percentile
        let pct = (filter_percentile * 100f32).floor();
        info!("calculating threshold at {pct}% percentile");
        calc_thresholds_per_base(
            &read_ids_to_base_mod_calls,
            filter_percentile,
            None,
            per_mod_thresholds,
        )?
    };

    sampled_reads_to_summary(
        read_ids_to_base_mod_calls,
        &threshold_caller,
        region,
        suppress_progress,
    )
}

pub fn sampled_reads_to_summary<'a>(
    read_ids_to_mod_calls: ReadIdsToBaseModProbs,
    threshold_caller: &MultipleThresholdModCaller,
    region: Option<&'a Region>,
    suppress_progress: bool,
) -> anyhow::Result<ModSummary<'a>> {
    let total_reads_used = read_ids_to_mod_calls.num_reads();
    let start_t = std::time::Instant::now();

    let pb = get_master_progress_bar(read_ids_to_mod_calls.num_reads());
    pb.set_message("compiling summary");
    if suppress_progress {
        pb.set_draw_target(indicatif::ProgressDrawTarget::hidden())
    }
    let read_summary_chunk = read_ids_to_mod_calls
        .inner
        .par_iter()
        .progress_with(pb)
        .map(|(_read_id, canonical_base_to_calls)| {
            let mut mod_call_counts = HashMap::new();
            let mut filtered_mod_call_counts = HashMap::new();
            let mut reads_with_mod_calls = HashMap::new();
            let mut observed_mods = HashMap::new();
            for (&canonical_base, base_modification_probs) in
                canonical_base_to_calls
            {
                *reads_with_mod_calls.entry(canonical_base).or_insert(0) += 1;
                let canonical_base_mod_counts = mod_call_counts
                    .entry(canonical_base)
                    .or_insert(HashMap::new());
                let canonical_base_filtered_mod_counts =
                    filtered_mod_call_counts
                        .entry(canonical_base)
                        .or_insert(HashMap::new());

                base_modification_probs
                    .iter()
                    .map(|bmp| {
                        // need the argmax base_mod_call here too so that we can
                        // add to the correct
                        // filtered category
                        // once the whole "ModCode" bits are refactored, this
                        // will no longer be a necessary
                        // match
                        let argmax_base_mod_call = bmp.argmax_base_mod_call();
                        let thresholded_call =
                            threshold_caller.call(&canonical_base, bmp);
                        observed_mods
                            .entry(canonical_base)
                            .or_insert(HashSet::new())
                            .extend(bmp.iter_probs().map(|(code, _)| *code));
                        (thresholded_call, argmax_base_mod_call)
                        // match (thresholded_call, base_mod_call) {
                        //     (Ok(bmc), Ok(arg_max_base_mod_call)) => {
                        //         // add the observed mod codes here so that we
                        // report on them even if they're
                        //         // never called
                        //         observed_mods.entry(canonical_base).
                        // or_insert(HashSet::new())
                        //             .extend(bmp.iter_probs().map(|(code, _)|
                        // *code));         Some((bmc,
                        // arg_max_base_mod_call))     }
                        //     (Err(e), Err(_)) => {
                        //         debug!(
                        //             "read {read_id} failed to make
                        // thresholded mod call {}",
                        //             e.to_string()
                        //         );
                        //         // expected failure
                        //         None
                        //     }
                        //     (Ok(_), Err(e)) | (Err(e), Ok(_)) => {
                        //         // logic error, until refactor the two errors
                        // here are the same
                        //         error!(
                        //             "both should error or neither should!
                        // {}",
                        // e.to_string()         );
                        //         None
                        //     }
                        // }
                    })
                    .for_each(|(threshold_call, argmax_call)| {
                        let agg = match (threshold_call, argmax_call) {
                            (BaseModCall::Canonical(_), _) => {
                                canonical_base_mod_counts
                                    .entry(BaseState::Canonical(canonical_base))
                                    .or_insert(0)
                            }
                            (BaseModCall::Modified(_, mod_code_repr), _) => {
                                canonical_base_mod_counts
                                    .entry(BaseState::Modified(mod_code_repr))
                                    .or_insert(0)
                            }
                            (
                                BaseModCall::Filtered,
                                BaseModCall::Canonical(_),
                            ) => canonical_base_filtered_mod_counts
                                .entry(BaseState::Canonical(canonical_base))
                                .or_insert(0),
                            (
                                BaseModCall::Filtered,
                                BaseModCall::Modified(_, mod_code_repr),
                            ) => canonical_base_filtered_mod_counts
                                .entry(BaseState::Modified(mod_code_repr))
                                .or_insert(0),
                            (BaseModCall::Filtered, BaseModCall::Filtered) => {
                                error!("should not get filtered argmax calls");
                                unreachable!(
                                    "should not get filtered argmax calls"
                                );
                            }
                        };
                        *agg += 1u64;
                    });
            }
            ReadSummaryChunk {
                reads_with_mod_calls,
                mod_call_counts,
                filtered_mod_call_counts,
                observed_mods,
            }
        })
        .reduce(|| ReadSummaryChunk::zero(), |a, b| a.op(b));

    let elap = start_t.elapsed();
    debug!("computing summary took {}s", elap.as_secs());

    let per_base_thresholds = threshold_caller
        .iter_thresholds()
        .map(|(b, t)| (*b, *t))
        .collect::<HashMap<DnaBase, f32>>();

    Ok(ModSummary::new(
        read_summary_chunk.reads_with_mod_calls,
        read_summary_chunk.mod_call_counts,
        read_summary_chunk.filtered_mod_call_counts,
        total_reads_used,
        per_base_thresholds,
        region,
        read_summary_chunk.observed_mods,
    ))
}

#[derive(Debug)]
struct ReadSummaryChunk {
    reads_with_mod_calls: HashMap<DnaBase, u64>,
    mod_call_counts: HashMap<DnaBase, HashMap<BaseState, u64>>,
    filtered_mod_call_counts: HashMap<DnaBase, HashMap<BaseState, u64>>,
    observed_mods: HashMap<DnaBase, HashSet<ModCodeRepr>>,
}

impl Moniod for ReadSummaryChunk {
    fn zero() -> Self {
        Self {
            reads_with_mod_calls: HashMap::new(),
            mod_call_counts: HashMap::new(),
            filtered_mod_call_counts: HashMap::new(),
            observed_mods: HashMap::new(),
        }
    }

    fn op(self, other: Self) -> Self {
        let mut mod_call_counts = self.mod_call_counts;
        let mut filtered_mod_call_counts = self.filtered_mod_call_counts;
        let mut total = self.reads_with_mod_calls;
        let mut observed_mods = self.observed_mods;

        total.op_mut(other.reads_with_mod_calls);
        mod_call_counts.op_mut(other.mod_call_counts);
        filtered_mod_call_counts.op_mut(other.filtered_mod_call_counts);
        observed_mods.op_mut(other.observed_mods);

        Self {
            reads_with_mod_calls: total,
            mod_call_counts,
            filtered_mod_call_counts,
            observed_mods,
        }
    }

    fn op_mut(&mut self, other: Self) {
        self.reads_with_mod_calls.op_mut(other.reads_with_mod_calls);
        self.mod_call_counts.op_mut(other.mod_call_counts);
        self.filtered_mod_call_counts.op_mut(other.filtered_mod_call_counts);
        self.observed_mods.op_mut(other.observed_mods);
    }

    fn len(&self) -> usize {
        todo!()
    }
}
