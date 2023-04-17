use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::anyhow;
use derive_new::new;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use crate::errs::RunError;
use crate::filter_thresholds::FilterThresholds;
use crate::mod_bam::BaseModCall;
use crate::mod_base_code::{DnaBase, ModCode};
use crate::monoid::Moniod;
use crate::reads_sampler::{
    get_sampled_read_ids_to_base_mod_calls, ReadIdsToBaseModCalls,
};
use crate::record_sampler::RecordSampler;
use crate::thresholds::{
    calc_thresholds_per_base, get_modbase_probs_from_bam, modbase_records,
};
use crate::util::{get_master_progress_bar, get_spinner, Region, Strand};

/// Count statistics from a modBAM.
#[derive(Debug, new, Eq, PartialEq)]
pub struct ModSummary<'a> {
    /// For each canonical base, how many reads had
    /// base modification calls for this base.
    pub reads_with_mod_calls: HashMap<DnaBase, u64>,
    /// For each canonical base, how many of each base modification
    /// code were observed and not filtered out.
    pub mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    /// For each canonical base, how many of each base modification
    /// code were observed but filtered out.
    pub filtered_mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    /// Total number of reads used in the summary. Usually a summary is computed
    /// on a sub-sample of the reads in a modBAM (or a sub-region).
    pub total_reads_used: usize,
    /// If a region is provided, this is a reference to that region.
    pub region: Option<&'a Region>,
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
    filter_thresholds: Option<FilterThresholds>,
) -> anyhow::Result<ModSummary<'a>> {
    let read_ids_to_base_mod_calls = get_sampled_read_ids_to_base_mod_calls(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
    )?;

    let filter_thresholds = if let Some(ft) = filter_thresholds {
        ft
    } else {
        // if sample_frac and num_reads are none, no sampling was performed
        // and we cannot estimate the filter threshold at the filter_percentile
        // so make a passthrough
        if sample_frac.is_none() && num_reads.is_none() {
            // no filtering
            FilterThresholds::new_passthrough()
        } else {
            // calculate the threshold at the given filter-percentile
            let pct = filter_percentile * 100f32;
            info!("calculating threshold at {pct}% percentile");
            calc_thresholds_per_base(
                &read_ids_to_base_mod_calls,
                filter_percentile,
                None,
            )?
        }
    };

    sampled_reads_to_summary(
        read_ids_to_base_mod_calls,
        &filter_thresholds,
        region,
    )
}

fn sampled_reads_to_summary<'a>(
    read_ids_to_mod_calls: ReadIdsToBaseModCalls,
    filter_thresholds: &FilterThresholds,
    region: Option<&'a Region>,
) -> anyhow::Result<ModSummary<'a>> {
    let total_reads_used = read_ids_to_mod_calls.num_reads();
    let start_t = std::time::Instant::now();

    let pb = get_master_progress_bar(read_ids_to_mod_calls.num_reads());
    pb.set_message("compiling summary");
    let read_summary_chunk = read_ids_to_mod_calls.inner
        .par_iter()
        .progress_with(pb)
        .map(|(_read_id, canonical_base_to_calls)| {
            let mut mod_call_counts = HashMap::new();
            let mut filtered_mod_call_counts = HashMap::new();
            let mut reads_with_mod_calls = HashMap::new();
            for (&canonical_base, base_mod_calls) in canonical_base_to_calls {
                *reads_with_mod_calls.entry(canonical_base).or_insert(0) += 1;
                let canonical_base_mod_counts = mod_call_counts
                    .entry(canonical_base)
                    .or_insert(HashMap::new());
                let canonical_base_filtered_mod_counts = filtered_mod_call_counts
                    .entry(canonical_base)
                    .or_insert(HashMap::new());

                let mod_code_for_canonical_base =
                    canonical_base.canonical_mod_code().unwrap();
                for base_mod_call in base_mod_calls {
                    let agg = match base_mod_call {
                        BaseModCall::Canonical(p) => {
                            if *p < filter_thresholds.get(&canonical_base) {
                                canonical_base_filtered_mod_counts
                                    .entry(mod_code_for_canonical_base)
                                    .or_insert(0)
                            } else {
                                canonical_base_mod_counts
                                    .entry(mod_code_for_canonical_base)
                                    .or_insert(0)
                            }
                        }
                        BaseModCall::Modified(p, mod_code) => {
                            if *p < filter_thresholds.get(&canonical_base) {
                                canonical_base_filtered_mod_counts
                                    .entry(*mod_code)
                                    .or_insert(0)
                            } else {
                                canonical_base_mod_counts
                                    .entry(*mod_code)
                                    .or_insert(0)
                            }
                        }
                        BaseModCall::Filtered => {
                            unreachable!(
                                "should not encounter filtered base mod calls here"
                            );
                        }
                    };
                    *agg += 1u64;
                }
            }
            ReadSummaryChunk {
                reads_with_mod_calls, mod_call_counts, filtered_mod_call_counts,
            }
        })
        .reduce(|| ReadSummaryChunk::zero(), |a, b| a.op(b));

    let elap = start_t.elapsed();
    debug!("computing summary took {}s", elap.as_secs());
    Ok(ModSummary::new(
        read_summary_chunk.reads_with_mod_calls,
        read_summary_chunk.mod_call_counts,
        read_summary_chunk.filtered_mod_call_counts,
        total_reads_used,
        region,
    ))
}

#[derive(Debug)]
struct ReadSummaryChunk {
    reads_with_mod_calls: HashMap<DnaBase, u64>,
    mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    filtered_mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
}

impl Moniod for ReadSummaryChunk {
    fn zero() -> Self {
        Self {
            reads_with_mod_calls: HashMap::new(),
            mod_call_counts: HashMap::new(),
            filtered_mod_call_counts: HashMap::new(),
        }
    }

    fn op(self, other: Self) -> Self {
        let mut mod_call_counts = self.mod_call_counts;
        let mut filtered_mod_call_counts = self.filtered_mod_call_counts;
        let mut total = self.reads_with_mod_calls;
        total.op_mut(other.reads_with_mod_calls);
        mod_call_counts.op_mut(other.mod_call_counts);
        filtered_mod_call_counts.op_mut(other.filtered_mod_call_counts);
        Self {
            reads_with_mod_calls: total,
            mod_call_counts,
            filtered_mod_call_counts,
        }
    }

    fn op_mut(&mut self, other: Self) {
        self.reads_with_mod_calls.op_mut(other.reads_with_mod_calls);
        self.mod_call_counts.op_mut(other.mod_call_counts);
        self.filtered_mod_call_counts
            .op_mut(other.filtered_mod_call_counts);
    }

    fn len(&self) -> usize {
        todo!()
    }
}
