use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::anyhow;
use derive_new::new;
use itertools::Itertools;
use log::{debug, error, info, warn};
use rust_htslib::bam;
use rust_htslib::bam::Read;

use crate::errs::RunError;
use crate::filter_thresholds::FilterThresholds;
use crate::mod_bam::BaseModCall;
use crate::mod_base_code::{DnaBase, ModCode};
use crate::reads_sampler::{
    get_sampled_read_ids_to_base_mod_calls, ReadIdsToBaseModCalls,
};
use crate::record_sampler::RecordSampler;
use crate::thresholds::{
    calc_thresholds_per_base, get_modbase_probs_from_bam, modbase_records,
};
use crate::util::{get_master_progress_bar, get_spinner, Region, Strand};

#[derive(Debug, new, Eq, PartialEq)]
pub struct ModSummary<'a> {
    pub reads_with_mod_calls: HashMap<DnaBase, usize>,
    pub mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub filtered_mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub total_reads_used: usize,
    pub region: Option<&'a Region>,
}

fn sampled_reads_to_summary<'a>(
    read_ids_to_mod_calls: ReadIdsToBaseModCalls,
    filter_thresholds: &FilterThresholds,
    region: Option<&'a Region>,
) -> anyhow::Result<ModSummary<'a>> {
    let total_reads_used = read_ids_to_mod_calls.num_reads();
    let mut reads_with_mod_calls = HashMap::new();
    let mut mod_call_counts = HashMap::new();
    let mut filtered_mod_call_counts = HashMap::new();
    let start_t = std::time::Instant::now();
    // todo this could be done in parallel
    for (_read_id, canonical_base_to_calls) in read_ids_to_mod_calls.inner {
        for (canonical_base, base_mod_calls) in canonical_base_to_calls {
            *reads_with_mod_calls.entry(canonical_base).or_insert(0) += 1;

            let canonical_base_mod_counts = mod_call_counts
                .entry(canonical_base)
                .or_insert(HashMap::new());
            let canonical_base_filtered_mod_counts = filtered_mod_call_counts
                .entry(canonical_base)
                .or_insert(HashMap::new());

            let mod_code_for_canonical_base =
                canonical_base.canonical_mod_code()?;
            for base_mod_call in base_mod_calls {
                let agg = match base_mod_call {
                    BaseModCall::Canonical(p) => {
                        if p < filter_thresholds.get(&canonical_base) {
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
                        if p < filter_thresholds.get(&canonical_base) {
                            canonical_base_filtered_mod_counts
                                .entry(mod_code)
                                .or_insert(0)
                        } else {
                            canonical_base_mod_counts
                                .entry(mod_code)
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
    }
    let elap = start_t.elapsed();
    info!("computing summary took {}s", elap.as_secs());
    Ok(ModSummary::new(
        reads_with_mod_calls,
        mod_call_counts,
        filtered_mod_call_counts,
        total_reads_used,
        region,
    ))
}

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
                threads,
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

// pub fn summarize_modbam_<T: AsRef<Path>>(
//     bam_fp: T,
//     threads: usize,
//     threshold: f32,
//     num_reads: Option<usize>,
// ) -> Result<ModSummary, RunError> {
//     let mut reader = bam::Reader::from_path(bam_fp)
//         .map_err(|e| RunError::new_input_error(e.to_string()))?;
//     reader.set_threads(threads).map_err(|e| {
//         RunError::new_failed(format!(
//             "failed to set threads on reader, {}",
//             e.to_string()
//         ))
//     })?;
//
//     let record_iter = modbase_records(reader.records());
//     let spinner = if let Some(n) = num_reads {
//         get_master_progress_bar(n)
//     } else {
//         get_spinner()
//     };
//
//     spinner.set_message("records processed");
//
//     let mut total_reads_used = 0;
//     let mut reads_with_mod_calls = HashMap::new();
//     let mut mod_call_counts = HashMap::new();
//     let mut filtered_mod_calls = HashMap::new();
//     for (i, modbase_info) in record_iter.enumerate() {
//         if modbase_info.is_empty() {
//             continue;
//         }
//
//         let (_converters, prob_iter) = modbase_info.into_iter_base_mod_probs();
//         for (canonical_base, strand, seq_pos_mod_probs) in prob_iter {
//             let canonical_base = match (DnaBase::parse(canonical_base), strand)
//             {
//                 (Err(_), _) => continue,
//                 (Ok(dna_base), Strand::Positive) => dna_base,
//                 (Ok(dna_base), Strand::Negative) => dna_base.complement(),
//             };
//             let count = reads_with_mod_calls.entry(canonical_base).or_insert(0);
//             *count += 1;
//             let mod_counts = mod_call_counts
//                 .entry(canonical_base)
//                 .or_insert(HashMap::new());
//             let filtered_counts = filtered_mod_calls
//                 .entry(canonical_base)
//                 .or_insert(HashMap::new());
//             for (_position, base_mod_probs) in
//                 seq_pos_mod_probs.pos_to_base_mod_probs
//             {
//                 let count = match base_mod_probs.base_mod_call() {
//                     BaseModCall::Canonical(p) => {
//                         if p > threshold {
//                             mod_counts
//                                 .entry(
//                                     canonical_base
//                                         .canonical_mod_code()
//                                         .unwrap(),
//                                 )
//                                 .or_insert(0)
//                         } else {
//                             filtered_counts
//                                 .entry(
//                                     canonical_base
//                                         .canonical_mod_code()
//                                         .unwrap(),
//                                 )
//                                 .or_insert(0)
//                         }
//                     }
//                     BaseModCall::Modified(p, mod_code) => {
//                         if p > threshold {
//                             mod_counts.entry(mod_code).or_insert(0)
//                         } else {
//                             filtered_counts.entry(mod_code).or_insert(0)
//                         }
//                     }
//                     BaseModCall::Filtered => {
//                         unreachable!("should not encounter filtered calls")
//                     }
//                 };
//                 *count += 1;
//             }
//         }
//         total_reads_used = i;
//         spinner.inc(1);
//         let done = num_reads.map(|n| i >= n).unwrap_or(false);
//         if done {
//             break;
//         }
//     }
//     spinner.finish_and_clear();
//
//     Ok(ModSummary {
//         reads_with_mod_calls,
//         mod_call_counts,
//         filtered_mod_calls,
//         total_reads_used: total_reads_used + 1,
//         region: None,
//     })
// }
