use std::sync::Arc;

use crate::dmr::bedmethyl::{aggregate_counts, BedMethylLine};
use crate::dmr::llr_model::{AggregatedCounts, ModificationCounts};
use crate::dmr::tabix::{ChromToSampleBMLines, MultiSampleIndex};
use crate::dmr::util::{DmrBatch, RegionOfInterest, RoiIter};
use crate::monoid::BorrowingMoniod;
use anyhow::{anyhow, bail};
use indicatif::ProgressBar;
use log::{debug, error};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

#[inline]
fn filter_sample_records<'a>(
    sample_records: &'a ChromToSampleBMLines,
    roi: &RegionOfInterest,
    sample_index: &MultiSampleIndex,
) -> FxHashMap<usize, Vec<&'a BedMethylLine>> {
    sample_records
        .get(&roi.dmr_interval.chrom)
        .map(|per_sample| {
            per_sample
                .iter()
                .map(|(sample, lines)| {
                    let overlaping_records = lines
                        .iter()
                        .filter(|record| {
                            roi.positions.contains(
                                &record.get_stranded_position(
                                    &sample_index.code_lookup,
                                ),
                            )
                        })
                        .collect::<Vec<&BedMethylLine>>();
                    (*sample, overlaping_records)
                })
                .collect::<FxHashMap<usize, Vec<&BedMethylLine>>>()
        })
        .unwrap_or_else(|| FxHashMap::default())
}

#[inline]
fn aggregate_counts_per_sample(
    per_sample_filtered_records: &FxHashMap<usize, Vec<&BedMethylLine>>,
    sample_index: &MultiSampleIndex,
) -> anyhow::Result<AggregatedCounts> {
    let combined_counts = per_sample_filtered_records
        .into_iter()
        .filter_map(|(sample, records)| {
            match aggregate_counts(&records, &sample_index.code_lookup) {
                Ok(counts) => Some(counts),
                Err(e) => {
                    debug!("sample {sample} failed to aggregate counts, {e}");
                    None
                }
            }
        })
        .collect::<Vec<AggregatedCounts>>();
    combined_counts
        .into_iter()
        .reduce(|a, b| a.op(&b))
        .ok_or_else(|| anyhow!("all samples failed"))
}

pub(super) fn get_modification_counts(
    sample_index: &MultiSampleIndex,
    dmr_batch: DmrBatch<Vec<RegionOfInterest>>,
) -> anyhow::Result<Vec<anyhow::Result<ModificationCounts>>> {
    // these are the bedmethyl records associated with the entire batch.
    // however, due to how tabix works, there will likely be additional
    // bedmethyl records that aren't part of any region, so we need to do
    // the filtering below.
    let (bedmethyl_lines_a, bedmethyl_lines_b) =
        sample_index.read_bedmethyl_group_by_chrom(&dmr_batch)?;

    let modification_counts_results = dmr_batch
        .dmr_chunks
        .into_par_iter()
        .map(|region_of_interest| {
            let filtered_a = filter_sample_records(
                &bedmethyl_lines_a,
                &region_of_interest,
                sample_index,
            );
            let filtered_b = filter_sample_records(
                &bedmethyl_lines_b,
                &region_of_interest,
                sample_index,
            );
            if filtered_a.is_empty() || filtered_b.is_empty() {
                let mut message = format!(
                    "missing bedMethy records for region {}, ",
                    &region_of_interest.dmr_interval
                );
                if filtered_a.is_empty() {
                    message.push_str("'a' has no records ");
                }
                if filtered_b.is_empty() {
                    message.push_str("'b' has no records")
                }
                bail!(message)
            } else {
                let control_counts =
                    aggregate_counts_per_sample(&filtered_a, &sample_index);
                let exp_counts =
                    aggregate_counts_per_sample(&filtered_b, &sample_index);
                match (control_counts, exp_counts) {
                    (Ok(control_counts), Ok(exp_counts)) => {
                        ModificationCounts::new(
                            control_counts,
                            exp_counts,
                            region_of_interest.dmr_interval,
                        )
                    }
                    (Err(e), Err(f)) => {
                        bail!(
                            "failed to aggregate control counts, {} and \
                             experimental counts, {}",
                            e.to_string(),
                            f.to_string()
                        )
                    }
                    (Err(e), _) => {
                        bail!(
                            "failed to aggregate control counts, {}",
                            e.to_string()
                        )
                    }
                    (_, Err(e)) => {
                        bail!(
                            "failed to aggregate experiment counts, {}",
                            e.to_string()
                        )
                    }
                }
            }
        })
        .collect::<Vec<Result<ModificationCounts, _>>>();

    Ok(modification_counts_results)
}

pub(super) fn run_pairwise_dmr(
    dmr_interval_iter: RoiIter,
    sample_index: Arc<MultiSampleIndex>,
    pool: rayon::ThreadPool,
    mut writer: Box<dyn std::io::Write>,
    pb: ProgressBar,
    header: bool,
    a_name: &str,
    b_name: &str,
    failure_counter: ProgressBar,
) -> anyhow::Result<usize> {
    if header {
        writer.write(ModificationCounts::header(a_name, b_name).as_bytes())?;
    }

    let (snd, rcv) = crossbeam_channel::bounded(1000);

    enum BatchResult {
        Results(Vec<anyhow::Result<ModificationCounts>>),
        BatchError(String, anyhow::Error, usize),
    }

    pool.spawn(move || {
        for batch in dmr_interval_iter {
            let batch_size = batch.dmr_chunks.len();
            let range_message = {
                let from = batch.dmr_chunks.iter().min_by(|a, b| a.cmp(&b));
                let to = batch.dmr_chunks.iter().max_by(|a, b| a.cmp(&b));
                match (from, to) {
                    (Some(s), Some(t)) => {
                        format!("{batch_size} intervals, from {} to {}", s, t)
                    }
                    _ => {
                        format!("{batch_size} intervals")
                    }
                }
            };
            match get_modification_counts(&sample_index, batch) {
                Ok(results) => {
                    let results = BatchResult::Results(results);
                    match snd.send(results) {
                        Ok(_) => {}
                        Err(e) => {
                            error!("failed to send results, {}", e.to_string())
                        }
                    }
                }
                Err(e) => {
                    let batch_error =
                        BatchResult::BatchError(range_message, e, batch_size);
                    match snd.send(batch_error) {
                        Ok(_) => {}
                        Err(e) => {
                            error!("failed to batch error, {}", e.to_string())
                        }
                    }
                }
            }
        }
    });

    let mut success_count = 0;
    for batch_result in rcv {
        match batch_result {
            BatchResult::Results(results) => {
                for result in results {
                    match result {
                        Ok(counts) => {
                            writer.write(counts.to_row()?.as_bytes())?;
                            success_count += 1;
                            pb.inc(1);
                        }
                        Err(e) => {
                            debug!("region failed, error: {}", e.to_string());
                            failure_counter.inc(1);
                        }
                    }
                }
            }
            BatchResult::BatchError(message, error, batch_size) => {
                debug!(
                    "failed entire dmr batch, {message}, {}",
                    error.to_string()
                );
                failure_counter.inc(batch_size as u64);
            }
        }
    }

    pb.finish_and_clear();

    Ok(success_count)
}
