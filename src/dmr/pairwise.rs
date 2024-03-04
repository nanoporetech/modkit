use std::sync::Arc;

use anyhow::bail;
use indicatif::ProgressBar;
use log::{debug, error};
use rayon::prelude::*;

use crate::dmr::bedmethyl::{aggregate_counts, BedMethylLine};
use crate::dmr::llr_model::ModificationCounts;
use crate::dmr::tabix::MultiSampleIndex;
use crate::dmr::util::{DmrBatch, RegionOfInterest, RoiIter};

pub(super) fn get_modification_counts(
    sample_index: &MultiSampleIndex,
    dmr_batch: DmrBatch<Vec<RegionOfInterest>>,
) -> anyhow::Result<Vec<anyhow::Result<ModificationCounts>>> {
    // these are the bedmethyl records associated with the entire batch.
    // however, due to how tabix works, there will likely be additional
    // bedmethyl records that aren't part of any region, so we need to do
    // the filtering below.
    let (bedmethyl_lines_a, bedmethyl_lines_b) =
        sample_index.read_bedmethyl_lines_collapse_on_chrom(&dmr_batch)?;

    let modification_counts_results = dmr_batch
        .dmr_chunks
        .into_par_iter()
        .map(|region_of_interest| {
            let filtered_a = bedmethyl_lines_a
                .get(&region_of_interest.dmr_interval.chrom)
                .map(|lines| {
                    lines
                        .iter()
                        .filter(|l| {
                            region_of_interest
                                .positions
                                .contains(&l.get_stranded_position())
                        })
                        .collect::<Vec<&BedMethylLine>>()
                })
                .unwrap_or_else(|| Vec::new());
            let filtered_b = bedmethyl_lines_b
                .get(&region_of_interest.dmr_interval.chrom)
                .map(|lines| {
                    lines
                        .iter()
                        .filter(|l| {
                            region_of_interest
                                .positions
                                .contains(&l.get_stranded_position())
                        })
                        .collect::<Vec<&BedMethylLine>>()
                })
                .unwrap_or_else(|| Vec::new());

            let control_counts = aggregate_counts(&filtered_a);
            let exp_counts = aggregate_counts(&filtered_b);
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
