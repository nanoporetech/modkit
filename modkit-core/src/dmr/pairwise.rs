use std::sync::Arc;

use crate::dmr::bedmethyl::{aggregate_counts, BedMethylLine};
use crate::dmr::llr_model::{AggregatedCounts, ModificationCounts};
use crate::dmr::tabix::{ChromToSampleBMLines, MultiSampleIndex};
use crate::dmr::util::{DmrBatch, RegionOfInterest, RoiIter};
use crate::errs::{MkError, MkResult};
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::BorrowingMoniod;
use indicatif::{MultiProgress, ProgressBar};
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
                .filter_map(|(sample, lines)| {
                    let overlapping_records = lines
                        .iter()
                        .filter(|record| {
                            roi.positions.contains(
                                &record.get_stranded_position(
                                    &sample_index.code_lookup,
                                ),
                            )
                        })
                        .collect::<Vec<&BedMethylLine>>();
                    if overlapping_records.is_empty() {
                        None
                    } else {
                        Some((*sample, overlapping_records))
                    }
                })
                .collect::<FxHashMap<usize, Vec<&BedMethylLine>>>()
        })
        .unwrap_or_else(|| FxHashMap::default())
}

#[inline]
fn aggregate_counts_per_sample(
    per_sample_filtered_records: &FxHashMap<usize, Vec<&BedMethylLine>>,
    sample_index: &MultiSampleIndex,
    single_code_op: Option<ModCodeRepr>,
) -> MkResult<AggregatedCounts> {
    // per_sample_filtered_records should always have non-zero length vectors
    let combined_counts = per_sample_filtered_records
        .values()
        .map(|records| {
            aggregate_counts(
                &records,
                &sample_index.code_lookup,
                single_code_op,
            )
        })
        .collect::<MkResult<Vec<AggregatedCounts>>>()?;
    combined_counts.into_iter().reduce(|a, b| a.op(&b)).ok_or_else(|| {
        // shouldn't really ever happen?
        debug!("all samples failed.. check the logs");
        MkError::DmrMissing
    })
}

/// Return type here is a little complicated:
/// The outermost result is to capture failure to read/IO the bedmethyl tables.
/// The results in the vector are due to the test of individual DMR intervals.
/// In general, if the outermost Result is Err, the whole program should fail.
/// It is up to the caller to decide whether to fail on the innermost Errs
pub(super) fn get_modification_counts(
    sample_index: &MultiSampleIndex,
    dmr_batch: DmrBatch<Vec<RegionOfInterest>>,
    single_code_op: Option<ModCodeRepr>,
) -> MkResult<Vec<Result<ModificationCounts, (MkError, Option<MkError>)>>> {
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
                debug!("{message}");
                Err((MkError::DmrMissing, None))
            } else {
                let control_counts = aggregate_counts_per_sample(
                    &filtered_a,
                    &sample_index,
                    single_code_op,
                );
                let exp_counts = aggregate_counts_per_sample(
                    &filtered_b,
                    &sample_index,
                    single_code_op,
                );
                match (control_counts, exp_counts) {
                    (Ok(control_counts), Ok(exp_counts)) => {
                        ModificationCounts::new(
                            control_counts,
                            exp_counts,
                            region_of_interest.dmr_interval,
                        )
                        .map_err(|e| (e, None))
                    }
                    (Err(e), Err(f)) => {
                        debug!(
                            "{}: failed to aggregate control counts, {} and \
                             experimental counts, {}",
                            region_of_interest.dmr_interval,
                            e.to_string(),
                            f.to_string()
                        );
                        Err((e, Some(f)))
                    }
                    (Err(e), _) => {
                        debug!(
                            "{}: failed to aggregate control counts, {}",
                            region_of_interest.dmr_interval,
                            e.to_string()
                        );
                        Err((e, None))
                    }
                    (_, Err(e)) => {
                        debug!(
                            "{}: failed to aggregate experiment counts, {}",
                            region_of_interest.dmr_interval,
                            e.to_string()
                        );
                        Err((e, None))
                    }
                }
            }
        })
        .collect::<Vec<Result<ModificationCounts, (MkError, Option<MkError>)>>>(
        );

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
    single_code_op: Option<ModCodeRepr>,
    failure_counter: ProgressBar,
    batch_failures: ProgressBar,
    multi_progress: MultiProgress,
) -> anyhow::Result<(usize, FxHashMap<String, usize>)> {
    if header {
        writer.write(ModificationCounts::header(a_name, b_name).as_bytes())?;
    }

    let (snd, rcv) = crossbeam_channel::bounded(1000);

    enum BatchResult {
        Results(Vec<Result<ModificationCounts, (MkError, Option<MkError>)>>),
        BatchError(String, MkError),
    }

    let pb_handle = multi_progress.clone();
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
            match get_modification_counts(&sample_index, batch, single_code_op)
            {
                Ok(results) => {
                    let results = BatchResult::Results(results);
                    match snd.send(results) {
                        Ok(_) => {}
                        Err(e) => {
                            pb_handle.suspend(|| {
                                error!(
                                    "failed to send results, {}",
                                    e.to_string()
                                );
                            });
                            break;
                        }
                    }
                }
                Err(e) => {
                    let batch_error = BatchResult::BatchError(range_message, e);
                    match snd.send(batch_error) {
                        Ok(_) => {}
                        Err(e) => {
                            pb_handle.suspend(|| {
                                error!(
                                    "failed to send batch error, {}",
                                    e.to_string()
                                );
                            });
                            break;
                        }
                    }
                }
            }
        }
    });

    let mut success_count = 0;
    let mut region_error_counts = FxHashMap::<String, usize>::default();
    let mut err: Option<MkError> = None;
    'rcv_loop: for batch_result in rcv {
        match batch_result {
            BatchResult::Results(results) => {
                for result in results {
                    match result {
                        Ok(counts) => {
                            writer.write(counts.to_row()?.as_bytes())?;
                            success_count += 1;
                            pb.inc(1);
                        }
                        Err((e, f)) => {
                            match (&e, f.as_ref()) {
                                (MkError::InvalidBedMethyl(message), _)
                                | (
                                    _,
                                    Some(MkError::InvalidBedMethyl(message)),
                                ) => {
                                    multi_progress.suspend(|| {
                                        error!(
                                            "encountered invalid bedMethyl \
                                             record(s), {message}, stopping"
                                        );
                                    });
                                    err = Some(e);
                                    break 'rcv_loop;
                                }
                                _ => {}
                            };
                            region_error_counts
                                .entry(e.to_string())
                                .and_modify(|e| *e = e.saturating_add(1))
                                .or_insert(1);
                            failure_counter.inc(1);
                        }
                    }
                }
            }
            BatchResult::BatchError(message, error) => {
                multi_progress.suspend(|| {
                    if let MkError::InvalidBedMethyl(m) = &error {
                        error!("failed entire dmr batch, {message}, {m}",);
                    } else {
                        error!("failed entire dmr batch, {message}, {error}",);
                    }
                });
                batch_failures.inc(1u64);
                err = Some(error);
                break 'rcv_loop;
            }
        }
    }

    pb.finish_and_clear();

    if let Some(e) = err {
        Err(e.into())
    } else {
        Ok((success_count, region_error_counts))
    }
}
