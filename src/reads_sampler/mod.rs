pub(crate) mod record_sampler;
pub(crate) mod sampling_schedule;

use crate::interval_chunks::IntervalChunks;
use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    get_master_progress_bar, get_subroutine_progress_bar, get_targets,
    get_ticker, Region,
};
use anyhow::anyhow;
use indicatif::{MultiProgress, ParallelProgressIterator};
use log::debug;
use rayon::prelude::*;
use record_sampler::RecordSampler;
use rust_htslib::bam::{self, Read};
use std::path::PathBuf;

pub(crate) fn get_sampled_read_ids_to_base_mod_probs<P: RecordProcessor>(
    bam_fp: &PathBuf,
    reader_threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    position_filter: Option<&StrandedPositionFilter>,
    only_mapped: bool,
    suppress_progress: bool,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let use_regions = bam::IndexedReader::from_path(&bam_fp).is_ok();
    if use_regions {
        let schedule = match (sample_frac, num_reads) {
            (_, Some(num_reads)) => SamplingSchedule::from_num_reads(
                bam_fp,
                num_reads,
                region,
                position_filter,
                !only_mapped,
            ),
            (Some(frac), _) => SamplingSchedule::from_sample_frac(
                bam_fp,
                frac as f32,
                region,
                position_filter,
                !only_mapped,
            ),
            (None, None) => SamplingSchedule::from_sample_frac(
                bam_fp,
                1.0,
                region,
                position_filter,
                !only_mapped,
            ),
        }?;
        let mut read_ids_to_base_mod_calls =
            sample_reads_base_mod_calls_over_regions::<P>(
                bam_fp,
                interval_size,
                region,
                edge_filter,
                collapse_method,
                position_filter,
                &schedule,
                only_mapped,
                suppress_progress,
            )?;
        let should_sample_unmapped = schedule.unmapped_count > 0
            || read_ids_to_base_mod_calls.len() < 100;
        if should_sample_unmapped && !only_mapped {
            debug!(
                "sampled {} mapped records, sampling unmapped records",
                read_ids_to_base_mod_calls.len()
            );
            let mut reader = bam::IndexedReader::from_path(bam_fp)?;
            reader.set_threads(reader_threads)?;
            reader.fetch(bam::FetchDefinition::Unmapped)?;
            let num_reads_unmapped =
                num_reads.map(|nr| nr - read_ids_to_base_mod_calls.len());
            let record_sampler = RecordSampler::new_from_options(
                sample_frac,
                num_reads_unmapped,
                seed,
            );
            let unmapped_read_ids_to_base_mod_calls = P::process_records(
                reader.records(),
                !suppress_progress,
                record_sampler,
                collapse_method,
                edge_filter,
                position_filter,
                only_mapped,
            )?;
            debug!(
                "sampled {} unmapped records",
                unmapped_read_ids_to_base_mod_calls.len()
            );
            read_ids_to_base_mod_calls
                .op_mut(unmapped_read_ids_to_base_mod_calls);
        }
        debug!("sampled {} records", read_ids_to_base_mod_calls.len());

        Ok(read_ids_to_base_mod_calls)
    } else {
        if region.is_some() {
            return Err(anyhow!("cannot use region without indexed BAM"));
        }
        if position_filter.is_some() {
            debug!("using include-bed with an indexed bam would improve performance");
        }
        let mut reader = bam::Reader::from_path(bam_fp)?;
        reader.set_threads(reader_threads)?;
        let record_sampler =
            RecordSampler::new_from_options(sample_frac, num_reads, seed);
        let read_ids_to_base_mod_probs = P::process_records(
            reader.records(),
            !suppress_progress,
            record_sampler,
            collapse_method,
            edge_filter,
            position_filter,
            only_mapped,
        )?;
        debug!("sampled {} records", read_ids_to_base_mod_probs.len());
        Ok(read_ids_to_base_mod_probs)
    }
}

/// Sample reads evenly over a specified region or over
/// an entire sorted, aligned BAM.
fn sample_reads_base_mod_calls_over_regions<P: RecordProcessor>(
    bam_fp: &PathBuf,
    interval_size: u32,
    region: Option<&Region>,
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    position_filter: Option<&StrandedPositionFilter>,
    sampling_schedule: &SamplingSchedule,
    only_mapped: bool,
    suppress_progress: bool,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let header = reader.header();

    let contigs = get_targets(header, region);

    // prog bar stuff
    let master_progress = MultiProgress::new();
    if suppress_progress {
        master_progress
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    let tid_progress =
        master_progress.add(get_master_progress_bar(contigs.len()));
    tid_progress.set_message("contigs");
    let sampled_items = master_progress.add(get_ticker());
    sampled_items.set_message("base mod calls sampled");
    // end prog bar stuff

    let mut aggregator = <P::Output as Moniod>::zero();
    for reference_record in contigs {
        let intervals = IntervalChunks::new(
            reference_record.start,
            reference_record.length,
            interval_size,
            reference_record.tid,
            None,
        )
        .filter(|(start, end)| {
            position_filter
                .as_ref()
                .map(|pf| {
                    pf.overlaps_not_stranded(
                        reference_record.tid,
                        *start as u64,
                        *end as u64,
                    )
                })
                .unwrap_or(true)
        })
        .collect::<Vec<(u32, u32)>>();
        // make the number of reads (if given) proportional to the length
        // of this reference
        let num_reads_for_contig =
            sampling_schedule.get_num_reads(reference_record.tid);
        if num_reads_for_contig == 0 {
            continue;
        }
        let total_interval_length = intervals
            .iter()
            .map(|(start, end)| end.checked_sub(*start).unwrap_or(0))
            .sum::<u32>();
        // progress bar stuff
        let interval_progress =
            master_progress.add(get_subroutine_progress_bar(intervals.len()));
        interval_progress
            .set_message(format!("processing {}", &reference_record.name));
        // end progress bar stuff

        let proc_outputs = intervals
            .into_par_iter()
            .progress_with(interval_progress)
            .filter_map(|(start, end)| {
                let n_reads_for_interval = sampling_schedule
                    .get_num_reads_for_interval(
                        &reference_record,
                        total_interval_length,
                        start,
                        end,
                    );
                let record_sampler =
                    RecordSampler::new_num_reads(n_reads_for_interval);
                match sample_reads_from_interval::<P>(
                    bam_fp,
                    reference_record.tid,
                    start,
                    end,
                    record_sampler,
                    collapse_method,
                    edge_filter,
                    position_filter,
                    only_mapped,
                ) {
                    Ok(res) => {
                        let sampled_count = res.size();
                        sampled_items.inc(sampled_count);
                        Some(res)
                    }
                    Err(e) => {
                        debug!(
                            "reference {} for interval {} to {} failed {}",
                            &reference_record.name,
                            start,
                            end,
                            e.to_string()
                        );
                        None
                    }
                }
            })
            .reduce(|| <P::Output as Moniod>::zero(), |a, b| a.op(b));
        aggregator.op_mut(proc_outputs);
        tid_progress.inc(1);
    }

    tid_progress.finish_and_clear();
    let _ = master_progress.clear();
    Ok(aggregator)
}

pub(crate) fn sample_reads_from_interval<P: RecordProcessor>(
    bam_fp: &PathBuf,
    chrom_tid: u32,
    start: u32,
    end: u32,
    record_sampler: RecordSampler,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    position_filter: Option<&StrandedPositionFilter>,
    only_mapped: bool,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod,
{
    let mut bam_reader = bam::IndexedReader::from_path(bam_fp)?;
    bam_reader.fetch(bam::FetchDefinition::Region(
        chrom_tid as i32,
        start as i64,
        end as i64,
    ))?;

    P::process_records(
        bam_reader.records(),
        false,
        record_sampler,
        collapse_method,
        edge_filter,
        position_filter,
        only_mapped,
    )
}
