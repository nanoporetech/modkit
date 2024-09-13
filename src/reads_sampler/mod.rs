use std::path::PathBuf;

use anyhow::anyhow;
use indicatif::{MultiProgress, ProgressBar};
use itertools::Itertools;
use log::debug;
use prettytable::row;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;

use crate::interval_chunks::{ChromCoordinates, ReferenceIntervalsFeeder};
use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::sampling_schedule::{
    CountOrSample, SamplingSchedule,
};
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    get_master_progress_bar, get_targets, get_ticker, ReferenceRecord, Region,
};
use record_sampler::RecordSampler;

pub(crate) mod record_sampler;
pub(crate) mod sampling_schedule;

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
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    suppress_progress: bool,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let use_regions = bam::IndexedReader::from_path(&bam_fp).is_ok();
    if use_regions {
        debug!(
            "found BAM index, sampling reads in {interval_size} base pair \
             chunks"
        );
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
                (reader_threads as f32 * 1.5).floor() as usize,
                region,
                edge_filter,
                collapse_method,
                position_filter,
                &schedule,
                only_mapped,
                suppress_progress,
            )?;
        let should_sample_unmapped =
            schedule.has_unmapped() || read_ids_to_base_mod_calls.len() < 100;
        if should_sample_unmapped && !only_mapped {
            debug!(
                "sampled {} mapped records, sampling unmapped records",
                read_ids_to_base_mod_calls.len()
            );
            let mut reader = bam::IndexedReader::from_path(bam_fp)?;
            reader.set_threads(reader_threads)?;
            reader.fetch(bam::FetchDefinition::Unmapped)?;
            let num_reads_unmapped = num_reads.map(|nr| {
                nr.checked_sub(read_ids_to_base_mod_calls.len()).unwrap_or(0)
            });
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
                false,
                None,
                None,
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
        debug!("did not find index to modBAM");
        if region.is_some() {
            return Err(anyhow!("cannot use region without indexed BAM"));
        }
        if position_filter.is_some() {
            debug!(
                "using include-bed with an indexed bam would improve \
                 performance"
            );
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
            false,
            None,
            None,
        )?;
        debug!("sampled {} records", read_ids_to_base_mod_probs.len());
        Ok(read_ids_to_base_mod_probs)
    }
}

/// Sample reads evenly over a specified region or over
/// an entire sorted, aligned BAM. Only uses primary alignments
fn sample_reads_base_mod_calls_over_regions<P: RecordProcessor>(
    bam_fp: &PathBuf,
    interval_size: u32,
    batch_size: usize,
    region: Option<&Region>,
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    sampling_schedule: &SamplingSchedule,
    only_mapped: bool,
    suppress_progress: bool,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let header = reader.header();

    let contigs = get_targets(header, region)
        .into_iter()
        .filter(|reference_record| {
            sampling_schedule.chrom_has_reads(reference_record.tid)
        })
        .collect::<Vec<ReferenceRecord>>();

    let contig_sizes = contigs
        .iter()
        .map(|rec| (rec.tid, rec.length))
        .collect::<FxHashMap<u32, u32>>();

    let feeder = ReferenceIntervalsFeeder::new(
        contigs,
        batch_size,
        interval_size,
        false,
        None,
        None,
    )?;

    // prog bar stuff
    let master_progress = MultiProgress::new();
    if suppress_progress {
        master_progress
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    let tid_progress =
        master_progress.add(get_master_progress_bar(feeder.total_length()));
    tid_progress.set_message("total bp");

    let sampled_items = master_progress.add(get_ticker());
    sampled_items.set_message("base mod calls sampled");
    // end prog bar stuff

    let mut aggregator = <P::Output as Moniod>::zero();
    let mut reads_sampled_per_chr = FxHashMap::default();
    let feeder = feeder.map(|x| x.unwrap());
    for super_batch in feeder {
        let total_batch_length =
            super_batch.iter().map(|c| c.total_length()).sum::<u64>();
        let super_batch_with_counts = sampling_schedule
            .accumulate_sample_counts(
                super_batch,
                &contig_sizes,
                &reads_sampled_per_chr,
                batch_size,
            );
        let (super_batch_result, chrom_counts_for_batch) =
            super_batch_with_counts
                .into_par_iter()
                .map(|multi_coords| {
                    run_batch::<P>(
                        bam_fp,
                        multi_coords,
                        sampling_schedule,
                        collapse_method,
                        edge_filter,
                        position_filter,
                        only_mapped,
                        false,
                        None,
                        &sampled_items,
                    )
                })
                .reduce(
                    || (<P::Output as Moniod>::zero(), FxHashMap::zero()),
                    |(a, x), (b, y)| (a.op(b), x.op(y)),
                );
        tid_progress.inc(total_batch_length);
        aggregator.op_mut(super_batch_result);
        reads_sampled_per_chr.op_mut(chrom_counts_for_batch);
    }
    log_sampled_reads(&reads_sampled_per_chr);

    Ok(aggregator)
}

fn run_batch<P: RecordProcessor>(
    bam_fp: &PathBuf,
    batch: Vec<(ChromCoordinates, CountOrSample)>,
    sampling_schedule: &SamplingSchedule,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    allow_non_primary: bool,
    kmer_size: Option<usize>,
    sampled_items_counter: &ProgressBar,
) -> (P::Output, FxHashMap<u32, usize>)
where
    P::Output: Moniod,
{
    batch
        .into_par_iter()
        .filter(|(cc, _)| sampling_schedule.chrom_has_reads(cc.chrom_tid))
        .filter(|(cc, _)| {
            position_filter
                .map(|pf| {
                    pf.overlaps_not_stranded(
                        cc.chrom_tid,
                        cc.start_pos as u64,
                        cc.end_pos as u64,
                    )
                })
                .unwrap_or(true)
        })
        .filter_map(|(cc, counts_or_sample)| {
            let record_sampler = match counts_or_sample {
                CountOrSample::Count(x) => RecordSampler::new_num_reads(x),
                CountOrSample::Sample(x) => {
                    RecordSampler::new_sample_frac(x as f64, None)
                }
                CountOrSample::All => RecordSampler::new_passthrough(),
            };

            match sample_reads_from_interval::<P>(
                bam_fp,
                cc.chrom_tid,
                cc.start_pos,
                cc.end_pos,
                None,
                record_sampler,
                collapse_method,
                edge_filter,
                position_filter,
                only_mapped,
                allow_non_primary,
                kmer_size,
            ) {
                Ok(res) => {
                    sampled_items_counter.inc(res.size());
                    Some((res, cc.chrom_tid))
                }
                Err(e) => {
                    debug!(
                        "reference {} for interval {} to {} failed {}",
                        cc.chrom_tid,
                        cc.start_pos,
                        cc.end_pos,
                        e.to_string()
                    );
                    None
                }
            }
        })
        .fold(
            || (<P::Output as Moniod>::zero(), FxHashMap::zero()),
            |(agg, mut counter), (out, chrom_tid)| {
                *counter.entry(chrom_tid).or_insert(0usize) += out.len();
                (agg.op(out), counter)
            },
        )
        .reduce(
            || (<P::Output as Moniod>::zero(), FxHashMap::zero()),
            |(a, x), (b, y)| (a.op(b), x.op(y)),
        )
}

pub(crate) fn sample_reads_from_interval<P: RecordProcessor>(
    bam_fp: &PathBuf,
    chrom_tid: u32,
    start: u32,
    end: u32,
    prev_end: Option<u32>,
    record_sampler: RecordSampler,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    allow_non_primary: bool,
    kmer_size: Option<usize>,
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
        allow_non_primary,
        prev_end,
        kmer_size,
    )
}

fn log_sampled_reads(sampled_reads_per_chr: &FxHashMap<u32, usize>) {
    let mut tab = prettytable::Table::new();
    tab.set_format(*prettytable::format::consts::FORMAT_CLEAN);
    let mut total = 0usize;
    sampled_reads_per_chr.iter().sorted_by(|(_, x), (_, y)| y.cmp(x)).for_each(
        |(chr, count)| {
            tab.add_row(row![chr, count]);
            total += *count;
        },
    );

    tab.add_row(row!["total", total]);
    debug!("final mapped reads sampled:\n{tab}");
}
