pub(crate) mod record_sampler;

use crate::interval_chunks::IntervalChunks;
use crate::mod_bam::CollapseMethod;
use crate::monoid::Moniod;
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    get_master_progress_bar, get_spinner, get_subroutine_progress_bar,
    get_targets, Region,
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
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let use_regions = bam::IndexedReader::from_path(&bam_fp).is_ok();
    if use_regions {
        let mut read_ids_to_base_mod_calls =
            sample_reads_base_mod_calls_over_regions::<P>(
                bam_fp,
                interval_size,
                sample_frac,
                num_reads,
                seed,
                region,
                collapse_method,
            )?;
        // sample unmapped reads iff we've sampled less than 90% of the number we've wanted to get
        // or 0 (from sample_frac).
        let should_sample_unmapped = num_reads
            .map(|nr| {
                let f = (nr as f32 - read_ids_to_base_mod_calls.len() as f32)
                    / nr as f32;
                let f = 1f32 - f;
                f <= 0.9f32
            })
            .unwrap_or(read_ids_to_base_mod_calls.len() < 100);
        if should_sample_unmapped {
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
                true,
                record_sampler,
                collapse_method,
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
        let mut reader = bam::Reader::from_path(bam_fp)?;
        reader.set_threads(reader_threads)?;
        let record_sampler =
            RecordSampler::new_from_options(sample_frac, num_reads, seed);
        let read_ids_to_base_mod_probs = P::process_records(
            reader.records(),
            true,
            record_sampler,
            collapse_method,
        )?;
        debug!("sampled {} records", read_ids_to_base_mod_probs.len());
        Ok(read_ids_to_base_mod_probs)
    }
}

// pub(crate) fn process_records_in_bam<P: RecordProcessor>(
//     bam_fp: &PathBuf,
//     interval_size: u32,
//     regions: &[Region],
//     collapse_method: Option<&CollapseMethod>,
// ) -> anyhow::Result<P::Output> {
//     let use_regions = bam::IndexedReader::from_path(bam_fp).is_ok();
//     if !use_regions && !regions.is_empty() {
//         return
//     }
//
//     unimplemented!()
// }

/// Sample reads evenly over a specified region or over
/// an entire sorted, aligned BAM.
fn sample_reads_base_mod_calls_over_regions<P: RecordProcessor>(
    bam_fp: &PathBuf,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
    collapse_method: Option<&CollapseMethod>,
) -> anyhow::Result<P::Output>
where
    P::Output: Moniod + WithRecords,
{
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let header = reader.header();
    // could be regions, plural here
    let references = get_targets(header, region);
    let total_length = references.iter().map(|r| r.length as u64).sum::<u64>();

    // prog bar stuff
    let master_progress = MultiProgress::new();
    let tid_progress =
        master_progress.add(get_master_progress_bar(references.len()));
    tid_progress.set_message("contigs");
    let sampled_items = master_progress.add(get_spinner());
    sampled_items.set_message("base mod calls sampled");
    // end prog bar stuff

    let mut aggregator = <P::Output as Moniod>::zero();
    for reference in references {
        let intervals = IntervalChunks::new(
            reference.start,
            reference.length,
            interval_size,
            reference.tid,
            None,
        )
        .collect::<Vec<(u32, u32)>>();
        // make the number of reads (if given) proportional to the length
        // of this reference
        let num_reads_for_reference = num_reads.map(|nr| {
            let f = reference.length as f64 / total_length as f64;
            let nr = nr as f64 * f;
            std::cmp::max(nr.floor() as usize, 1usize)
        });

        // progress bar stuff
        let interval_progress =
            master_progress.add(get_subroutine_progress_bar(intervals.len()));
        interval_progress
            .set_message(format!("processing {}", &reference.name));
        // end progress bar stuff

        let proc_outputs = intervals
            .into_par_iter()
            .progress_with(interval_progress)
            .filter_map(|(start, end)| {
                let n_reads_for_interval = num_reads_for_reference.map(|nr| {
                    let f = (end - start) as f64 / reference.length as f64;
                    let nr = nr as f64 * f;
                    std::cmp::max(nr.floor() as usize, 1usize)
                });
                let record_sampler = RecordSampler::new_from_options(
                    sample_frac,
                    n_reads_for_interval,
                    seed,
                );
                match sample_reads_from_interval::<P>(
                    bam_fp,
                    reference.tid,
                    start,
                    end,
                    record_sampler,
                    collapse_method,
                ) {
                    Ok(res) => {
                        let sampled_count = res.size();
                        sampled_items.inc(sampled_count);
                        Some(res)
                    }
                    Err(e) => {
                        debug!(
                            "reference {} for interval {} to {} failed {}",
                            &reference.name,
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

    // sample_read_base_mod_calls(
    //     bam_reader.records(),
    //     false,
    //     record_sampler,
    //     collapse_method,
    // )
    P::process_records(
        bam_reader.records(),
        false,
        record_sampler,
        collapse_method,
    )
}

// fn sample_read_base_mod_calls<T: Read>(
//     records: bam::Records<T>,
//     with_progress: bool,
//     mut record_sampler: RecordSampler,
//     collapse_method: Option<&CollapseMethod>,
// ) -> anyhow::Result<ReadIdsToBaseModProbs> {
//     let spinner = if with_progress {
//         Some(record_sampler.get_progress_bar())
//     } else {
//         None
//     };
//     let mod_base_info_iter = filter_records_iter(records);
//     let mut read_ids_to_mod_base_probs = ReadIdsToBaseModProbs::zero();
//     for (record, mod_base_info) in mod_base_info_iter {
//         match record_sampler.ask() {
//             Indicator::Use => {
//                 let record_name = get_query_name_string(&record)
//                     .unwrap_or("FAILED_UTF_DECODE".to_string());
//                 if mod_base_info.is_empty() {
//                     // add count of unused/no calls
//                     read_ids_to_mod_base_probs
//                         .add_read_without_probs(&record_name);
//                     continue;
//                 }
//
//                 let (_, base_mod_probs_iter) =
//                     mod_base_info.into_iter_base_mod_probs();
//                 for (raw_canonical_base, strand, seq_pos_base_mod_probs) in
//                     base_mod_probs_iter
//                 {
//                     let canonical_base =
//                         match (DnaBase::parse(raw_canonical_base), strand) {
//                             (Err(_), _) => continue,
//                             (Ok(dna_base), Strand::Positive) => dna_base,
//                             (Ok(dna_base), Strand::Negative) => {
//                                 dna_base.complement()
//                             }
//                         };
//                     let mod_probs = seq_pos_base_mod_probs
//                         .pos_to_base_mod_probs
//                         .into_iter()
//                         .map(|(_q_pos, base_mod_probs)| {
//                             if let Some(method) = collapse_method {
//                                 base_mod_probs.into_collapsed(method)
//                             } else {
//                                 base_mod_probs
//                             }
//                         })
//                         .collect::<Vec<BaseModProbs>>();
//                     read_ids_to_mod_base_probs.add_mod_probs_for_read(
//                         &record_name,
//                         canonical_base,
//                         mod_probs,
//                     );
//                 }
//                 if let Some(pb) = &spinner {
//                     pb.inc(1);
//                 }
//             }
//             Indicator::Skip => continue,
//             Indicator::Done => break,
//         }
//     }
//
//     if let Some(pb) = &spinner {
//         pb.finish_and_clear();
//     }
//
//     Ok(read_ids_to_mod_base_probs)
// }
