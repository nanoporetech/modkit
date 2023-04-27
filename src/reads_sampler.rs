use crate::interval_chunks::IntervalChunks;
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;
use crate::record_sampler::{Indicator, RecordSampler};
use crate::util::{
    get_master_progress_bar, get_query_name_string, get_spinner,
    get_subroutine_progress_bar, get_targets, record_is_secondary, Region,
    Strand,
};
use anyhow::anyhow;
use indicatif::{MultiProgress, ParallelProgressIterator};
use log::{debug, error};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

/// Read IDs mapped to their base modification calls, organized by the
/// canonical base. This data structure contains essentially all of the
/// same data as in the records themselves, but with the query position
/// and the alternative probabilities removed (i.e. it only has the
/// probability of the called modification).
pub(crate) struct ReadIdsToBaseModCalls {
    // mapping of read id to canonical base mapped to a vec
    // of base mod calls on that canonical base
    pub(crate) inner: HashMap<String, HashMap<DnaBase, Vec<BaseModCall>>>,
}

impl ReadIdsToBaseModCalls {
    fn add_read_without_calls(&mut self, read_id: &str) {
        self.inner
            .entry(read_id.to_owned())
            .or_insert(HashMap::new());
    }

    fn add_mod_calls_for_read(
        &mut self,
        read_id: &str,
        canonical_base: DnaBase,
        mod_calls: Vec<BaseModCall>,
    ) {
        let read_id_entry = self
            .inner
            .entry(read_id.to_owned())
            .or_insert(HashMap::new());
        let added = read_id_entry.insert(canonical_base, mod_calls);
        if added.is_some() {
            error!("double added base mod calls, potentially a logic error!")
        }
    }

    fn size(&self) -> u64 {
        let s = self
            .inner
            .iter()
            .map(|(_, base_mod_calls)| {
                base_mod_calls.values().map(|vs| vs.len()).sum::<usize>()
            })
            .sum::<usize>();
        s as u64
    }

    pub(crate) fn num_reads(&self) -> usize {
        self.inner.len()
    }

    #[inline]
    pub(crate) fn probs_per_base(&self) -> HashMap<DnaBase, Vec<f32>> {
        let pb = get_master_progress_bar(self.inner.len());
        pb.set_message("aggregating per-base modification probabilities");
        self.inner
            .par_iter()
            .progress_with(pb)
            .map(|(_, canonical_base_to_base_mod_calls)| {
                canonical_base_to_base_mod_calls
                    .iter()
                    .map(|(canonical_base, base_mod_calls)| {
                        let probs = base_mod_calls
                            .iter()
                            .filter_map(|bmc| match bmc {
                                BaseModCall::Modified(f, _) => Some(*f),
                                BaseModCall::Canonical(f) => Some(*f),
                                _ => None,
                            })
                            .collect::<Vec<f32>>();
                        (*canonical_base, probs)
                    })
                    .collect::<HashMap<DnaBase, Vec<f32>>>()
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }

    pub(crate) fn probs_per_base_mod_call(&self) -> HashMap<char, Vec<f64>> {
        self.inner
            .par_iter()
            .filter_map(|(_, base_mod_calls)| {
                let grouped = base_mod_calls
                    .iter()
                    .map(|(base, base_mod_calls)| {
                        let canonical_code = base
                            .canonical_mod_code()
                            .map(|c| c.char())
                            .unwrap_or(base.char());
                        base_mod_calls
                            .iter()
                            .filter_map(|bmc| match bmc {
                                BaseModCall::Modified(p, code) => {
                                    Some((code.char(), *p))
                                }
                                BaseModCall::Canonical(p) => {
                                    Some((canonical_code, *p))
                                }
                                BaseModCall::Filtered => None,
                            })
                            .fold(
                                HashMap::<char, Vec<f64>>::new(),
                                |mut acc, (base, p)| {
                                    acc.entry(base)
                                        .or_insert(Vec::new())
                                        .push(p as f64);
                                    acc
                                },
                            )
                    })
                    .reduce(|a, b| a.op(b));
                grouped
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }
}

impl Moniod for ReadIdsToBaseModCalls {
    fn zero() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    fn op(self, other: Self) -> Self {
        let mut acc = self.inner;
        for (read_id, base_mod_calls) in other.inner {
            if acc.contains_key(&read_id) {
                continue;
            } else {
                acc.insert(read_id, base_mod_calls);
            }
        }

        Self { inner: acc }
    }

    fn op_mut(&mut self, other: Self) {
        for (read_id, base_mod_calls) in other.inner {
            if self.inner.contains_key(&read_id) {
                continue;
            } else {
                self.inner.insert(read_id, base_mod_calls);
            }
        }
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Entry point to sampling base modification calls from a BAM without regard
/// to their mapping.
pub(crate) fn get_sampled_read_ids_to_base_mod_calls(
    bam_fp: &PathBuf,
    reader_threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
) -> anyhow::Result<ReadIdsToBaseModCalls> {
    let use_regions = bam::IndexedReader::from_path(&bam_fp).is_ok();
    if use_regions {
        let mut read_ids_to_base_mod_calls =
            sample_reads_base_mod_calls_over_regions(
                bam_fp,
                interval_size,
                sample_frac,
                num_reads,
                seed,
                region,
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
            let unmapped_read_ids_to_base_mod_calls =
                sample_read_base_mod_calls(
                    reader.records(),
                    true,
                    record_sampler,
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
        let read_ids_to_base_mod_calls =
            sample_read_base_mod_calls(reader.records(), true, record_sampler)?;
        debug!("sampled {} records", read_ids_to_base_mod_calls.len());
        Ok(read_ids_to_base_mod_calls)
    }
}

/// Sample reads evenly over a specified region or over
/// an entire sorted, aligned BAM.
fn sample_reads_base_mod_calls_over_regions(
    bam_fp: &PathBuf,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
) -> anyhow::Result<ReadIdsToBaseModCalls> {
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

    let mut aggregator = ReadIdsToBaseModCalls::zero();
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

        let read_ids_to_mod_calls = intervals
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
                match sample_reads_from_interval(
                    bam_fp,
                    reference.tid,
                    start,
                    end,
                    record_sampler,
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
            .reduce(|| ReadIdsToBaseModCalls::zero(), |a, b| a.op(b));
        aggregator.op_mut(read_ids_to_mod_calls);
        tid_progress.inc(1);
    }

    tid_progress.finish_and_clear();
    Ok(aggregator)
}

fn filter_records_iter<T: Read>(
    records: bam::Records<T>,
) -> impl Iterator<Item = (bam::Record, ModBaseInfo)> + '_ {
    records
        // skip records that fail to parse htslib (todo this could be cleaned up)
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        .filter_map(|record| {
            ModBaseInfo::new_from_record(&record).ok().and_then(
                |mod_base_info| {
                    if mod_base_info.is_empty() {
                        None
                    } else {
                        Some((record, mod_base_info))
                    }
                },
            )
        })
}

fn sample_reads_from_interval(
    bam_fp: &PathBuf,
    chrom_tid: u32,
    start: u32,
    end: u32,
    record_sampler: RecordSampler,
) -> anyhow::Result<ReadIdsToBaseModCalls> {
    let mut bam_reader = bam::IndexedReader::from_path(bam_fp)?;
    bam_reader.fetch(bam::FetchDefinition::Region(
        chrom_tid as i32,
        start as i64,
        end as i64,
    ))?;

    sample_read_base_mod_calls(bam_reader.records(), false, record_sampler)
}

fn sample_read_base_mod_calls<T: Read>(
    records: bam::Records<T>,
    with_progress: bool,
    mut record_sampler: RecordSampler,
) -> anyhow::Result<ReadIdsToBaseModCalls> {
    let spinner = if with_progress {
        Some(record_sampler.get_progress_bar())
    } else {
        None
    };
    let mod_base_info_iter = filter_records_iter(records);
    let mut read_ids_to_mod_base_probs = ReadIdsToBaseModCalls::zero();
    let mut warned = HashSet::new();
    for (record, mod_base_info) in mod_base_info_iter {
        match record_sampler.ask() {
            Indicator::Use => {
                let record_name = get_query_name_string(&record)
                    .unwrap_or("FAILED_UTF_DECODE".to_string());
                if mod_base_info.is_empty() {
                    read_ids_to_mod_base_probs
                        .add_read_without_calls(&record_name);
                    continue;
                }

                for (raw_canonical_base, strand, seq_pos_base_mod_probs) in
                    mod_base_info.iter_seq_base_mod_probs()
                {
                    let canonical_base =
                        match (DnaBase::parse(*raw_canonical_base), strand) {
                            (Err(_), _) => continue,
                            (Ok(dna_base), Strand::Positive) => dna_base,
                            (Ok(dna_base), Strand::Negative) => {
                                dna_base.complement()
                            }
                        };
                    let mod_probs = seq_pos_base_mod_probs
                        .pos_to_base_mod_probs
                        .iter()
                        .filter_map(|(_q_pos, base_mod_probs)| {
                            match base_mod_probs.base_mod_call() {
                                Ok(base_mod_call) => Some(base_mod_call),
                                Err(e) => {
                                    let warning = e.to_string();
                                    match warned.contains(&warning) {
                                        true => {}
                                        false => {
                                            debug!("{warning}");
                                            warned.insert(warning);
                                        }
                                    }
                                    None
                                }
                            }
                        })
                        .collect::<Vec<BaseModCall>>();
                    read_ids_to_mod_base_probs.add_mod_calls_for_read(
                        &record_name,
                        canonical_base,
                        mod_probs,
                    );
                }
                if let Some(pb) = &spinner {
                    pb.inc(1);
                }
            }
            Indicator::Skip => continue,
            Indicator::Done => break,
        }
    }

    if let Some(pb) = &spinner {
        pb.finish_and_clear();
    }

    Ok(read_ids_to_mod_base_probs)
}
