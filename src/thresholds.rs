use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result as AnyhowResult};
use derive_new::new;
use itertools::Itertools;
use log::{debug, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read as _};

use crate::commands::SampleModBaseProbs;
use crate::errs::RunError;
use crate::filter_thresholds::FilterThresholds;
use crate::mod_bam::{BaseModCall, BaseModProbs, ModBaseInfo};
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;
use crate::reads_sampler::{
    get_sampled_read_ids_to_base_mod_calls, ReadIdsToBaseModCalls,
};
use crate::util;
use crate::util::{record_is_secondary, AlignedPairs, Region, Strand};

fn percentile_linear_interp(xs: &[f32], q: f32) -> AnyhowResult<f32> {
    if xs.len() < 2 {
        Err(anyhow!("not enough data points to calculate percentile",))
    } else {
        if q > 1.0 {
            return Err(anyhow!("quantile must be less than 1.0 got {q}"));
        }
        assert!(q <= 1.0);
        let l = xs.len() as f32;
        let left = (l * q).floor();
        let right = (l * q).ceil();
        let g = (l * q).fract();
        let y0 = xs[left as usize];
        let y1 = xs[right as usize];
        let y = y0 * (1f32 - g) + y1 * g;
        Ok(y)
    }
}

pub fn modbase_records<T: bam::Read>(
    records: bam::Records<T>,
) -> impl Iterator<Item = (ModBaseInfo, AlignedPairs, String)> + '_ {
    records
        // skip records that fail to parse htslib (todo this could be cleaned up)
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        .filter_map(|record| {
            ModBaseInfo::new_from_record(&record).ok().and_then(|mbi| {
                if mbi.is_empty() {
                    None
                } else {
                    let record_name = util::get_query_name_string(&record)
                        .unwrap_or("failed_utf".to_string());
                    let aligned_pairs =
                        util::get_aligned_pairs_forward(&record)
                            .filter_map(|ap| ap.ok())
                            .collect::<AlignedPairs>();
                    Some((mbi, aligned_pairs, record_name))
                }
            })
        })
}

// #[derive(new)]
// pub struct SampledBaseModProbs {
//     pub(crate) base_mod_calls: HashMap<DnaBase, Vec<BaseModCall>>,
//     pub(crate) records_to_bases: HashMap<String, HashSet<DnaBase>>,
// }
//
// impl SampledBaseModProbs {
//     fn calc_threshold_from_probs(
//         probs: &mut [f32],
//         filter_percentile: f32,
//     ) -> AnyhowResult<f32> {
//         probs.sort_by(|x, y| x.partial_cmp(y).unwrap());
//
//         percentile_linear_interp(&probs, filter_percentile)
//             .with_context(|| format!("didn't sample enough data, try a larger fraction of another seed"))
//     }
//
//     pub(crate) fn get_per_base_thresholds(
//         &self,
//         filter_percentile: f32,
//     ) -> AnyhowResult<FilterThresholds> {
//         let per_base_thresholds = self
//             .base_mod_calls
//             .iter()
//             .map(|(dna_base, base_mod_calls)| {
//                 let mut probs = base_mod_calls
//                     .iter()
//                     .filter_map(|bmc| match bmc {
//                         BaseModCall::Canonical(f) => Some(*f),
//                         BaseModCall::Modified(f, _) => Some(*f),
//                         _ => None,
//                     })
//                     .collect::<Vec<f32>>();
//                 match Self::calc_threshold_from_probs(
//                     &mut probs,
//                     filter_percentile,
//                 ) {
//                     Ok(thresh) => Ok((*dna_base, thresh)),
//                     Err(e) => Err(e),
//                 }
//             })
//             .collect::<AnyhowResult<HashMap<DnaBase, f32>>>()?;
//         Ok(FilterThresholds::new(0f32, per_base_thresholds))
//     }
// }

// impl Moniod for SampledBaseModProbs {
//     fn zero() -> Self {
//         Self {
//             base_mod_calls: HashMap::new(),
//             records_to_bases: HashMap::new(),
//         }
//     }
//
//     fn op(self, other: Self) -> Self {
//         let mut records_to_bases = self.records_to_bases;
//         for (record_id, bases) in other.records_to_bases {
//             records_to_bases
//                 .entry(record_id)
//                 .or_insert(HashSet::new())
//                 .extend(&mut bases.into_iter());
//         }
//
//         let mut base_mod_calls = self.base_mod_calls;
//         for (dna_base, mut calls) in other.base_mod_calls {
//             base_mod_calls
//                 .entry(dna_base)
//                 .or_insert(Vec::new())
//                 .append(&mut calls);
//         }
//         Self {
//             base_mod_calls,
//             records_to_bases,
//         }
//     }
//
//     fn op_mut(&mut self, other: Self) {
//         for (dna_base, mut calls) in other.base_mod_calls {
//             self.base_mod_calls
//                 .entry(dna_base)
//                 .or_insert(Vec::new())
//                 .append(&mut calls);
//         }
//         for (record_id, bases) in other.records_to_bases {
//             self.records_to_bases
//                 .entry(record_id)
//                 .or_insert(HashSet::new())
//                 .extend(&mut bases.into_iter());
//         }
//     }
//
//     fn len(&self) -> usize {
//         self.base_mod_calls
//             .values()
//             .map(|vs| vs.len())
//             .sum::<usize>()
//     }
// }

// pub struct ModbaseProbSampler {
//     record_sampler: RecordSampler,
// }

// impl ModbaseProbSampler {
//     pub fn new(
//         sample_frac: Option<f64>,
//         num_reads: Option<usize>,
//         seed: Option<u64>,
//     ) -> Self {
//         Self {
//             record_sampler: RecordSampler::new_from_options(
//                 sample_frac,
//                 num_reads,
//                 seed,
//             ),
//         }
//     }
//
//     #[inline]
//     fn filter_call_in_region(
//         ref_start: Option<u64>,
//         ref_end: Option<u64>,
//         aligned_pairs: &AlignedPairs,
//         query_pos: &usize,
//         base_mod_probs: &BaseModProbs,
//     ) -> Option<BaseModCall> {
//         let r_pos = aligned_pairs.get(query_pos);
//         match (ref_start, ref_end, r_pos) {
//             (Some(st), Some(ed), Some(r_pos)) => {
//                 if *r_pos >= st && *r_pos < ed {
//                     Some(base_mod_probs.base_mod_call())
//                 } else {
//                     None
//                 }
//             }
//             (Some(_), Some(_), None) => None,
//             _ => Some(base_mod_probs.base_mod_call()),
//         }
//     }
//
//     pub fn sample_modbase_probs<T: bam::Read>(
//         &mut self,
//         records: bam::Records<T>,
//         with_progress: bool,
//         ref_start: Option<u64>,
//         ref_end: Option<u64>,
//     ) -> AnyhowResult<SampledBaseModProbs> {
//         let spinner = if with_progress {
//             Some(self.record_sampler.get_progress_bar())
//         } else {
//             None
//         };
//         let mod_base_info_iter = modbase_records(records);
//         // mapping of each base to their probabilities, i.e. C <> [(p_m, p_h), ..., etc]
//         let mut can_base_probs = HashMap::new();
//         // mapping of record_ids to which canonical bases they have calls for, i.e. read_id <> {C, A}
//         let mut records_to_can_bases = HashMap::new();
//         for (modbase_info, aligned_pairs, record_name) in mod_base_info_iter {
//             match self.record_sampler.ask() {
//                 Indicator::Use => {
//                     if modbase_info.is_empty() {
//                         records_to_can_bases
//                             .insert(record_name, HashSet::new());
//                         continue;
//                     }
//                     for (raw_canonical_base, strand, seq_pos_base_mod_probs) in
//                         modbase_info.iter_seq_base_mod_probs()
//                     {
//                         let canonical_base =
//                             match (DnaBase::parse(*raw_canonical_base), strand)
//                             {
//                                 (Err(_), _) => continue,
//                                 (Ok(dna_base), Strand::Positive) => dna_base,
//                                 (Ok(dna_base), Strand::Negative) => {
//                                     dna_base.complement()
//                                 }
//                             };
//                         let mut mod_probs = seq_pos_base_mod_probs
//                             .pos_to_base_mod_probs
//                             .iter()
//                             .filter_map(|(q_pos, base_mod_probs)| {
//                                 Self::filter_call_in_region(
//                                     ref_start,
//                                     ref_end,
//                                     &aligned_pairs,
//                                     q_pos,
//                                     &base_mod_probs,
//                                 )
//                             })
//                             .collect::<Vec<BaseModCall>>();
//                         let probs = can_base_probs
//                             .entry(canonical_base)
//                             .or_insert(Vec::new());
//                         probs.append(&mut mod_probs);
//                         records_to_can_bases
//                             .entry(record_name.to_owned())
//                             .or_insert(HashSet::new())
//                             .insert(canonical_base);
//                     }
//                     if let Some(pb) = &spinner {
//                         pb.inc(1);
//                     }
//                 }
//                 Indicator::Skip => continue,
//                 Indicator::Done => break,
//             }
//         }
//         if let Some(pb) = &spinner {
//             pb.finish_and_clear();
//             info!("sampled {} records", records_to_can_bases.len());
//             for (can_base, probs) in can_base_probs.iter() {
//                 info!("Sampled {} {} calls", probs.len(), can_base.char(),);
//             }
//         }
//
//         Ok(SampledBaseModProbs::new(
//             can_base_probs,
//             records_to_can_bases,
//         ))
//     }
// }
//
// impl IntervalProcessor for ModbaseProbSampler {
//     type Output = SampledBaseModProbs;
//
//     fn new(
//         sample_frac: Option<f64>,
//         num_reads: Option<usize>,
//         seed: Option<u64>,
//     ) -> Self {
//         Self::new(sample_frac, num_reads, seed)
//     }
//
//     // todo rename process records in interval?
//     fn process_records<T: bam::Read>(
//         &mut self,
//         records: bam::Records<T>,
//         chrom_length: u32,
//         total_length: u64,
//         start: u64,
//         end: u64,
//     ) -> AnyhowResult<Self::Output> {
//         let num_reads = self.record_sampler.num_reads.map(|nr| {
//             let f = chrom_length as f64 / total_length as f64;
//             let nr = nr as f64 * f;
//             std::cmp::max(nr.floor() as usize, 1usize)
//         });
//         let sample_frac = self.record_sampler.sample_frac;
//         let seed = self.record_sampler.seed;
//         self.record_sampler =
//             RecordSampler::new_from_options(sample_frac, num_reads, seed);
//         self.sample_modbase_probs(records, false, Some(start), Some(end))
//     }
//
//     fn label() -> &'static str {
//         "base modification calls"
//     }
// }

pub struct Percentiles {
    qs: Vec<(f32, f32)>,
}

impl Percentiles {
    pub fn new(
        probs: &mut [f32],
        desired_percentiles: &[f32],
    ) -> AnyhowResult<Self> {
        probs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let qs = desired_percentiles
            .iter()
            .map(|q| percentile_linear_interp(&probs, *q).map(|p| (*q, p)))
            .collect::<Result<Vec<(f32, f32)>, _>>()?;
        Ok(Self { qs })
    }

    pub fn report(&self) -> String {
        let mut s = String::new();
        s.push_str("q\tp\n");
        for (q, p) in &self.qs {
            s.push_str(&format!("{:.2}\t{:.3}\n", q, p));
        }

        s
    }
}

pub fn calc_thresholds_per_base(
    read_ids_to_base_mod_calls: &ReadIdsToBaseModCalls,
    filter_percentile: f32,
    threads: usize,
    default_threshold: Option<f32>,
) -> AnyhowResult<FilterThresholds> {
    debug!("calculating per base thresholds");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    // todo thread pools should only be used at the commands.rs level
    pool.install(|| {
        let st = std::time::Instant::now();
        let mut probs_per_base = read_ids_to_base_mod_calls.probs_per_base();
        debug!("probs per base took {:?}s", st.elapsed().as_secs());
        let st = std::time::Instant::now();
        let filter_thresholds = probs_per_base
            .iter_mut()
            .map(|(canonical_base, probs)| {
                probs.par_sort_by(|a, b| a.partial_cmp(b).unwrap());
                percentile_linear_interp(&probs, filter_percentile)
                    .map(|t| (*canonical_base, t))
            })
            .collect::<AnyhowResult<HashMap<DnaBase, f32>>>()?;
        debug!("filter thresholds took {}s", st.elapsed().as_secs());
        info!("calculated thresholds:");
        for (dna_base, thresh) in filter_thresholds.iter() {
            info!("{}: {}", dna_base.char(), thresh);
        }
        Ok(FilterThresholds::new(
            default_threshold.unwrap_or(0f32),
            filter_thresholds,
        ))
    })
}

pub fn calc_threshold_from_bam(
    bam_fp: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    filter_percentile: f32,
    seed: Option<u64>,
    region: Option<&Region>,
) -> AnyhowResult<f32> {
    // todo implement per-base thresholds
    let can_base_probs = get_modbase_probs_from_bam(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
    )?;
    let mut probs = can_base_probs
        .values()
        .flatten()
        .map(|p| *p)
        .collect::<Vec<f32>>();
    probs.sort_by(|x, y| x.partial_cmp(y).unwrap());

    let threshold = percentile_linear_interp(&probs, filter_percentile)
        .with_context(|| format!("didn't sample enough data, try a larger fraction of another seed"))?;
    Ok(threshold)
}

pub fn get_modbase_probs_from_bam(
    bam_fp: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
) -> AnyhowResult<HashMap<DnaBase, Vec<f32>>> {
    get_sampled_read_ids_to_base_mod_calls(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
    )
    .map(|x| x.probs_per_base())
}
