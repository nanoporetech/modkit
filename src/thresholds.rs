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

pub(crate) fn calc_thresholds_per_base(
    read_ids_to_base_mod_calls: &ReadIdsToBaseModCalls,
    filter_percentile: f32,
    default_threshold: Option<f32>,
) -> AnyhowResult<FilterThresholds> {
    debug!("calculating per base thresholds");
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
