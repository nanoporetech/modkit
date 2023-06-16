use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result as AnyhowResult};

use log::{debug, info};
use rayon::prelude::*;

use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::Region;

fn percentile_linear_interp(xs: &[f32], q: f32) -> AnyhowResult<f32> {
    if xs.len() < 2 {
        Err(anyhow!(
            "not enough data points (got {}) to calculate percentile",
            xs.len()
        ))
    } else {
        if q > 1.0 {
            return Err(anyhow!("quantile must be less than 1.0 got {q}"));
        }
        if q == 1.0f32 {
            Ok(xs[xs.len() - 1])
        } else {
            assert!(q < 1.0);
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
}

pub struct Percentiles {
    pub(crate) qs: Vec<(f32, f32)>,
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
    read_ids_to_base_mod_calls: &ReadIdsToBaseModProbs,
    filter_percentile: f32,
    default_threshold: Option<f32>,
    per_mod_thresholds: Option<HashMap<ModCode, f32>>,
) -> AnyhowResult<MultipleThresholdModCaller> {
    debug!("calculating per base thresholds");
    let st = std::time::Instant::now();
    let mut probs_per_base = read_ids_to_base_mod_calls.mle_probs_per_base();
    debug!("probs per base took {:?}s", st.elapsed().as_secs());

    let st = std::time::Instant::now();
    let filter_thresholds = probs_per_base
        .iter_mut()
        .map(|(canonical_base, probs)| {
            probs.par_sort_by(|a, b| a.partial_cmp(b).unwrap());
            percentile_linear_interp(&probs, filter_percentile)
                .with_context(|| {
                    format!(
                        "failed to calculate threshold for base {}",
                        canonical_base.char()
                    )
                })
                .map(|t| (*canonical_base, t))
        })
        .collect::<AnyhowResult<HashMap<DnaBase, f32>>>()?;
    debug!("filter thresholds took {}s", st.elapsed().as_secs());

    let mut threshold_message = "calculated thresholds:".to_string();
    for (dna_base, thresh) in filter_thresholds.iter() {
        threshold_message.push_str(&format!(
            " {}: {}",
            dna_base.char(),
            thresh
        ));
    }
    info!("{threshold_message}");

    Ok(MultipleThresholdModCaller::new(
        filter_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        default_threshold.unwrap_or(0f32),
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
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    position_filter: Option<&StrandedPositionFilter>,
    only_mapped: bool,
    suppress_progress: bool,
) -> AnyhowResult<HashMap<DnaBase, f32>> {
    let mut can_base_probs = get_modbase_probs_from_bam(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
        collapse_method,
        edge_filter,
        position_filter,
        only_mapped,
        suppress_progress,
    )?;
    can_base_probs
        .iter_mut()
        .map(|(dna_base, mod_base_probs)| {
            mod_base_probs.par_sort_by(|x, y| x.partial_cmp(y).unwrap());
            let threshold =
                percentile_linear_interp(&mod_base_probs, filter_percentile)?;
            Ok((*dna_base, threshold))
        })
        .collect()
}

pub fn get_modbase_probs_from_bam(
    bam_fp: &PathBuf,
    threads: usize,
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
) -> AnyhowResult<HashMap<DnaBase, Vec<f32>>> {
    get_sampled_read_ids_to_base_mod_probs::<ReadIdsToBaseModProbs>(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
        collapse_method,
        edge_filter,
        position_filter,
        only_mapped,
        suppress_progress,
    )
    .map(|x| x.mle_probs_per_base())
}
