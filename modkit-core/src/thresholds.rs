use std::collections::HashMap;
use std::path::PathBuf;

use crate::errs::{MkError, MkResult};
use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs_with_reference;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::Region;
use anyhow::{Context, Result as AnyhowResult};
use log::{debug, info};
use rayon::prelude::*;

pub(crate) fn percentile_linear_interp(xs: &[f32], q: f32) -> MkResult<f32> {
    if xs.len() < 2 {
        Err(MkError::PercentileNotEnoughDatapoints(xs.len()))
    } else {
        if q > 1.0 {
            return Err(MkError::PercentileInvalidQuantile(q));
        }
        if q == 1.0f32 {
            Ok(xs[xs.len() - 1])
        } else {
            assert!(q < 1.0);
            let l = (xs.len() - 1usize) as f32;
            let left = (l * q).floor();
            let right = (l * q).ceil() as usize;
            let g = (l * q).fract();
            let y0 = xs[left as usize];
            let y1 = xs[right];
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
        probs.par_sort_by(|x, y| x.partial_cmp(y).unwrap());
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

pub(crate) fn log_calculated_thresholds(
    per_base_thresholds: &HashMap<DnaBase, f32>,
) {
    let mut threshold_message = "calculated thresholds:".to_string();
    for (dna_base, thresh) in per_base_thresholds.iter() {
        threshold_message.push_str(&format!(
            " {}: {}",
            dna_base.char(),
            thresh
        ));
    }
    info!("{threshold_message}");
}

pub fn calc_thresholds_per_base(
    read_ids_to_base_mod_calls: &ReadIdsToBaseModProbs,
    filter_percentile: f32,
    default_threshold: Option<f32>,
    per_mod_thresholds: Option<HashMap<ModCodeRepr, f32>>,
) -> AnyhowResult<MultipleThresholdModCaller> {
    debug!("calculating per base thresholds");
    let st = std::time::Instant::now();
    let (probs_per_base, explicit_canonical_probs) =
        read_ids_to_base_mod_calls.mle_probs_per_base();
    debug!("probs per base took {:?}s", st.elapsed().as_secs());

    let st = std::time::Instant::now();
    let filter_thresholds = calculate_threshold_with_fallback(
        probs_per_base,
        &explicit_canonical_probs,
        filter_percentile,
    )?;
    debug!("filter thresholds took {}s", st.elapsed().as_secs());

    log_calculated_thresholds(&filter_thresholds);

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
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    suppress_progress: bool,
) -> AnyhowResult<HashMap<DnaBase, f32>> {
    calc_threshold_from_bam_with_reference(
        bam_fp,
        None,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        filter_percentile,
        seed,
        region,
        edge_filter,
        collapse_method,
        position_filter,
        only_mapped,
        suppress_progress,
    )
}

pub(crate) fn calc_threshold_from_bam_with_reference(
    bam_fp: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    filter_percentile: f32,
    seed: Option<u64>,
    region: Option<&Region>,
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    suppress_progress: bool,
) -> AnyhowResult<HashMap<DnaBase, f32>> {
    let (can_base_probs, explicit_can_probs) =
        get_modbase_probs_from_bam_with_reference(
            bam_fp,
            reference_fasta,
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
    calculate_threshold_with_fallback(
        can_base_probs,
        &explicit_can_probs,
        filter_percentile,
    )
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
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
    suppress_progress: bool,
) -> AnyhowResult<(HashMap<DnaBase, Vec<f32>>, HashMap<DnaBase, f32>)> {
    get_modbase_probs_from_bam_with_reference(
        bam_fp,
        None,
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
}

pub(crate) fn get_modbase_probs_from_bam_with_reference(
    bam_fp: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    threads: usize,
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
) -> AnyhowResult<(HashMap<DnaBase, Vec<f32>>, HashMap<DnaBase, f32>)> {
    get_sampled_read_ids_to_base_mod_probs_with_reference::<
        ReadIdsToBaseModProbs,
    >(
        bam_fp,
        reference_fasta,
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

pub(crate) fn calculate_threshold_with_fallback(
    mut probs_per_base: HashMap<DnaBase, Vec<f32>>,
    max_explicit_canonical_prob: &HashMap<DnaBase, f32>,
    filter_percentile: f32,
) -> anyhow::Result<HashMap<DnaBase, f32>> {
    probs_per_base
        .iter_mut()
        .map(|(dna_base, mod_base_probs)| {
            mod_base_probs.par_sort_by(|x, y| x.partial_cmp(y).unwrap());
            let threshold =
                percentile_linear_interp(&mod_base_probs, filter_percentile)
                    .with_context(|| {
                        format!(
                            "failed to calculate threshold for base {}",
                            dna_base.char()
                        )
                    })?;
            let min_canonical_prob = max_explicit_canonical_prob
                .get(dna_base)
                .copied()
                .unwrap_or_else(|| {
                    debug!(
                        "no explicit unmodified probs collected for {dna_base}"
                    );
                    0.9941406f32
                });
            if threshold < min_canonical_prob {
                Ok((*dna_base, threshold))
            } else {
                debug!(
                    "estimated threshold too high {threshold}, using \
                     {min_canonical_prob} (maximum explicit canonical \
                     probability)"
                );
                Ok((*dna_base, min_canonical_prob))
            }
        })
        .collect()
}

#[cfg(test)]
mod thresolds_tests {
    use common_macros::hash_map;

    use crate::{
        mod_base_code::DnaBase, thresholds::calculate_threshold_with_fallback,
    };

    use super::percentile_linear_interp;

    #[test]
    fn test_thresholds_oob() {
        let xs = (0..10usize).map(|i| i as f32).collect::<Vec<f32>>();
        percentile_linear_interp(&xs, 0.95)
            .expect("should calculate percentile");
    }

    #[test]
    fn test_thresholds_fallback() {
        let explicit_max_probs = hash_map! {
            DnaBase::C => 0.95,
            DnaBase::A => 0.95,
        };

        let c_probs = vec![1.0f32; 100];
        let a_probs = c_probs.clone();

        let probs_per_base = hash_map! {
            DnaBase::C => c_probs,
            DnaBase::A => a_probs,
        };

        let thresholds = calculate_threshold_with_fallback(
            probs_per_base,
            &explicit_max_probs,
            0.1,
        )
        .expect("should calculate thresholds");

        for (_, t) in thresholds {
            assert_eq!(t, 0.95);
        }
    }
}
