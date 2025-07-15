use std::collections::HashMap;
use std::path::PathBuf;

use crate::errs::{MkError, MkResult};
use crate::mod_bam::{BaseModProbs, CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::Region;
use anyhow::{Context, Result as AnyhowResult};
use derive_new::new;
use log::{debug, info};
use rayon::prelude::*;
use sortedlist_rs::SortedList;

#[derive(new)]
pub(crate) struct BaseThresholdInfo {
    pub(crate) threshold: f32,
    pub(crate) pseudo_probs: Vec<f32>,
}

#[inline]
fn calc_mean(xs: &[f32]) -> f32 {
    assert!(!xs.is_empty());
    let tot = xs.iter().sum::<f32>();
    tot / (xs.len() as f32)
}

fn calc_std(xs: &[f32], mean: f32) -> f32 {
    let var =
        xs.iter().map(|x| (*x - mean).powi(2)).sum::<f32>() / (xs.len() as f32);
    var.sqrt()
}

pub(crate) trait Percenileable {
    fn nelems(&self) -> usize;
    fn get(&self, idx: usize) -> f32;
    fn last_elem(&self) -> f32;
    fn mean(&self) -> f32;
    fn std(&self, mean: f32) -> f32;
}

impl Percenileable for &[f32] {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx]
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self[self.len() - 1]
    }

    fn mean(&self) -> f32 {
        calc_mean(&self)
    }

    fn std(&self, mean: f32) -> f32 {
        calc_std(&self, mean)
    }
}

impl Percenileable for [f32] {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx]
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self[self.len() - 1]
    }

    fn mean(&self) -> f32 {
        calc_mean(&self)
    }

    fn std(&self, mean: f32) -> f32 {
        calc_std(&self, mean)
    }
}

impl Percenileable for Vec<f32> {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx]
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self[self.len() - 1]
    }

    fn mean(&self) -> f32 {
        calc_mean(&self)
    }

    fn std(&self, mean: f32) -> f32 {
        calc_std(&self, mean)
    }
}

impl Percenileable for &mut [f32] {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx]
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self[self.len() - 1]
    }

    fn mean(&self) -> f32 {
        calc_mean(self)
    }

    fn std(&self, mean: f32) -> f32 {
        calc_std(&self, mean)
    }
}

impl Percenileable for &mut Vec<f32> {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx]
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self[self.len() - 1]
    }

    #[inline]
    fn mean(&self) -> f32 {
        calc_mean(&self)
    }

    #[inline]
    fn std(&self, mean: f32) -> f32 {
        calc_std(&self, mean)
    }
}

impl Percenileable for &SortedList<BaseModProbs> {
    fn nelems(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, idx: usize) -> f32 {
        self[idx].argmax_base_mod_call().prob()
    }

    #[inline]
    fn last_elem(&self) -> f32 {
        self.last().unwrap().argmax_base_mod_call().prob()
    }

    #[inline]
    fn mean(&self) -> f32 {
        assert!(!self.is_empty());
        let tot =
            self.iter().map(|x| x.argmax_base_mod_call().prob()).sum::<f32>();
        tot / (self.len() as f32)
    }

    #[inline]
    fn std(&self, mean: f32) -> f32 {
        let var = self
            .iter()
            .map(|x| (x.argmax_base_mod_call().prob() - mean).powi(2))
            .sum::<f32>()
            / (self.len() as f32);
        var.sqrt()
    }
}

pub(crate) fn percentile_linear_interp<T: Percenileable + ?Sized>(
    xs: &T,
    q: f32,
) -> MkResult<f32> {
    if xs.nelems() < 2 {
        Err(MkError::PercentileNotEnoughDatapoints(xs.nelems()))
    } else {
        if q > 1.0 {
            return Err(MkError::PercentileInvalidQuantile(q));
        }
        if q == 1.0f32 {
            Ok(xs.last_elem())
        } else {
            assert!(q < 1.0);
            let l = (xs.nelems() - 1usize) as f32;
            let left = (l * q).floor();
            let right = (l * q).ceil() as usize;
            let g = (l * q).fract();
            let y0 = xs.get(left as usize);
            let y1 = xs.get(right);
            let y = y0 * (1f32 - g) + y1 * g;
            Ok(y)
        }
    }
}

// pub(crate) fn percentile_linear_interp<T: Index<usize, Output=f32> +
// ?Sized>(xs: &T, q: f32) -> MkResult<f32> {     if xs.len() < 2 {
//         Err(MkError::PercentileNotEnoughDatapoints(xs.len()))
//     } else {
//         if q > 1.0 {
//             return Err(MkError::PercentileInvalidQuantile(q));
//         }
//         if q == 1.0f32 {
//             Ok(xs[xs.len() - 1])
//         } else {
//             assert!(q < 1.0);
//             let l = (xs.len() - 1usize) as f32;
//             let left = (l * q).floor();
//             let right = (l * q).ceil() as usize;
//             let g = (l * q).fract();
//             let y0 = xs[left as usize];
//             let y1 = xs[right];
//             let y = y0 * (1f32 - g) + y1 * g;
//             Ok(y)
//         }
//     }
// }

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
    let (mut probs_per_base, explicit_canonical_probs) =
        read_ids_to_base_mod_calls.mle_probs_per_base();
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
                .map(|t| {
                    let min_canonical_prob = explicit_canonical_probs
                        .get(canonical_base)
                        .copied()
                        // qual 254
                        .unwrap_or(0.9941406f32);
                    if t < min_canonical_prob {
                        (*canonical_base, t)
                    } else {
                        debug!(
                            "estimated threshold too high {t}, using \
                             {min_canonical_prob}"
                        );
                        (*canonical_base, min_canonical_prob)
                    }
                })
        })
        .collect::<AnyhowResult<HashMap<DnaBase, f32>>>()?;
    debug!("filter thresholds took {}s", st.elapsed().as_secs());

    log_calculated_thresholds(&filter_thresholds);

    Ok(MultipleThresholdModCaller::new(
        filter_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        default_threshold.unwrap_or(0f32),
    ))
}

pub(crate) fn calc_threshold_from_bam(
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
    num_pseudos: Option<usize>,
) -> anyhow::Result<HashMap<DnaBase, BaseThresholdInfo>> {
    let (mut can_base_probs, min_explicit_can_probs) =
        get_modbase_probs_from_bam(
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
            mod_base_probs.par_sort_by(|a, b| a.partial_cmp(b).unwrap());
            #[cfg(debug_assertions)]
            {
                let is_sorted = mod_base_probs.is_sorted();
                assert!(is_sorted);
            }
            let threshold =
                percentile_linear_interp(&mod_base_probs, filter_percentile)?;
            let min_canonical_prob = min_explicit_can_probs
                .get(dna_base)
                .copied()
                .unwrap_or(0.9941406f32); // highest bin
            let threshold = if threshold < min_canonical_prob {
                threshold
            } else {
                debug!(
                    "estimated threshold too high {threshold} for {dna_base}, \
                     using {min_canonical_prob}"
                );
                min_canonical_prob
            };
            // these are safe because the percentile calc would fail if there
            // were 1 element or fewer
            let pseudo_probs = if let Some(num_pseudos) = num_pseudos {
                ndarray::Array1::linspace(0.1, 0.9, num_pseudos)
                    .to_vec()
                    .into_iter()
                    .map(|pct| {
                        percentile_linear_interp(&mod_base_probs, pct).unwrap()
                    })
                    .collect::<Vec<f32>>()
            } else {
                Vec::new()
            };
            Ok((*dna_base, BaseThresholdInfo::new(threshold, pseudo_probs)))
        })
        .collect()
}

/// Returns a tuple where the first element is of each primary sequence base to
/// a vector of call probabilities and the second item is a mapping of each
/// canonical base to the minimum explicit canonical probability observed
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

#[cfg(test)]
mod thresolds_tests {
    use super::percentile_linear_interp;

    #[test]
    fn test_thresholds_oob() {
        let xs = (0..10usize).map(|i| i as f32).collect::<Vec<f32>>();
        percentile_linear_interp(&xs, 0.95)
            .expect("should calculate percentile");
    }
}
