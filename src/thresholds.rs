use std::path::PathBuf;

use anyhow::{anyhow, Context, Result as AnyhowResult};
use log::info;
use rust_htslib::bam::{self, Read as _};

use crate::errs::RunError;
use crate::interval_processor::{
    run_sampled_region_processor, Indicator, IntervalProcessor, RecordSampler,
};
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::util::{record_is_secondary, Region};

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
) -> impl Iterator<Item = ModBaseInfo> + '_ {
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
                    Some(mbi)
                }
            })
        })
}

pub struct ModbaseProbSampler {
    record_sampler: RecordSampler,
}

impl ModbaseProbSampler {
    pub fn new(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self {
        Self {
            record_sampler: RecordSampler::new_from_options(
                sample_frac,
                num_reads,
                seed,
            ),
        }
    }

    pub fn sample_modbase_probs<T: bam::Read>(
        &mut self,
        records: bam::Records<T>,
        with_progress: bool,
    ) -> AnyhowResult<Vec<f32>> {
        let spinner = if with_progress {
            Some(self.record_sampler.get_progress_bar())
        } else {
            None
        };
        let mod_base_info_iter = modbase_records(records);
        let mut probs = Vec::new();
        let mut record_count = 0usize;
        for modbase_info in mod_base_info_iter {
            match self.record_sampler.ask() {
                Indicator::Use => {
                    for (_canonical_base, _strand, seq_pos_base_mod_probs) in
                        modbase_info.iter_seq_base_mod_probs()
                    {
                        let mut mod_probs = seq_pos_base_mod_probs
                            .pos_to_base_mod_probs
                            .iter()
                            .map(|(_pos, base_mod_probs)| match base_mod_probs
                                .base_mod_call()
                            {
                                BaseModCall::Modified(p, _) => p,
                                BaseModCall::Canonical(p) => p,
                                BaseModCall::Filtered => {
                                    unreachable!(
                                        "should not encounter filtered calls"
                                    )
                                }
                            })
                            .collect::<Vec<f32>>();
                        probs.append(&mut mod_probs);
                        if let Some(pb) = &spinner {
                            pb.inc(1);
                        }
                    }
                    record_count += 1;
                }
                Indicator::Skip => continue,
                Indicator::Done => break,
            }
        }
        if let Some(pb) = &spinner {
            pb.finish_and_clear();
            info!(
                "Sampled {} probabilities from {} records",
                probs.len(),
                record_count
            );
        }

        Ok(probs)
    }
}

impl IntervalProcessor for ModbaseProbSampler {
    type Output = f32;

    fn new(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self {
        Self::new(sample_frac, num_reads, seed)
    }

    fn process_records<T: bam::Read>(
        &mut self,
        records: bam::Records<T>,
        chrom_length: u32,
        total_length: u64,
    ) -> AnyhowResult<Vec<Self::Output>> {
        let num_reads = self.record_sampler.num_reads.map(|nr| {
            let f = chrom_length as f64 / total_length as f64;
            let nr = nr as f64 * f;
            std::cmp::max(nr.floor() as usize, 1usize)
        });
        let sample_frac = self.record_sampler.sample_frac;
        let seed = self.record_sampler.seed;
        self.record_sampler =
            RecordSampler::new_from_options(sample_frac, num_reads, seed);
        self.sample_modbase_probs(records, false)
    }

    fn label() -> &'static str {
        "probabilities"
    }
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
            s.push_str(&format!("{:.2}\t{:.2}\n", q, p));
        }

        s
    }
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
    let mut probs = get_modbase_probs_from_bam(
        bam_fp,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        seed,
        region,
    )?;
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
) -> AnyhowResult<Vec<f32>> {
    // todo make this return association of base to probs
    let use_regions = bam::IndexedReader::from_path(&bam_fp).is_ok();
    if use_regions {
        run_sampled_region_processor::<ModbaseProbSampler>(
            bam_fp,
            threads,
            interval_size,
            sample_frac,
            num_reads,
            seed,
            region,
        )
    } else {
        if region.is_some() {
            return Err(anyhow!("cannot use region without indexed BAM"));
        }
        let mut reader = bam::Reader::from_path(bam_fp).map_err(|e| {
            RunError::new_input_error(format!(
                "failed to open bam, {}",
                e.to_string()
            ))
        })?;
        reader
            .set_threads(threads)
            .map_err(|e| RunError::new_input_error(e.to_string()))?;
        let mut sampler = ModbaseProbSampler::new(sample_frac, num_reads, seed);
        sampler.sample_modbase_probs(reader.records(), true)
    }
}
