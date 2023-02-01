use crate::errs::{InputError, RunError};
use crate::mod_bam::{
    base_mod_probs_from_record, get_canonical_bases_with_mod_calls,
    BaseModCall, DeltaListConverter,
};
use crate::util::record_is_secondary;
use anyhow::Context;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rust_htslib::bam::{Read, Reader};
use std::path::Path;

fn percentile_linear_interp(xs: &[f32], q: f32) -> Result<f32, InputError> {
    if xs.len() < 2 {
        Err(InputError::new(
            "not enough data points to calculate percentile",
        ))
    } else {
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

pub fn sample_modbase_probs(
    reader: &mut Reader,
    seed: Option<u64>,
    frac: f64,
) -> Result<Vec<f32>, RunError> {
    let mut rng = if let Some(seed) = seed {
        SeedableRng::seed_from_u64(seed)
    } else {
        StdRng::from_entropy()
    };

    let record_iter = reader
        .records()
        // skip records that fail to parse htslib
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        // sample
        .filter_map(|record| {
            if rng.gen_bool(frac) {
                Some(record)
            } else {
                None
            }
        })
        // pull out the canonical bases in the MM tags, drop records that fail to parse
        .filter_map(|record| {
            get_canonical_bases_with_mod_calls(&record)
                .map(|bases| (bases, record))
                .ok()
        });

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} [{elapsed_precise}] {pos} {msg}",
        )
        .unwrap()
        .tick_strings(&[
            "▹▹▹▹▹",
            "▸▹▹▹▹",
            "▹▸▹▹▹",
            "▹▹▸▹▹",
            "▹▹▹▸▹",
            "▹▹▹▹▸",
            "▪▪▪▪▪",
        ]),
    );
    spinner.set_message("records sampled");

    let mut probs = Vec::new();
    let mut record_count = 0usize;
    for (canonical_bases, record) in record_iter {
        for canonical_base in canonical_bases {
            let converter = DeltaListConverter::new_from_record(
                &record,
                canonical_base.char(),
            )?;
            let seq_pos_base_mod_probs =
                base_mod_probs_from_record(&record, &converter)?;
            let mut mod_probs = seq_pos_base_mod_probs
                .iter()
                .map(|(_pos, base_mod_probs)| {
                    match base_mod_probs.base_mod_call() {
                        BaseModCall::Modified(p, _) => p,
                        BaseModCall::Canonical(p) => p,
                        BaseModCall::Filtered => {
                            panic!("should not encounter filtered calls")
                        }
                    }
                })
                .collect::<Vec<f32>>();
            probs.append(&mut mod_probs);
        }
        spinner.inc(1);
        record_count += 1;
    }
    spinner.finish_and_clear();
    info!(
        "Sampled {} probabilities from {} records",
        probs.len(),
        record_count
    );
    Ok(probs)
}

pub struct Percentiles {
    qs: Vec<(f32, f32)>,
}

impl Percentiles {
    pub fn new(
        probs: &mut [f32],
        desired_percentiles: &[f32],
    ) -> Result<Self, InputError> {
        probs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let qs = desired_percentiles
            .iter()
            .map(|q| percentile_linear_interp(&probs, *q).map(|p| (*q, p)))
            .collect::<Result<Vec<(f32, f32)>, InputError>>()?;
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

pub fn calc_threshold_from_bam<T: AsRef<Path>>(
    bam_fp: &T,
    threads: usize,
    frac: f64,
    filter_percentile: f32,
    seed: Option<u64>,
) -> anyhow::Result<f32, InputError> {
    let mut bam = Reader::from_path(bam_fp).map_err(|e| e.to_string())?;
    bam.set_threads(threads).map_err(|e| e.to_string())?;

    let mut probs = sample_modbase_probs(&mut bam, seed, frac)
        .map_err(|e| e.to_string())?;
    probs.sort_by(|x, y| x.partial_cmp(y).unwrap());

    let threshold = percentile_linear_interp(&probs, filter_percentile)
        .with_context(|| format!("didn't sample enough data, try a larger fraction of another seed"))
        .map_err(|e| e.to_string())?;

    Ok(threshold)
}
