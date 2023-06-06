use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCode, ParseChar};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::calc_threshold_from_bam;
use crate::util::Region;
use anyhow::{anyhow, bail, Context};
use log::{debug, info};
use std::collections::HashMap;
use std::path::PathBuf;

pub(crate) fn parse_per_mod_thresholds(
    raw_per_mod_thresholds: &[String],
) -> anyhow::Result<HashMap<ModCode, f32>> {
    let per_mod_thresholds = raw_per_mod_thresholds
        .iter()
        .map(|raw| parse_raw_threshold::<ModCode>(raw))
        .collect::<anyhow::Result<HashMap<ModCode, f32>>>()?;
    per_mod_thresholds.iter().for_each(|(mod_code, thresh)| {
        info!("using threshold {thresh} for mod-code {}", mod_code.char());
    });
    Ok(per_mod_thresholds)
}

pub(crate) fn parse_thresholds(
    raw_base_thresholds: &[String],
    per_mod_thresholds: Option<HashMap<ModCode, f32>>,
) -> anyhow::Result<MultipleThresholdModCaller> {
    let (default, per_base_thresholds) =
        parse_per_base_thresholds(raw_base_thresholds)?;
    Ok(MultipleThresholdModCaller::new(
        per_base_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        default.unwrap_or(0f32),
    ))
}

pub(crate) fn get_threshold_from_options(
    in_bam: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: usize,
    no_filtering: bool,
    filter_percentile: f32,
    seed: Option<u64>,
    region: Option<&Region>,
    per_mod_thresholds: Option<HashMap<ModCode, f32>>,
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    suppress_progress: bool,
) -> anyhow::Result<MultipleThresholdModCaller> {
    if no_filtering {
        info!("not performing filtering");
        return Ok(MultipleThresholdModCaller::new_passthrough());
    }
    let (sample_frac, num_reads) = match sample_frac {
        Some(f) => {
            let pct = f * 100f64;
            info!("sampling {pct}% of reads");
            (Some(f), None)
        }
        None => {
            info!("sampling {num_reads} reads from BAM");
            (None, Some(num_reads))
        }
    };
    let per_base_thresholds = calc_threshold_from_bam(
        in_bam,
        threads,
        interval_size,
        sample_frac,
        num_reads,
        filter_percentile,
        seed,
        region,
        edge_filter,
        collapse_method,
        suppress_progress,
    )?;

    for (dna_base, threshold) in per_base_thresholds.iter() {
        debug!(
            "estimated pass threshold {threshold} for primary sequence base {}",
            dna_base.char()
        );
    }

    Ok(MultipleThresholdModCaller::new(
        per_base_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        0f32,
    ))
}

fn parse_raw_threshold<T: ParseChar>(raw: &str) -> anyhow::Result<(T, f32)> {
    let parts = raw.split(':').collect::<Vec<&str>>();
    if parts.len() != 2 {
        bail!(
            "encountered illegal per-base threshold {raw}, should \
                be <base>:<threshold>, e.g. C:0.75"
        )
    }
    let raw_base = parts[0]
        .chars()
        .nth(0)
        .ok_or(anyhow!("failed to parse canonical base {}", &parts[0]))?;
    let base = T::parse_char(raw_base)
        .context(format!("failed to parse base {}", raw_base))?;
    let threshold_value = parts[1]
        .parse::<f32>()
        .context(format!("failed to parse threshold value {}", &parts[1]))?;
    Ok((base, threshold_value))
}

fn parse_per_base_thresholds(
    raw_thresholds: &[String],
) -> anyhow::Result<(Option<f32>, HashMap<DnaBase, f32>)> {
    if raw_thresholds.is_empty() {
        return Err(anyhow!("no thresholds provided"));
    }
    if raw_thresholds.len() == 1 {
        let raw = &raw_thresholds[0];
        if raw.contains(':') {
            let (dna_base, threshold) = parse_raw_threshold::<DnaBase>(raw)?;
            info!("using threshold {} for base {}", threshold, dna_base.char());
            let per_base_threshold = vec![(dna_base, threshold)]
                .into_iter()
                .collect::<HashMap<DnaBase, f32>>();
            Ok((None, per_base_threshold))
        } else {
            let default_threshold = raw.parse::<f32>().context(format!(
                "failed to parse user defined threshold {raw}"
            ))?;
            Ok((Some(default_threshold), HashMap::new()))
        }
    } else {
        let mut default: Option<f32> = None;
        let mut per_base_thresholds = HashMap::new();
        for raw_threshold in raw_thresholds {
            if raw_threshold.contains(':') {
                let (dna_base, threshold) =
                    parse_raw_threshold::<DnaBase>(raw_threshold)?;
                info!(
                    "using threshold {} for base {}",
                    threshold,
                    dna_base.char()
                );
                let repeated = per_base_thresholds.insert(dna_base, threshold);
                if repeated.is_some() {
                    bail!("repeated threshold for base {}", dna_base.char())
                }
            } else {
                if let Some(_) = default {
                    bail!("default threshold encountered more than once")
                }
                let default_threshold =
                    raw_threshold.parse::<f32>().context(format!(
                        "failed to parse default threshold {raw_threshold}"
                    ))?;
                info!("setting default threshold to {}", default_threshold);
                default = Some(default_threshold);
            }
        }
        Ok((default, per_base_thresholds))
    }
}
