use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{anyhow, bail, Context};
use itertools::Itertools;
use log::{debug, info, warn};
use rust_htslib::bam::{self, Header};

use crate::adjust::OverlappingRegexOffset;
use crate::fasta::MotifLocationsLookup;
use crate::mod_bam::{prob_to_qual, CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::motifs::motif_bed::RegexMotif;
use crate::position_filter::StrandedPositionFilter;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::calc_threshold_from_bam;
use crate::util::{create_out_directory, Region};

pub(crate) fn parse_sampling_fraction(raw: &str) -> Result<f64, String> {
    let fraction = raw.parse::<f64>().map_err(|_| {
        format!(
            "sampling fraction must be a finite number in the inclusive \
             range [0, 1]; got '{raw}'"
        )
    })?;
    if fraction.is_finite() && (0.0..=1.0).contains(&fraction) {
        Ok(fraction)
    } else {
        Err(format!(
            "sampling fraction must be a finite number in the inclusive \
             range [0, 1]; got '{raw}'"
        ))
    }
}

pub fn parse_per_mod_thresholds(
    raw_per_mod_thresholds: &[String],
) -> anyhow::Result<HashMap<ModCodeRepr, f32>> {
    let per_mod_thresholds = raw_per_mod_thresholds
        .iter()
        .map(|raw| {
            let parts = raw.split(":").collect::<Vec<&str>>();
            if parts.len() != 2 {
                Err(anyhow!(
                    "encountered illegal per-mod threshold: {raw}. Should be \
                     mod_code:threshold e.g. h:0.8"
                ))
            } else {
                ModCodeRepr::parse(parts[0]).and_then(|x| {
                    parts[1].parse::<f32>().map(|t| (x, t)).map_err(|e| {
                        anyhow!(
                            "failed to parse per-mod threshold value {}, {}",
                            &parts[1],
                            e.to_string()
                        )
                    })
                })
            }
        })
        .collect::<anyhow::Result<HashMap<ModCodeRepr, f32>>>()?;
    per_mod_thresholds.iter().for_each(|(mod_code, thresh)| {
        info!("parsed user-input threshold {thresh} for mod-code {}", mod_code);
    });
    Ok(per_mod_thresholds)
}

pub fn parse_thresholds_values(
    raw_base_thresholds: &[String],
) -> anyhow::Result<(f32, [f32; 4])> {
    let (default_threshold, per_base_thresholds) =
        parse_per_base_thresholds(raw_base_thresholds)?;
    if !per_base_thresholds.is_empty() {
        for (base, threshold) in per_base_thresholds.iter() {
            info!(
                "using user-specified filter threshold {threshold} (qual: {}) \
                 for {base}",
                prob_to_qual(*threshold)
            );
        }
    }
    if let Some(default_threshold) = default_threshold {
        info!(
            "user-specified pass threshold {default_threshold} (qual: {})",
            prob_to_qual(default_threshold)
        );
    } else {
        let bases_with_thresholds = per_base_thresholds
            .keys()
            .map(|x| format!("{}", x.char()))
            .join(",");
        info!(
            "no default pass threshold was provided, so base modifications at \
             primary sequence bases other than {bases_with_thresholds} will \
             not be filtered"
        );
    }

    let base_thresholds = [
        *per_base_thresholds
            .get(&DnaBase::A)
            .or(default_threshold.as_ref())
            .unwrap_or(&0f32),
        *per_base_thresholds
            .get(&DnaBase::C)
            .or(default_threshold.as_ref())
            .unwrap_or(&0f32),
        *per_base_thresholds
            .get(&DnaBase::G)
            .or(default_threshold.as_ref())
            .unwrap_or(&0f32),
        *per_base_thresholds
            .get(&DnaBase::T)
            .or(default_threshold.as_ref())
            .unwrap_or(&0f32),
    ];

    Ok((default_threshold.unwrap_or(0f32), base_thresholds))
}

pub fn parse_thresholds(
    raw_base_thresholds: &[String],
    per_mod_thresholds: Option<HashMap<ModCodeRepr, f32>>,
) -> anyhow::Result<MultipleThresholdModCaller> {
    let (default, per_base_thresholds) =
        parse_per_base_thresholds(raw_base_thresholds)?;
    if default.is_none() {
        let bases_with_thresholds = per_base_thresholds
            .keys()
            .map(|x| format!("{}", x.char()))
            .join(",");
        info!(
            "no default pass threshold was provided, so base modifications at \
             primary sequence bases other than {bases_with_thresholds} will \
             not be filtered"
        );
    }

    Ok(MultipleThresholdModCaller::new(
        per_base_thresholds,
        per_mod_thresholds.unwrap_or_else(HashMap::new),
        default.unwrap_or(0f32),
    ))
}

pub fn get_threshold_from_options(
    in_bam: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: usize,
    no_filtering: bool,
    filter_percentile: f32,
    seed: Option<u64>,
    region: Option<&Region>,
    per_mod_thresholds: Option<HashMap<ModCodeRepr, f32>>,
    edge_filter: Option<&EdgeFilter>,
    collapse_method: Option<&CollapseMethod>,
    position_filter: Option<&StrandedPositionFilter<()>>,
    only_mapped: bool,
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
            info!("attempting to sample {num_reads} reads");
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
        position_filter,
        only_mapped,
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

fn parse_raw_threshold(raw: &str) -> anyhow::Result<(DnaBase, f32)> {
    let parts = raw.split(':').collect::<Vec<&str>>();
    if parts.len() != 2 {
        bail!(
            "encountered illegal per-base threshold {raw}, should be \
             <base>:<threshold>, e.g. C:0.75"
        )
    }
    let raw_base = parts[0]
        .chars()
        .nth(0)
        .ok_or(anyhow!("failed to parse canonical base {}", &parts[0]))?;
    let base = DnaBase::parse(raw_base)
        .context(format!("failed to parse base {}", raw_base))?;
    let threshold_value = parts[1]
        .parse::<f32>()
        .context(format!("failed to parse threshold value {}", &parts[1]))?;
    Ok((base, threshold_value))
}

pub(crate) fn parse_per_base_thresholds(
    raw_thresholds: &[String],
) -> anyhow::Result<(Option<f32>, HashMap<DnaBase, f32>)> {
    if raw_thresholds.is_empty() {
        bail!("no thresholds provided")
    }
    if raw_thresholds.len() == 1 {
        let raw = &raw_thresholds[0];
        if raw.contains(':') {
            let (dna_base, threshold) = parse_raw_threshold(raw)?;
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
                let (dna_base, threshold) = parse_raw_threshold(raw_threshold)?;
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
                default = Some(default_threshold);
            }
        }
        Ok((default, per_base_thresholds))
    }
}

pub(crate) fn parse_raw_thresholds_string_with_default(
    raw_thresholds: &[String],
) -> anyhow::Result<[f32; 4]> {
    parse_per_base_thresholds(raw_thresholds).map(|(x, per_base)| {
        let default = x.unwrap_or(0f32);
        [
            per_base.get(&DnaBase::A).copied().unwrap_or(default),
            per_base.get(&DnaBase::C).copied().unwrap_or(default),
            per_base.get(&DnaBase::G).copied().unwrap_or(default),
            per_base.get(&DnaBase::T).copied().unwrap_or(default),
        ]
    })
}

pub fn using_stream(raw: &str) -> bool {
    raw == "-" || raw == "stdin" || raw == "stdout"
}

pub fn get_serial_reader(
    raw: &str,
) -> rust_htslib::errors::Result<bam::Reader> {
    if using_stream(raw) {
        bam::Reader::from_stdin()
    } else {
        bam::Reader::from_path(raw)
    }
}

pub fn get_bam_writer(
    raw: &str,
    header: &Header,
    output_sam: bool,
) -> anyhow::Result<bam::Writer> {
    let format = if output_sam { bam::Format::Sam } else { bam::Format::Bam };
    if using_stream(raw) {
        bam::Writer::from_stdout(&header, format).map_err(|e| {
            anyhow!(
                "failed to make stdout {format:?} writer, {}",
                e.to_string()
            )
        })
    } else {
        create_out_directory(&raw)?;
        bam::Writer::from_path(&raw, &header, format).map_err(|e| {
            anyhow!("failed to make {format:?} writer, {}", e.to_string())
        })
    }
}

pub fn parse_edge_filter_input(
    raw: &str,
    inverted: bool,
) -> anyhow::Result<EdgeFilter> {
    if raw.contains(',') {
        let parts = raw.split(',').collect::<Vec<&str>>();
        if parts.len() != 2 {
            bail!(
                "illegal edge filter input {raw}, should be \
                 start_trim,end_trim (e.g. 4,5)"
            )
        }
        let start_trim = parts[0].parse::<usize>().context(format!(
            "failed to parse edge filter start trim {raw}, should be a number"
        ))?;
        let end_trim = parts[1].parse::<usize>().context(format!(
            "failed to parse edge filter end trim {raw}, should be a number"
        ))?;
        info!(
            "filtering out base modification calls {start_trim} bases from \
             the start and {end_trim} bases from the end of each read"
        );
        Ok(EdgeFilter::new(start_trim, end_trim, inverted))
    } else {
        let trim = raw.parse::<usize>().context(format!(
            "failed to parse edge filter input {raw}, should be a number"
        ))?;

        info!(
            "filtering out base modification calls {trim} bases from the \
             start and end of each read"
        );
        Ok(EdgeFilter::new(trim, trim, inverted))
    }
}

pub(crate) fn calculate_chunk_size(
    chunk_size: Option<usize>,
    interval_size: u32,
    threads: usize,
) -> usize {
    if let Some(chunk_size) = chunk_size {
        if chunk_size < threads {
            warn!(
                "chunk size {chunk_size} is less than number of threads ({}), \
                 this will limit parallelism",
                threads
            );
        }
        chunk_size
    } else {
        let cs = (threads as f32 * 1.5).floor() as usize;
        info!(
            "calculated chunk size: {cs}, interval size {}, processing {} \
             positions concurrently",
            interval_size,
            cs * interval_size as usize
        );
        cs
    }
}

pub fn parse_forward_motifs(
    input_motifs: &Option<Vec<String>>,
    cpg: bool,
) -> anyhow::Result<Option<Vec<OverlappingRegexOffset>>> {
    match (input_motifs.as_ref(), cpg) {
        (Some(raw_motifs), _) => {
            let motifs = RegexMotif::from_raw_parts(&raw_motifs, cpg)?;
            let motifs = motifs
                .into_iter()
                .map(|mot| {
                    let offset = mot.forward_offset();
                    OverlappingRegexOffset::new(mot.forward_pattern, offset)
                })
                .collect::<Vec<OverlappingRegexOffset>>();
            Ok(Some(motifs))
        }
        (None, true) => {
            let motifs = RegexMotif::from_raw_parts(&[], true).unwrap();
            let motifs = motifs
                .into_iter()
                .map(|mot| {
                    let offset = mot.forward_offset();
                    OverlappingRegexOffset::new(mot.forward_pattern, offset)
                })
                .collect::<Vec<OverlappingRegexOffset>>();
            Ok(Some(motifs))
        }
        (None, false) => Ok(None),
    }
}

pub(crate) fn get_motif_lookup_from_parts(
    regex_motifs: Option<Vec<RegexMotif>>,
    reference_fasta: Option<&PathBuf>,
    combine_strands: bool,
    mask: bool,
    preload_references: bool,
) -> anyhow::Result<Option<MotifLocationsLookup>> {
    if let Some(motifs) = regex_motifs {
        let fasta_fp = reference_fasta.as_ref().ok_or(anyhow!(
            "reference fasta is required for using --modified-bases, --motif \
             or --cpg options"
        ))?;
        if combine_strands {
            if motifs.iter().any(|rm| !rm.is_palendrome()) {
                bail!(
                    "cannot combine strands with a motif that is not a \
                     palindrome"
                )
            }
            debug!("combining + and - strand counts");
        }

        Ok(Some(MotifLocationsLookup::from_paths(
            fasta_fp,
            mask,
            None,
            motifs,
            preload_references,
        )?))
    } else {
        Ok(None)
    }
}

pub(crate) fn parse_raw_motifs(
    raw_motif_parts: &[String],
    cpg: bool,
    allow_single_base: bool,
) -> anyhow::Result<Vec<RegexMotif>> {
    if raw_motif_parts.len() % 2 != 0 {
        bail!("illegal number of parts for motif")
    }
    if raw_motif_parts.is_empty() {
        bail!("zero motif parts")
    }
    match RegexMotif::from_raw_parts(raw_motif_parts, cpg) {
        Ok(regex_motifs) => {
            for motif in regex_motifs.iter() {
                if motif.length() == 1 && !allow_single_base {
                    bail!(
                        "single base motifs not supported, use \
                         --modified-bases"
                    )
                }
            }
            Ok(regex_motifs)
        }
        e @ _ => e,
    }
}

#[cfg(test)]
mod tests {
    use super::parse_sampling_fraction;

    #[test]
    fn parse_sampling_fraction_matrix() {
        for (raw, expected) in [
            ("0", 0.0),
            ("-0.0", -0.0),
            ("0.25", 0.25),
            ("1e-3", 0.001),
            ("1", 1.0),
            ("1e0", 1.0),
        ] {
            assert_eq!(parse_sampling_fraction(raw).unwrap(), expected);
        }
        assert!(parse_sampling_fraction("-0.0").unwrap().is_sign_negative());

        for raw in [
            "-0.0000001",
            "1.0000001",
            "NaN",
            "nan",
            "inf",
            "+inf",
            "-inf",
            "1e309",
            "not-a-number",
            "",
        ] {
            assert!(
                parse_sampling_fraction(raw).is_err(),
                "expected {raw:?} to be rejected"
            );
        }
    }
}
