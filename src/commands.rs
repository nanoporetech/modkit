use std::collections::{HashMap, HashSet};
use std::io::BufWriter;
use std::num::ParseFloatError;
use std::path::PathBuf;
use std::thread;

use crate::adjust::{adjust_modbam, record_is_valid};
use anyhow::{anyhow, Context, Result as AnyhowResult};
use clap::{Args, Subcommand, ValueEnum};
use crossbeam_channel::bounded;
use histo_fp::Histogram;
use indicatif::{
    MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle,
};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::Read;

use crate::errs::{InputError, RunError};
use crate::interval_chunks::IntervalChunks;
use crate::logging::init_logging;
use crate::mod_bam::{
    format_mm_ml_tag, CollapseMethod, ModBaseInfo, RawModCode, SkipMode,
    ML_TAGS, MM_TAGS,
};
use crate::mod_base_code::{DnaBase, ModCode, ParseChar};
use crate::mod_pileup::{process_region, ModBasePileup, PileupNumericOptions};
use crate::motif_bed::{motif_bed, MotifLocations, RegexMotif};
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs;
use crate::summarize::{summarize_modbam, ModSummary};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::{calc_threshold_from_bam, Percentiles};
use crate::util;
use crate::util::{add_modkit_pg_records, get_spinner, get_targets, Region};
use crate::writers::{
    BedGraphWriter, BedMethylWriter, MultiTableWriter, OutWriter, SampledProbs,
    TableWriter, TsvWriter,
};

#[derive(Subcommand)]
pub enum Commands {
    /// Tabulates base modification calls across genomic positions. This command
    /// produces a bedMethyl formatted file. Schema and description of fields can
    /// be found in the README.
    Pileup(ModBamPileup),
    /// Performs various operations on BAM files containing base modification
    /// information, such as converting base modification codes and ignoring
    /// modification calls. Produces a BAM output file.
    AdjustMods(Adjust),
    /// Renames Mm/Ml to tags to MM/ML. Also allows changing the the mode flag from
    /// silent '.' to explicitly '?' or '.'.
    UpdateTags(Update),
    /// Calculate an estimate of the base modification probability distribution.
    SampleProbs(SampleModBaseProbs),
    /// Summarize the mod tags present in a BAM and get basic statistics. The default
    /// output is a totals table (designated by '#' lines) and a modification calls
    /// table. Descriptions of the columns can be found in the README.
    Summary(ModSummarize),
    /// Call mods from a modbam, creates a new modbam with probabilities set to 100%
    /// if a base modification is called or 0% if called canonical.
    CallMods(CallMods),
    /// Create BED file with all locations of a sequence motif.
    /// Example: modkit motif-bed CG 0
    MotifBed(MotifBed),
}

impl Commands {
    pub fn run(&self) -> Result<(), String> {
        match self {
            Self::AdjustMods(x) => x.run().map_err(|e| e.to_string()),
            Self::Pileup(x) => x.run().map_err(|e| e.to_string()),
            Self::SampleProbs(x) => x.run().map_err(|e| e.to_string()),
            Self::Summary(x) => x.run().map_err(|e| e.to_string()),
            Self::MotifBed(x) => x.run().map_err(|e| e.to_string()),
            Self::UpdateTags(x) => x.run(),
            Self::CallMods(x) => x.run().map_err(|e| e.to_string()),
        }
    }
}

type CliResult<T> = Result<T, RunError>;

fn get_sampling_options(
    no_sampling: bool,
    sampling_frac: Option<f64>,
    num_reads: usize,
) -> (Option<f64>, Option<usize>) {
    match (no_sampling, sampling_frac, num_reads) {
        // Both None tells RecordSampler to use passthrough
        // see `RecordSampler::new_from_options`
        (true, _, _) => {
            info!("not subsampling, using all reads");
            (None, None)
        }
        (false, Some(frac), _) => {
            let pct = frac * 100f64;
            info!("sampling {pct}% of reads");
            (sampling_frac, None)
        }
        (false, None, num_reads) => {
            info!("sampling {num_reads} reads from BAM");
            (None, Some(num_reads))
        }
    }
}

fn get_threshold_from_options(
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
    collapse_method: Option<&CollapseMethod>,
) -> AnyhowResult<MultipleThresholdModCaller> {
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
        collapse_method,
    )?;

    Ok(MultipleThresholdModCaller::new(
        per_base_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        0f32,
    ))
}

fn parse_raw_threshold<T: ParseChar>(raw: &str) -> AnyhowResult<(T, f32)> {
    let parts = raw.split(':').collect::<Vec<&str>>();
    if parts.len() != 2 {
        return Err(anyhow!(
            "encountered illegal per-base threshold {raw}, should \
                be <base>:<threshold>, e.g. C:0.75"
        ));
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
) -> AnyhowResult<(Option<f32>, HashMap<DnaBase, f32>)> {
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
                    return Err(anyhow!(
                        "repeated threshold for base {}",
                        dna_base.char()
                    ));
                }
            } else {
                if let Some(_) = default {
                    return Err(anyhow!(
                        "default threshold encountered more than once"
                    ));
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

fn parse_per_mod_thresholds(
    raw_per_mod_thresholds: &[String],
) -> AnyhowResult<HashMap<ModCode, f32>> {
    let per_mod_thresholds = raw_per_mod_thresholds
        .iter()
        .map(|raw| parse_raw_threshold::<ModCode>(raw))
        .collect::<anyhow::Result<HashMap<ModCode, f32>>>()?;
    per_mod_thresholds.iter().for_each(|(mod_code, thresh)| {
        info!("using threshold {thresh} for mod-code {}", mod_code.char());
    });
    Ok(per_mod_thresholds)
}

fn parse_thresholds(
    raw_base_thresholds: &[String],
    per_mod_thresholds: Option<HashMap<ModCode, f32>>,
) -> AnyhowResult<MultipleThresholdModCaller> {
    let (default, per_base_thresholds) =
        parse_per_base_thresholds(raw_base_thresholds)?;
    Ok(MultipleThresholdModCaller::new(
        per_base_thresholds,
        per_mod_thresholds.unwrap_or(HashMap::new()),
        default.unwrap_or(0f32),
    ))
}

#[derive(Args)]
pub struct Adjust {
    /// BAM file to collapse mod call from.
    in_bam: PathBuf,
    /// File path to new BAM file to be created.
    out_bam: PathBuf,
    /// Output debug logs to file at this path.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Modified base code to ignore/remove, see
    /// https://samtools.github.io/hts-specs/SAMtags.pdf for details on
    /// the modified base codes.
    #[arg(long, conflicts_with = "convert", default_value_t = 'h')]
    ignore: char,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Fast fail, stop processing at the first invalid sequence record. Default
    /// behavior is to continue and report failed/skipped records at the end.
    #[arg(short, long = "ff", default_value_t = false)]
    fail_fast: bool,
    /// Convert one mod-tag to another, summing the probabilities together if
    /// the retained mod tag is already present.
    #[arg(group = "prob_args", long, action = clap::ArgAction::Append, num_args = 2)]
    convert: Option<Vec<char>>,
}

impl Adjust {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let fp = &self.in_bam;
        let out_fp = &self.out_bam;
        let mut reader = bam::Reader::from_path(fp)?;
        let threads = self.threads;
        reader.set_threads(threads)?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);
        let mut out_bam =
            bam::Writer::from_path(out_fp, &header, bam::Format::Bam)?;

        let methods = if let Some(convert) = &self.convert {
            let mut conversions = HashMap::new();
            for chunk in convert.chunks(2) {
                debug_assert_eq!(chunk.len(), 2);
                let from: RawModCode = chunk[0];
                let to: RawModCode = chunk[1];
                let froms = conversions.entry(to).or_insert(HashSet::new());
                froms.insert(from);
            }
            for (to_code, from_codes) in conversions.iter() {
                info!(
                    "Converting {} to {}",
                    from_codes.iter().collect::<String>(),
                    to_code
                )
            }
            conversions
                .into_iter()
                .map(|(to_mod_code, from_mod_codes)| {
                    let method = CollapseMethod::Convert {
                        to: to_mod_code,
                        from: from_mod_codes,
                    };

                    method
                })
                .collect::<Vec<CollapseMethod>>()
        } else {
            info!(
                "Removing mod base {} from {}, new bam {}",
                self.ignore,
                fp.to_str().unwrap_or("???"),
                out_fp.to_str().unwrap_or("???")
            );
            let method = CollapseMethod::ReDistribute(self.ignore);
            vec![method]
        };

        adjust_modbam(
            &mut reader,
            &mut out_bam,
            &methods,
            None,
            self.fail_fast,
            "Adjusting modBAM",
        )?;
        Ok(())
    }
}

#[derive(Args)]
pub struct ModBamPileup {
    // running args
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,
    /// Output file (or directory with --bedgraph option) to write results into.
    /// Specify "-" or "stdout" to direct output to stdout.
    out_bed: String,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Process only the specified region of the BAM when performing pileup.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    region: Option<String>,

    // processing args
    /// Number of threads to use while processing chunks concurrently.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Interval chunk size to process concurrently. Smaller interval chunk
    /// sizes will use less memory but incur more overhead.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,

    // sampling args
    /// Sample this many reads when estimating the filtering threshold. Reads will
    /// be sampled evenly across aligned genome. If a region is specified, either with
    /// the --region option or the --sample-region option, then reads will be sampled
    /// evenly across the region given. This option is useful for large BAM files.
    /// In practice, 10-50 thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold.
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Sample this fraction of the reads when estimating the filter-percentile.
    /// In practice, 50-100 thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold. See filtering.md for
    /// details on filtering.
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Set a random seed for deterministic running, the default is non-deterministic.
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Do not perform any filtering, include all mod base calls in output. See
    /// filtering.md for details on filtering.
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will filter
    /// out the 10% lowest confidence modification calls.
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per-base. Global filter threshold
    /// can be specified with by a decimal number (e.g. 0.75). Per-base thresholds
    /// can be specified by colon-separated values, for example C:0.75 specifies a
    /// threshold value of 0.75 for cytosine modification calls. Additional
    /// per-base thresholds can be specified by repeating the option: for example
    /// --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a single
    /// base option and a default for all other bases with:
    /// --filter-threshold A:0.70 --filter-threshold 0.9 will specify a threshold
    /// value of 0.70 for adenosine and 0.9 for all other base modification calls.
    #[arg(
    long,
    group = "thresholds",
    action = clap::ArgAction::Append,
    alias = "pass_threshold"
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent of the
    /// threshold for the primary sequence base or the default. For example, to set
    /// the pass threshold for 5hmC to 0.8 use `--mod-threshold h:0.8`. The pass
    /// threshold will still be estimated as usual and used for canonical cytosine and
    /// 5mC unless the `--filter-threshold` option is also passed. See the online
    /// documentation for more details.
    #[arg(
    long,
    alias = "mod-threshold",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Specify a region for sampling reads from when estimating the threshold probability.
    /// If this option is not provided, but --region is provided, the genomic interval
    /// passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    sample_region: Option<String>,
    /// Interval chunk size to process concurrently when estimating the threshold
    /// probability, can be larger than the pileup processing interval.
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,

    // collapsing and combining args
    /// Ignore a modified base class  _in_situ_ by redistributing base modification
    /// probability equally across other options. For example, if collapsing 'h',
    /// with 'm' and canonical options, half of the probability of 'h' will be added to
    /// both 'm' and 'C'. A full description of the methods can be found in
    /// collapse.md.
    #[arg(long, group = "combine_args", hide_short_help = true)]
    ignore: Option<char>,
    /// Force allow implicit-canonical mode. By default modkit does not allow
    /// pileup with the implicit mode ('.', or silent). The `update-tags`
    /// subcommand is provided to update tags to the new mode. This option allows
    /// the interpretation of implicit mode tags: residues without modified
    /// base probability will be interpreted as being the non-modified base.
    /// We do not recommend using this option.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = false,
        hide_short_help = true
    )]
    force_allow_implicit: bool,
    /// Only output counts at CpG motifs. Requires a reference sequence to be
    /// provided.
    #[arg(long, requires = "reference_fasta", default_value_t = false)]
    cpg: bool,
    /// Reference sequence in FASTA format. Required for CpG motif filtering.
    #[arg(long = "ref")]
    reference_fasta: Option<PathBuf>,
    /// Respect soft masking in the reference FASTA.
    #[arg(
        long,
        short = 'k',
        requires = "reference_fasta",
        default_value_t = false,
        hide_short_help = true
    )]
    mask: bool,
    /// Optional preset options for specific applications.
    /// traditional: Prepares bedMethyl analogous to that generated from other technologies
    /// for the analysis of 5mC modified bases. Shorthand for --cpg --combine-strands
    /// --ignore h.
    #[arg(
        long,
        requires = "reference_fasta",
        conflicts_with_all = ["combine_mods", "cpg", "combine_strands", "ignore"],
    )]
    preset: Option<Presets>,
    /// Combine base modification calls, all counts of modified bases are summed together. See
    /// collapse.md for details.
    #[arg(
        long,
        default_value_t = false,
        group = "combine_args",
        hide_short_help = true
    )]
    combine_mods: bool,
    /// When performing CpG analysis, sum the counts from the positive and
    /// negative strands into the counts for the positive strand.
    #[arg(long, requires = "cpg", default_value_t = false)]
    combine_strands: bool,

    // output args
    /// For bedMethyl output, separate columns with only tabs. The default is
    /// to use tabs for the first 10 fields and spaces thereafter. The
    /// default behavior is more likely to be compatible with genome viewers.
    /// Enabling this option may make it easier to parse the output with
    /// tabular data handlers that expect a single kind of separator.
    #[arg(
        long,
        conflicts_with = "bedgraph",
        default_value_t = false,
        hide_short_help = true
    )]
    only_tabs: bool,
    /// Output bedGraph format, see https://genome.ucsc.edu/goldenPath/help/bedgraph.html.
    /// For this setting, specify a directory for output files to be make in.
    /// Two files for each modification will be produced, one for the positive strand
    /// and one for the negative strand. So for 5mC (m) and 5hmC (h) there will be 4 files
    /// produced.
    #[arg(
        long,
        conflicts_with = "only_tabs",
        default_value_t = false,
        hide_short_help = true
    )]
    bedgraph: bool,
    /// Prefix to prepend on bedgraph output file names. Without this option the files
    /// will be <mod_code>_<strand>.bedgraph
    #[arg(long, requires = "bedgraph")]
    prefix: Option<String>,
}

impl ModBamPileup {
    fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        // do this first so we fail when the file isn't readable
        let header = bam::IndexedReader::from_path(&self.in_bam)
            .map(|reader| reader.header().to_owned())?;
        let region = if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Some(Region::parse_str(raw_region, &header)?)
        } else {
            None
        };
        let sampling_region = if let Some(raw_region) = &self.sample_region {
            info!("parsing sample region {raw_region}");
            Some(Region::parse_str(raw_region, &header)?)
        } else {
            None
        };

        let (pileup_options, combine_strands, threshold_collapse_method) =
            match self.preset {
                Some(Presets::traditional) => {
                    info!("ignoring mod code {}", ModCode::h.char());
                    info!(
                        "NOTICE, in the next version of modkit the 'traditional' preset \
                         will perform --combine-mods instead of --ignore h"
                    );
                    (
                        PileupNumericOptions::Collapse(
                            CollapseMethod::ReDistribute('h'),
                        ),
                        true,
                        Some(CollapseMethod::ReDistribute('h')),
                    )
                }
                None => {
                    let (options, collapse_method) =
                        match (self.combine_mods, &self.ignore) {
                            (false, None) => {
                                (PileupNumericOptions::Passthrough, None)
                            }
                            (true, _) => (PileupNumericOptions::Combine, None),
                            (_, Some(raw_mod_code)) => {
                                info!("ignoring mod code {}", raw_mod_code);
                                let method =
                                    CollapseMethod::ReDistribute(*raw_mod_code);
                                (
                                    PileupNumericOptions::Collapse(
                                        method.clone(),
                                    ),
                                    Some(method),
                                )
                            }
                        };
                    (options, self.combine_strands, collapse_method)
                }
            };

        // setup the writer here so we fail before doing any work (if there are problems).
        let out_fp_str = self.out_bed.clone();
        let mut writer: Box<dyn OutWriter<ModBasePileup>> = if self.bedgraph {
            Box::new(BedGraphWriter::new(&out_fp_str, self.prefix.as_ref())?)
        } else {
            match out_fp_str.as_str() {
                "stdout" | "-" => {
                    let writer = BufWriter::new(std::io::stdout());
                    Box::new(BedMethylWriter::new(writer, self.only_tabs))
                }
                _ => {
                    let fh = std::fs::File::create(out_fp_str)
                        .context("failed to make output file")?;
                    let writer = BufWriter::new(fh);
                    Box::new(BedMethylWriter::new(writer, self.only_tabs))
                }
            }
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .with_context(|| "failed to make threadpool")?;

        let per_mod_thresholds =
            if let Some(raw_per_mod_thresholds) = &self.mod_thresholds {
                Some(parse_per_mod_thresholds(raw_per_mod_thresholds)?)
            } else {
                None
            };
        let threshold_caller =
            if let Some(raw_threshold) = &self.filter_threshold {
                parse_thresholds(raw_threshold, per_mod_thresholds)?
            } else {
                pool.install(|| {
                    get_threshold_from_options(
                        &self.in_bam,
                        self.threads,
                        self.sampling_interval_size,
                        self.sampling_frac,
                        self.num_reads,
                        self.no_filtering,
                        self.filter_percentile,
                        self.seed,
                        sampling_region.as_ref().or(region.as_ref()),
                        per_mod_thresholds,
                        threshold_collapse_method.as_ref(),
                    )
                })?
            };

        if !self.no_filtering {
            for (base, threshold) in threshold_caller.iter_thresholds() {
                let base = base.char();
                match (threshold * 100f32).ceil() as usize {
                    0..=60 => error!(
                "Threshold of {threshold} for base {base} is very low. Consider increasing the \
                filter-percentile or specifying a higher threshold."),
                    61..=70 => warn!(
                "Threshold of {threshold} for base {base} is low. Consider increasing the \
                filter-percentile or specifying a higher threshold."
            ),
                    _ => info!("Using filter threshold {} for {base}.", threshold),
                }
            }
            for (base, threshold) in threshold_caller.iter_mod_thresholds() {
                let base = base.char();
                match (threshold * 100f32).ceil() as usize {
                    0..=60 => error!(
                "Threshold of {threshold} for mod code {base} is very low. Consider increasing the \
                filter-percentile or specifying a higher threshold."),
                    61..=70 => warn!(
                "Threshold of {threshold} for mod code {base} is low. Consider increasing the \
                filter-percentile or specifying a higher threshold."
            ),
                    _ => info!("Using filter threshold {} for mod code {base}.", threshold),
                }
            }
        }

        let tids = get_targets(&header, region.as_ref());
        let use_cpg_motifs = self.cpg
            || self
                .preset
                .map(|preset| match preset {
                    Presets::traditional => true,
                })
                .unwrap_or(false);
        let (motif_locations, tids) = if use_cpg_motifs {
            let fasta_fp = self
                .reference_fasta
                .as_ref()
                .ok_or(anyhow!("reference fasta is required for CpG"))?;
            let regex_motif = RegexMotif::parse_string("CG", 0).unwrap();
            debug!("filtering output to only CpG motifs");
            if combine_strands {
                debug!("combining + and - strand counts");
            }
            let names_to_tid = tids
                .iter()
                .map(|target| (target.name.as_str(), target.tid))
                .collect::<HashMap<&str, u32>>();
            let motif_locations = pool.install(|| {
                MotifLocations::from_fasta(
                    fasta_fp,
                    regex_motif,
                    &names_to_tid,
                    self.mask,
                )
            })?;
            let filtered_tids = motif_locations.filter_reference_records(tids);
            (Some(motif_locations), filtered_tids)
        } else {
            (None, tids)
        };

        let (snd, rx) = bounded(1_000); // todo figure out sane default for this?
        let in_bam_fp = self.in_bam.clone();
        let interval_size = self.interval_size;

        let master_progress = MultiProgress::new();
        let sty = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.green/yellow} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");
        let tid_progress = master_progress
            .add(ProgressBar::new(tids.len() as u64))
            .with_style(sty.clone());
        tid_progress.set_message("contigs");
        let write_progress = master_progress.add(get_spinner());
        write_progress.set_message("rows written");

        let force_allow = self.force_allow_implicit;

        let interval_style = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");

        thread::spawn(move || {
            pool.install(|| {
                for target in tids {
                    let intervals = IntervalChunks::new(
                        target.start,
                        target.length,
                        interval_size,
                        target.tid,
                        motif_locations.as_ref(),
                    )
                    .collect::<Vec<(u32, u32)>>();
                    let n_intervals = intervals.len();
                    let interval_progress = master_progress.add(
                        ProgressBar::new(n_intervals as u64)
                            .with_style(interval_style.clone()),
                    );
                    interval_progress
                        .set_message(format!("processing {}", &target.name));
                    let mut result: Vec<Result<ModBasePileup, String>> = vec![];
                    let (res, _) = rayon::join(
                        || {
                            intervals
                                .into_par_iter()
                                .progress_with(interval_progress)
                                .map(|(start, end)| {
                                    process_region(
                                        &in_bam_fp,
                                        target.tid,
                                        start,
                                        end,
                                        &threshold_caller,
                                        &pileup_options,
                                        force_allow,
                                        combine_strands,
                                        motif_locations.as_ref(),
                                    )
                                })
                                .collect::<Vec<Result<ModBasePileup, String>>>()
                        },
                        || {
                            result.into_iter().for_each(|mod_base_pileup| {
                                snd.send(mod_base_pileup)
                                    .expect("failed to send")
                            });
                        },
                    );
                    result = res;
                    result.into_iter().for_each(|pileup| {
                        snd.send(pileup).expect("failed to send")
                    });
                    tid_progress.inc(1);
                }
                tid_progress.finish_and_clear();
            });
        });

        for result in rx.into_iter() {
            match result {
                Ok(mod_base_pileup) => {
                    let rows_written = writer.write(mod_base_pileup)?;
                    write_progress.inc(rows_written);
                }
                Err(message) => {
                    debug!("> unexpected error {message}");
                }
            }
        }
        let rows_processed = write_progress.position();
        write_progress.finish_and_clear();
        info!("Done, processed {rows_processed} rows.");
        Ok(())
    }
}

fn parse_percentiles(
    raw_percentiles: &str,
) -> Result<Vec<f32>, ParseFloatError> {
    if raw_percentiles.contains("..") {
        todo!("handle parsing ranges")
    } else {
        raw_percentiles
            .split(',')
            .map(|x| x.parse::<f32>())
            .collect()
    }
}

#[derive(Args)]
pub struct SampleModBaseProbs {
    /// Input BAM with modified base tags. If a index is found
    /// reads will be sampled evenly across the length of the
    /// reference sequence.
    in_bam: PathBuf,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Percentiles to calculate, a space separated list of floats.
    #[arg(short, long, default_value_t=String::from("0.1,0.5,0.9"))]
    percentiles: String,
    /// Directory to deposit result tables into. Required for model probability
    /// histogram output. Creates two files probabilities.tsv and probabilities.txt
    /// The .txt contains ASCII-histograms and the .tsv contains tab-separated variable
    /// data represented by the histograms.
    #[arg(short = 'o', long)]
    out_dir: Option<PathBuf>,
    /// Label to prefix output files with. E.g. 'foo' will output
    /// foo_thresholds.tsv, foo_probabilities.tsv, and foo_probabilities.txt.
    #[arg(long, requires = "out_dir")]
    prefix: Option<String>,
    /// Overwrite results if present.
    #[arg(long, requires = "out_dir", default_value_t = false)]
    force: bool,
    /// Ignore a modified base class  _in_situ_ by redistributing base modification
    /// probability equally across other options. For example, if collapsing 'h',
    /// with 'm' and canonical options, half of the probability of 'h' will be added to
    /// both 'm' and 'C'. A full description of the methods can be found in
    /// collapse.md.
    #[arg(long, hide_short_help = true)]
    ignore: Option<char>,

    // probability histogram options
    /// Output histogram of base modification prediction probabilities.
    #[arg(long = "hist", requires = "out_dir", default_value_t = false)]
    histogram: bool,
    /// Number of buckets for the histogram, if used.
    #[arg(long, requires = "histogram", default_value_t = 128)]
    buckets: u64,

    /// Max number of reads to use, especially recommended when using a large
    /// BAM without an index. If an indexed BAM is provided, the reads will be
    /// sampled evenly over the length of the aligned reference. If a region is
    /// passed with the --region option, they will be sampled over the genomic
    /// region.
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Instead of using a defined number of reads, specify a fraction of reads
    /// to sample, for example 0.1 will sample 1/10th of the reads.
    #[arg(group = "sampling_options", short = 'f', long)]
    sampling_frac: Option<f64>,
    /// No sampling, use all of the reads to calculate the filter thresholds.
    #[arg(long, group = "sampling_options", default_value_t = false)]
    no_sampling: bool,
    /// Random seed for deterministic running, the default is non-deterministic.
    #[arg(short, requires = "sampling_frac", long)]
    seed: Option<u64>,

    /// Process only the specified region of the BAM when collecting probabilities.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// Interval chunk size to process concurrently. Smaller interval chunk
    /// sizes will use less memory but incur more overhead. Only used when
    /// sampling probs from an indexed bam.
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
}

impl SampleModBaseProbs {
    fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let reader = bam::Reader::from_path(&self.in_bam)?;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let region = if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Some(Region::parse_str(raw_region, reader.header())?)
        } else {
            None
        };

        let (sample_frac, num_reads) = get_sampling_options(
            self.no_sampling,
            self.sampling_frac,
            self.num_reads,
        );

        let collapse_method = if let Some(raw_mod_code_to_ignore) = self.ignore
        {
            let _ = ModCode::parse_raw_mod_code(raw_mod_code_to_ignore)?;
            Some(CollapseMethod::ReDistribute(raw_mod_code_to_ignore))
        } else {
            None
        };

        let desired_percentiles = parse_percentiles(&self.percentiles)
            .with_context(|| {
                format!("failed to parse percentiles: {}", &self.percentiles)
            })?;

        pool.install(|| {
            let read_ids_to_base_mod_calls =
                get_sampled_read_ids_to_base_mod_probs(
                    &self.in_bam,
                    self.threads,
                    self.interval_size,
                    sample_frac,
                    num_reads,
                    self.seed,
                    region.as_ref(),
                    collapse_method.as_ref(),
                )?;

            let histograms = if self.histogram {
                let mod_call_probs =
                    read_ids_to_base_mod_calls.mle_probs_per_base_mod();
                Some(
                    mod_call_probs
                        .iter()
                        .map(|(base, calls)| {
                            let mut hist =
                                Histogram::with_buckets(self.buckets, Some(0));
                            for prob in calls {
                                hist.add(*prob)
                            }
                            (*base, hist)
                        })
                        .collect::<HashMap<char, Histogram>>(),
                )
            } else {
                None
            };

            let percentiles = read_ids_to_base_mod_calls
                .mle_probs_per_base()
                .into_iter()
                .map(|(canonical_base, mut probs)| {
                    Percentiles::new(&mut probs, &desired_percentiles)
                        .with_context(|| {
                            format!(
                                "failed to calculate threshold for base {}",
                                canonical_base.char()
                            )
                        })
                        .map(|percs| (canonical_base.char(), percs))
                })
                .collect::<AnyhowResult<HashMap<char, Percentiles>>>()?;

            let sampled_probs =
                SampledProbs::new(histograms, percentiles, self.prefix.clone());

            let mut writer: Box<dyn OutWriter<SampledProbs>> =
                if let Some(p) = &self.out_dir {
                    sampled_probs.check_path(p, self.force)?;
                    Box::new(MultiTableWriter::new(p.clone()))
                } else {
                    Box::new(TsvWriter::new_stdout())
                };

            writer.write(sampled_probs)?;

            Ok(())
        })
    }
}

#[derive(Args)]
pub struct ModSummarize {
    /// Input ModBam file.
    in_bam: PathBuf,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Output summary as a tab-separated variables stdout instead of a table.
    #[arg(long = "tsv", default_value_t = false)]
    tsv_format: bool,

    // sampling options
    /// Max number of reads to use for estimating the filter threshold and
    /// generating the summary, especially recommended when using a large
    /// BAM without an index. If an indexed BAM is provided, the reads will
    /// be sampled evenly over the length of the aligned reference. If a
    /// region is passed with the --region option, they will be sampled
    /// over the genomic region.
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Instead of using a defined number of reads, specify a fraction of reads
    /// to sample when estimating the filter threshold. For example 0.1 will
    /// sample 1/10th of the reads.
    #[arg(group = "sampling_options", short = 'f', long)]
    sampling_frac: Option<f64>,
    /// No sampling, use all of the reads to calculate the filter thresholds and
    /// generating the summary.
    #[arg(long, group = "sampling_options", default_value_t = false)]
    no_sampling: bool,
    /// Sets a random seed for deterministic running (when using --sample-frac),
    /// the default is non-deterministic.
    #[arg(short, requires = "sampling_frac", long)]
    seed: Option<u64>,

    // threshold options
    /// Do not perform any filtering, include all base modification calls in the
    /// summary. See filtering.md for details on filtering.
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will filter
    /// out the 10% lowest confidence base modification calls.
    #[arg(group = "thresholds", short = 'p', long, default_value_t = 0.1)]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per-base. Global filter threshold
    /// can be specified with by a decimal number (e.g. 0.75). Per-base thresholds
    /// can be specified by colon-separated values, for example C:0.75 specifies a
    /// threshold value of 0.75 for cytosine modification calls. Additional
    /// per-base thresholds can be specified by repeating the option: for example
    /// --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a single
    /// base option and a default for all other bases with:
    /// --filter-threshold A:0.70 --filter-threshold 0.9 will specify a threshold
    /// value of 0.70 for adenosine and 0.9 for all other base modification calls.
    #[arg(
        long,
        group = "thresholds",
        action = clap::ArgAction::Append
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent of the
    /// threshold for the primary sequence base or the default. For example, to set
    /// the pass threshold for 5hmC to 0.8 use `--mod-threshold h:0.8`. The pass
    /// threshold will still be estimated as usual and used for canonical cytosine and
    /// 5mC unless the `--filter-threshold` option is also passed. See the online
    /// documentation for more details.
    #[arg(
    long,
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Ignore a modified base class  _in_situ_ by redistributing base modification
    /// probability equally across other options. For example, if collapsing 'h',
    /// with 'm' and canonical options, half of the probability of 'h' will be added to
    /// both 'm' and 'C'. A full description of the methods can be found in
    /// collapse.md.
    #[arg(long, group = "combine_args", hide_short_help = true)]
    ignore: Option<char>,

    /// Process only the specified region of the BAM when collecting probabilities.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// When using regions, interval chunk size to process concurrently.
    /// Smaller interval chunk sizes will use less memory but incur more
    /// overhead.
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
}

impl ModSummarize {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let reader = bam::Reader::from_path(&self.in_bam)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        let region = if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Some(Region::parse_str(&raw_region, reader.header())?)
        } else {
            None
        };

        let (sample_frac, num_reads) = get_sampling_options(
            self.no_sampling,
            self.sampling_frac,
            self.num_reads,
        );

        let per_mod_thresholds =
            if let Some(raw_per_mod_thresholds) = &self.mod_thresholds {
                Some(parse_per_mod_thresholds(raw_per_mod_thresholds)?)
            } else {
                None
            };

        let filter_thresholds =
            if let Some(raw_thresholds) = &self.filter_threshold {
                info!("parsing user defined thresholds");
                Some(parse_thresholds(
                    raw_thresholds,
                    per_mod_thresholds.clone(),
                )?)
            } else if self.no_filtering {
                info!("not performing filtering");
                Some(MultipleThresholdModCaller::new_passthrough())
            } else {
                None
            };

        let collapse_method = if let Some(raw_mod_code_to_ignore) = self.ignore
        {
            let _ = ModCode::parse_raw_mod_code(raw_mod_code_to_ignore)?;
            Some(CollapseMethod::ReDistribute(raw_mod_code_to_ignore))
        } else {
            None
        };

        let mod_summary = pool.install(|| {
            summarize_modbam(
                &self.in_bam,
                self.threads,
                self.interval_size,
                sample_frac,
                num_reads,
                self.seed,
                region.as_ref(),
                self.filter_percentile,
                filter_thresholds,
                per_mod_thresholds,
                collapse_method.as_ref(),
            )
        })?;

        let mut writer: Box<dyn OutWriter<ModSummary>> = if self.tsv_format {
            Box::new(TsvWriter::new_stdout())
        } else {
            Box::new(TableWriter::new())
        };
        writer.write(mod_summary)?;
        Ok(())
    }
}

#[derive(Args)]
pub struct MotifBed {
    /// Input FASTA file
    fasta: PathBuf,
    /// Motif to search for within FASTA, e.g. CG
    motif: String,
    /// Offset within motif, e.g. 0
    offset: usize,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
}

impl MotifBed {
    fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(None);
        motif_bed(&self.fasta, &self.motif, self.offset, self.mask)
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[allow(non_camel_case_types)]
enum Presets {
    traditional,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[allow(non_camel_case_types)]
enum ModMode {
    ambiguous,
    implicit,
}

impl ModMode {
    fn to_skip_mode(self) -> SkipMode {
        match self {
            Self::ambiguous => SkipMode::Ambiguous,
            Self::implicit => SkipMode::ProbModified,
        }
    }
}

#[derive(Args)]
pub struct Update {
    /// BAM file to update modified base tags in.
    in_bam: PathBuf,
    /// File path to new BAM file to be created.
    out_bam: PathBuf,
    /// Mode, change mode to this value, options {'ambiguous', 'implicit'}.
    /// See spec at: https://samtools.github.io/hts-specs/SAMtags.pdf.
    /// 'ambiguous' ('?') means residues without explicit modification
    /// probabilities will not be assumed canonical or modified. 'implicit'
    /// means residues without explicit modification probabilities are
    /// assumed to be canonical.
    #[arg(short, long, value_enum)]
    mode: Option<ModMode>,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Output debug logs to file at this path.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
}

fn update_mod_tags(
    mut record: bam::Record,
    new_mode: Option<SkipMode>,
) -> CliResult<bam::Record> {
    let _ok = record_is_valid(&record)?;
    let mod_base_info = ModBaseInfo::new_from_record(&record)?;
    let mm_style = mod_base_info.mm_style;
    let ml_style = mod_base_info.ml_style;

    let mut mm_agg = String::new();
    let mut ml_agg = Vec::new();

    let (converters, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();
    for (base, strand, mut seq_pos_mod_probs) in mod_prob_iter {
        let converter = converters.get(&base).unwrap();
        if let Some(mode) = new_mode {
            seq_pos_mod_probs.skip_mode = mode;
        }
        let (mm, mut ml) =
            format_mm_ml_tag(seq_pos_mod_probs, strand, converter);
        mm_agg.push_str(&mm);
        ml_agg.extend_from_slice(&mut ml);
    }
    record
        .remove_aux(mm_style.as_bytes())
        .expect("failed to remove MM tag");
    record
        .remove_aux(ml_style.as_bytes())
        .expect("failed to remove ML tag");
    let mm = Aux::String(&mm_agg);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml_agg;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record
        .push_aux(MM_TAGS[0].as_bytes(), mm)
        .expect("failed to add MM tag");
    record
        .push_aux(ML_TAGS[0].as_bytes(), ml)
        .expect("failed to add ML tag");

    Ok(record)
}

impl Update {
    fn run(&self) -> Result<(), String> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let fp = &self.in_bam;
        let out_fp = &self.out_bam;
        let threads = self.threads;
        let mut reader =
            bam::Reader::from_path(fp).map_err(|e| e.to_string())?;
        reader.set_threads(threads).map_err(|e| e.to_string())?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);

        let mut out_bam =
            bam::Writer::from_path(out_fp, &header, bam::Format::Bam)
                .map_err(|e| e.to_string())?;
        let spinner = get_spinner();

        spinner.set_message("Updating ModBAM");
        let mut total = 0usize;
        let mut total_failed = 0usize;
        let mut total_skipped = 0usize;

        for (i, result) in reader.records().enumerate() {
            if let Ok(record) = result {
                let record_name = util::get_query_name_string(&record)
                    .unwrap_or("???".to_owned());
                match update_mod_tags(
                    record,
                    self.mode.map(|m| m.to_skip_mode()),
                ) {
                    Err(RunError::BadInput(InputError(err)))
                    | Err(RunError::Failed(err)) => {
                        debug!("read {} failed, {}", record_name, err);
                        total_failed += 1;
                    }
                    Err(RunError::Skipped(_reason)) => {
                        total_skipped += 1;
                    }
                    Ok(record) => {
                        if let Err(err) = out_bam.write(&record) {
                            debug!("failed to write {}", err);
                            total_failed += 1;
                        } else {
                            spinner.inc(1);
                            total = i;
                        }
                    }
                }
            } else {
                total_failed += 1;
            }
        }

        spinner.finish_and_clear();

        info!(
            "done, {} records processed, {} failed, {} skipped",
            total, total_failed, total_skipped
        );
        Ok(())
    }
}

#[derive(Args)]
pub struct CallMods {
    // running args
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,
    /// Output BAM filepath.
    out_bam: PathBuf,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    // /// Process only the specified region of the BAM when performing transformation.
    // /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    // #[arg(long)] todo(arand)
    // region: Option<String>,
    /// Fast fail, stop processing at the first invalid sequence record. Default
    /// behavior is to continue and report failed/skipped records at the end.
    #[arg(long = "ff", default_value_t = false)]
    fail_fast: bool,

    // processing args
    /// Number of threads to use while processing chunks concurrently.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    // /// Interval chunk size to process concurrently. Smaller interval chunk
    // /// sizes will use less memory but incur more overhead. Only used when
    // /// provided an indexed BAM.
    // #[arg( todo(arand)
    // short = 'i',
    // long,
    // default_value_t = 100_000,
    // hide_short_help = true
    // )]
    // interval_size: u32,

    // sampling args
    /// Sample this many reads when estimating the filtering threshold. If alignments are
    /// present reads will be sampled evenly across aligned genome. If a region is
    /// specified, either with the --region option or the --sample-region option, then
    /// reads will be sampled evenly across the region given. This option is useful for
    /// large BAM files. In practice, 10-50 thousand reads is sufficient to estimate the
    /// model output distribution and determine the filtering threshold.
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Sample this fraction of the reads when estimating the filter-percentile.
    /// In practice, 50-100 thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold. See filtering.md for
    /// details on filtering.
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Set a random seed for deterministic running, the default is non-deterministic.
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Specify a region for sampling reads from when estimating the threshold probability.
    /// If this option is not provided, but --region is provided, the genomic interval
    /// passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    sample_region: Option<String>,
    /// Interval chunk size to process concurrently when estimating the threshold
    /// probability, can be larger than the pileup processing interval.
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,

    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will filter
    /// out the 10% lowest confidence modification calls.
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per primary base. A global filter
    /// threshold can be specified with by a decimal number (e.g. 0.75). Per-base
    /// thresholds can be specified by colon-separated values, for example C:0.75
    /// specifies a threshold value of 0.75 for cytosine modification calls. Additional
    /// per-base thresholds can be specified by repeating the option: for example
    /// --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a single
    /// base option and a default for all other bases with:
    /// --filter-threshold A:0.70 --filter-threshold 0.9 will specify a threshold
    /// value of 0.70 for adenosine and 0.9 for all other base modification calls.
    #[arg(
    long,
    group = "thresholds",
    action = clap::ArgAction::Append,
    alias = "pass_threshold"
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent of the
    /// threshold for the primary sequence base or the default. For example, to set
    /// the pass threshold for 5hmC to 0.8 use `--mod-threshold h:0.8`. The pass
    /// threshold will still be estimated as usual and used for canonical cytosine and
    /// 5mC unless the `--filter-threshold` option is also passed. See the online
    /// documentation for more details.
    #[arg(
    long = "mod-threshold",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Don't filter base modification calls, assign each base modification to the
    /// highest probability prediction.
    #[arg(long, default_value_t = false)]
    no_filtering: bool,
}

impl CallMods {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

        let mut reader = bam::Reader::from_path(&self.in_bam)?;
        let threads = self.threads;
        reader.set_threads(threads)?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);
        let mut out_bam =
            bam::Writer::from_path(&self.out_bam, &header, bam::Format::Bam)?;

        let per_mod_thresholds =
            if let Some(raw_per_mod_thresholds) = &self.mod_thresholds {
                Some(parse_per_mod_thresholds(raw_per_mod_thresholds)?)
            } else {
                None
            };

        let sampling_region = if let Some(raw_region) = &self.sample_region {
            info!("parsing sample region {raw_region}");
            Some(Region::parse_str(raw_region, &reader.header())?)
        } else {
            None
        };

        let caller = if let Some(raw_threshold) = &self.filter_threshold {
            parse_thresholds(raw_threshold, per_mod_thresholds)?
        } else {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .with_context(|| "failed to make threadpool")?;
            pool.install(|| {
                get_threshold_from_options(
                    &self.in_bam,
                    self.threads,
                    self.sampling_interval_size,
                    self.sampling_frac,
                    self.num_reads,
                    false,
                    self.filter_percentile,
                    self.seed,
                    sampling_region.as_ref(),
                    per_mod_thresholds,
                    None,
                )
            })?
        };

        adjust_modbam(
            &mut reader,
            &mut out_bam,
            &[],
            Some(&caller),
            self.fail_fast,
            "Calling Mods",
        )?;

        Ok(())
    }
}
