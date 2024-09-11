use std::collections::{HashMap, HashSet};
use std::num::ParseFloatError;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result as AnyhowResult};
use clap::{Args, Subcommand, ValueEnum};
use itertools::Itertools;
use log::{debug, info};
use rust_htslib::bam::{
    self,
    record::{Aux, AuxArray},
    Read,
};
use rust_htslib::tpool;

use crate::adjust::adjust_modbam;
use crate::command_utils::{
    get_bam_writer, get_serial_reader, get_threshold_from_options,
    parse_edge_filter_input, parse_per_mod_thresholds, parse_thresholds,
    using_stream,
};
use crate::dmr::subcommands::BedMethylDmr;
use crate::entropy::subcommand::MethylationEntropy;
use crate::errs::{InputError, RunError};
use crate::extract::subcommand::ExtractMods;
use crate::find_motifs::subcommand::EntryMotifs;
use crate::localise::subcommand::EntryLocalize;
use crate::logging::init_logging;
use crate::mod_bam::{
    format_mm_ml_tag, CollapseMethod, ModBaseInfo, SkipMode, ML_TAGS, MM_TAGS,
};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::monoid::Moniod;
use crate::pileup::subcommand::{DuplexModBamPileup, ModBamPileup};
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::get_sampled_read_ids_to_base_mod_probs;
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::record_processor::RecordProcessor;
use crate::repair_tags::RepairTags;
use crate::stats::subcommand::EntryStats;
use crate::summarize::{sampled_reads_to_summary, ModSummary};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::{calc_thresholds_per_base, Percentiles};
use crate::util;
use crate::util::{add_modkit_pg_records, get_targets, get_ticker, Region};
use crate::validate::subcommand::ValidateFromModbam;
use crate::writers::{
    MultiTableWriter, OutWriter, SampledProbs, TableWriter, TsvWriter,
};

#[derive(Subcommand)]
pub enum Commands {
    /// Tabulates base modification calls across genomic positions. This
    /// command produces a bedMethyl formatted file. Schema and description
    /// of fields can be found in the README.
    Pileup(ModBamPileup),
    /// Performs various operations on BAM files containing base modification
    /// information, such as converting base modification codes and ignoring
    /// modification calls. Produces a BAM output file.
    AdjustMods(Adjust),
    /// Renames Mm/Ml to tags to MM/ML. Also allows changing the mode flag
    /// from silent '.' to explicitly '?' or '.'.
    UpdateTags(Update),
    /// Calculate an estimate of the base modification probability
    /// distribution.
    SampleProbs(SampleModBaseProbs),
    /// Summarize the mod tags present in a BAM and get basic statistics. The
    /// default output is a totals table (designated by '#' lines) and a
    /// modification calls table. Descriptions of the columns can be found
    /// in the README.
    Summary(ModSummarize),
    /// Call mods from a modbam, creates a new modbam with probabilities set to
    /// 100% if a base modification is called or 0% if called canonical.
    CallMods(CallMods),
    /// Extract read-level base modification information from a modBAM into a
    /// tab-separated values table.
    Extract(ExtractMods),
    /// Repair MM and ML tags in one bam with the correct tags from another. To
    /// use this command, both modBAMs _must_ be sorted by read name. The
    /// "donor" modBAM's reads must be a superset of the acceptor's reads.
    /// Extra reads in the donor are allowed, and multiple reads with the
    /// same name (secondary, etc.) are allowed in the acceptor. Reads with
    /// an empty SEQ field cannot be repaired and will be rejected. Reads
    /// where there is an ambiguous alignment of the acceptor to the
    /// donor will be rejected (and logged). See the full documentation for
    /// details.
    Repair(RepairTags),
    /// Perform DMR test on a set of regions. Output a BED file of regions
    /// with the score column indicating the magnitude of the difference. Find
    /// the schema and description of fields can in the README as well as a
    /// description of the model and method. See subcommand help for
    /// additional details.
    #[clap(subcommand)]
    Dmr(BedMethylDmr),
    /// Tabulates double-stranded base modification patters (such as
    /// hemi-methylation) across genomic motif positions. This command
    /// produces a bedMethyl file, the schema can be found in the online
    /// documentation.
    PileupHemi(DuplexModBamPileup),
    /// Validate results from a set of mod-BAM files and associated BED files
    /// containing the ground truth modified base status at reference
    /// positions.
    Validate(ValidateFromModbam),
    /// Find sequence motifs in a bedMethyl pileup that are enriched for base
    /// modification.
    #[clap(subcommand)]
    Motif(EntryMotifs),
    /// Use a mod-BAM to calculate methylation entropy over genomic windows.
    Entropy(MethylationEntropy),
    /// Run localise
    Localise(EntryLocalize),
    Stats(EntryStats),
}

impl Commands {
    pub fn run(&self) -> AnyhowResult<()> {
        match self {
            Self::AdjustMods(x) => x.run(),
            Self::Pileup(x) => x.run(),
            Self::SampleProbs(x) => x.run(),
            Self::Summary(x) => x.run(),
            Self::UpdateTags(x) => x.run(),
            Self::CallMods(x) => x.run(),
            Self::Extract(x) => x.run(),
            Self::Repair(x) => x.run(),
            Self::Dmr(x) => x.run(),
            Self::PileupHemi(x) => x.run(),
            Self::Validate(x) => x.run(),
            Self::Motif(x) => x.run(),
            Self::Entropy(x) => x.run(),
            Self::Localise(x) => x.run(),
            Self::Stats(x) => x.run(),
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

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct Adjust {
    /// Input BAM file, can be a path to a file or one of
    /// `-` or `stdin` to specify a stream from standard input.
    in_bam: String,
    /// File path to new BAM file to be created. Can be a path to a file or one
    /// of `-` or `stdin` to specify a stream from standard output.
    out_bam: String,
    /// Output debug logs to file at this path.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Modified base code to ignore/remove, see
    /// https://samtools.github.io/hts-specs/SAMtags.pdf for details on
    /// the modified base codes.
    #[arg(long, conflicts_with = "convert")]
    ignore: Option<String>,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Fast fail, stop processing at the first invalid sequence record.
    /// Default behavior is to continue and report failed/skipped records
    /// at the end.
    #[arg(short, long = "ff", default_value_t = false)]
    fail_fast: bool,
    /// Convert one mod-tag to another, summing the probabilities together if
    /// the retained mod tag is already present.
    #[arg(long, action = clap::ArgAction::Append, num_args = 2, conflicts_with_all=["ignore", "filter_probs"])]
    convert: Option<Vec<String>>,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,
    /// Output SAM format instead of BAM.
    #[arg(long, default_value_t = false)]
    output_sam: bool,

    // filtering options
    // sampling args
    /// Filter out the lowest confidence base modification probabilities.
    #[arg(long, default_value_t = false)]
    filter_probs: bool,
    /// Sample approximately this many reads when estimating the filtering
    /// threshold. If alignments are present reads will be sampled evenly
    /// across aligned genome. If a region is specified, either with the
    /// --region option or the --sample-region option, then reads will be
    /// sampled evenly across the region given. This option is useful for
    /// large BAM files. In practice, 10-50 thousand reads is sufficient to
    /// estimate the model output distribution and determine the filtering
    /// threshold.
    #[arg(
        short = 'n',
        requires = "filter_probs",
        long,
        default_value_t = 10_042,
        hide_short_help = true
    )]
    num_reads: usize,
    /// Specify a region for sampling reads from when estimating the threshold
    /// probability. If this option is not provided, but --region is
    /// provided, the genomic interval passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long, requires = "filter_probs", hide_short_help = true)]
    sample_region: Option<String>,
    /// Interval chunk size to process concurrently when estimating the
    /// threshold probability, can be larger than the pileup processing
    /// interval.
    #[arg(
        long,
        requires = "filter_probs",
        default_value_t = 1_000_000,
        hide_short_help = true
    )]
    sampling_interval_size: u32,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[arg(short = 'p', requires = "filter_probs", long, default_value_t = 0.1)]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per primary base. A global
    /// filter threshold can be specified with by a decimal number (e.g.
    /// 0.75). Per-base thresholds can be specified by colon-separated
    /// values, for example C:0.75 specifies a threshold value of 0.75 for
    /// cytosine modification calls. Additional per-base thresholds can be
    /// specified by repeating the option: for example --filter-threshold
    /// C:0.75 --filter-threshold A:0.70 or specify a single base option
    /// and a default for all other bases with: --filter-threshold A:0.70
    /// --filter-threshold 0.9 will specify a threshold value of 0.70 for
    /// adenine and 0.9 for all other base modification calls.
    #[arg(
        long,
        conflicts_with="filter_percentile",
        requires="filter_probs",
        action = clap::ArgAction::Append,
        alias = "pass_threshold",
        hide_short_help = true,
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[arg(
        requires="filter_probs",
        long = "mod-threshold",
        action = clap::ArgAction::Append,
        hide_short_help = true,
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Only use base modification probabilities from bases that are aligned
    /// when estimating the filter threshold (i.e. ignore soft-clipped, and
    /// inserted bases).
    #[arg(
        long,
        default_value_t = false,
        hide_short_help = true,
        conflicts_with = "filter_percentile",
        requires = "filter_probs",
        hide_short_help = true
    )]
    only_mapped: bool,

    /// Hide the progress bar.
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
}

impl Adjust {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let io_threadpool = tpool::ThreadPool::new(self.threads as u32)?;
        let mut reader = get_serial_reader(self.in_bam.as_str())?;
        reader.set_thread_pool(&io_threadpool)?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);
        let mut bam_writer =
            get_bam_writer(&self.out_bam, &header, self.output_sam)?;
        bam_writer.set_thread_pool(&io_threadpool)?;

        let methods = if let Some(convert) = &self.convert {
            let convert = convert
                .iter()
                .map(|s| ModCodeRepr::parse(s.as_str()))
                .collect::<anyhow::Result<Vec<ModCodeRepr>>>()?;
            let mut conversions = HashMap::new();
            for chunk in convert.chunks(2) {
                debug_assert_eq!(chunk.len(), 2);
                let from: ModCodeRepr = chunk[0];
                let to: ModCodeRepr = chunk[1];
                conversions.entry(to).or_insert(HashSet::new()).insert(from);
            }
            for (to_code, from_codes) in conversions.iter() {
                info!(
                    "Converting {} to {}",
                    from_codes.iter().sorted().join(","),
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
            if let Some(ignore_base) = self.ignore.as_ref() {
                let ignore_base = ModCodeRepr::parse(ignore_base.as_str())?;
                info!(
                    "Removing mod base {} from {}, new bam {}",
                    ignore_base,
                    {
                        if using_stream(&self.in_bam) {
                            "stdin"
                        } else {
                            self.in_bam.as_str()
                        }
                    },
                    {
                        if using_stream(&self.out_bam) {
                            "stdout"
                        } else {
                            self.out_bam.as_str()
                        }
                    }
                );
                let method = CollapseMethod::ReDistribute(ignore_base);
                vec![method]
            } else {
                Vec::new()
            }
        };

        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;

        let methods = if edge_filter.is_none()
            && methods.is_empty()
            && !self.filter_probs
        {
            bail!(
                "no edge-filter, ignore, or convert was provided, no work to \
                 do. Provide --edge-filter, --ignore, --filter-probs, or \
                 --convert option to use `modkit adjust-mods`"
            )
        } else {
            debug_assert!(methods.len() <= 2 || !self.filter_probs);
            methods
        };

        let caller = if self.filter_probs {
            let per_mod_thresholds =
                if let Some(raw_per_mod_thresholds) = &self.mod_thresholds {
                    Some(parse_per_mod_thresholds(raw_per_mod_thresholds)?)
                } else {
                    None
                };

            let sampling_region = if let Some(raw_region) = &self.sample_region
            {
                info!("parsing sample region {raw_region}");
                Some(Region::parse_str(raw_region, &reader.header())?)
            } else {
                None
            };

            let caller = if let Some(raw_threshold) = &self.filter_threshold {
                parse_thresholds(raw_threshold, per_mod_thresholds)?
            } else {
                if using_stream(&self.in_bam) {
                    bail!(
                        "must specify all thresholds with --filter-threshold \
                         and (optionally) --mod-threshold when using stdin \
                         stream"
                    )
                }
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(self.threads)
                    .build()
                    .with_context(|| "failed to make threadpool")?;
                pool.install(|| {
                    get_threshold_from_options(
                        &Path::new(&self.in_bam).to_path_buf(),
                        self.threads,
                        self.sampling_interval_size,
                        None,
                        self.num_reads,
                        false,
                        self.filter_percentile,
                        None,
                        sampling_region.as_ref(),
                        per_mod_thresholds,
                        edge_filter.as_ref(),
                        methods.get(0),
                        None,
                        self.only_mapped,
                        self.suppress_progress,
                    )
                })?
            };
            Some(caller)
        } else {
            None
        };

        adjust_modbam(
            &mut reader,
            &mut bam_writer,
            &methods,
            caller.as_ref(),
            edge_filter.as_ref(),
            self.fail_fast,
            "Adjusting modBAM, records processed",
            self.suppress_progress,
            self.filter_probs,
        )?;
        Ok(())
    }
}

fn parse_percentiles(
    raw_percentiles: &str,
) -> Result<Vec<f32>, ParseFloatError> {
    if raw_percentiles.contains("..") {
        todo!("handle parsing ranges")
    } else {
        raw_percentiles.split(',').map(|x| x.parse::<f32>()).collect()
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct SampleModBaseProbs {
    /// Input BAM with modified base tags. If a index is found
    /// reads will be sampled evenly across the length of the
    /// reference sequence. Can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    in_bam: String,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Percentiles to calculate, a space separated list of floats.
    #[arg(short, long, default_value_t=String::from("0.1,0.5,0.9"))]
    percentiles: String,
    /// Directory to deposit result tables into. Required for model probability
    /// histogram output.
    #[arg(short = 'o', long)]
    out_dir: Option<PathBuf>,
    /// Label to prefix output files with.
    #[arg(long, requires = "out_dir")]
    prefix: Option<String>,
    /// Overwrite results if present.
    #[arg(long, requires = "out_dir", default_value_t = false)]
    force: bool,
    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[arg(long, hide_short_help = true)]
    ignore: Option<String>,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,

    // probability histogram options
    /// Output histogram of base modification prediction probabilities.
    #[arg(long = "hist", requires = "out_dir", default_value_t = false)]
    histogram: bool,

    /// Approximate maximum number of reads to use, especially recommended when
    /// using a large BAM without an index. If an indexed BAM is provided, the
    /// reads will be sampled evenly over the length of the aligned reference.
    /// If a region is passed with the --region option, they will be sampled
    /// over the genomic region. Actual number of reads used may deviate
    /// slightly from this number.
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
    /// Random seed for deterministic running, the default is
    /// non-deterministic, only used when no BAM index is provided.
    #[arg(short, requires = "sampling_frac", long)]
    seed: Option<u64>,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    /// Only used when sampling probs from an indexed bam.
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
    /// Only sample base modification probabilities that are aligned
    /// to the positions in this BED file. (alias: include-positions)
    #[arg(long, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Only use base modification probabilities that are aligned (i.e. ignore
    /// soft-clipped, and inserted bases).
    #[arg(long, default_value_t = false)]
    only_mapped: bool,
}

impl SampleModBaseProbs {
    fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if let Some(p) = self.out_dir.as_ref() {
            SampledProbs::check_files(
                p,
                self.prefix.as_ref(),
                self.force,
                self.histogram,
            )?;
        }

        let mut reader = get_serial_reader(&self.in_bam)?;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let region = if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Some(Region::parse_str(raw_region, reader.header())?)
        } else {
            None
        };
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;

        let (sample_frac, num_reads) = get_sampling_options(
            self.no_sampling,
            self.sampling_frac,
            self.num_reads,
        );

        let targets = get_targets(reader.header(), region.as_ref());
        let position_filter = self
            .include_bed
            .as_ref()
            .map(|bed_fp| {
                let chrom_to_tid = targets
                    .iter()
                    .map(|reference_record| {
                        (reference_record.name.as_str(), reference_record.tid)
                    })
                    .collect::<HashMap<&str, u32>>();
                StrandedPositionFilter::from_bed_file(
                    bed_fp,
                    &chrom_to_tid,
                    self.suppress_progress,
                )
            })
            .transpose()?;

        let collapse_method =
            if let Some(raw_mod_code_to_ignore) = self.ignore.as_ref() {
                let mod_code_to_ignore =
                    ModCodeRepr::parse(raw_mod_code_to_ignore)?;
                Some(CollapseMethod::ReDistribute(mod_code_to_ignore))
            } else {
                None
            };

        let desired_percentiles = parse_percentiles(&self.percentiles)
            .with_context(|| {
                format!("failed to parse percentiles: {}", &self.percentiles)
            })?;

        pool.install(|| {
            let read_ids_to_base_mod_calls = if using_stream(&self.in_bam) {
                reader.set_threads(self.threads)?;
                let record_sampler = RecordSampler::new_from_options(
                    sample_frac,
                    num_reads,
                    self.seed,
                );
                let read_ids_to_base_mod_probs =
                    ReadIdsToBaseModProbs::process_records(
                        reader.records(),
                        !self.suppress_progress,
                        record_sampler,
                        collapse_method.as_ref(),
                        edge_filter.as_ref(),
                        position_filter.as_ref(),
                        self.only_mapped || position_filter.is_some(),
                        false,
                        None,
                        None,
                    )?;
                debug!("sampled {} records", read_ids_to_base_mod_probs.len());
                read_ids_to_base_mod_probs
            } else {
                drop(reader);
                get_sampled_read_ids_to_base_mod_probs::<ReadIdsToBaseModProbs>(
                    &Path::new(&self.in_bam).to_path_buf(),
                    self.threads,
                    self.interval_size,
                    sample_frac,
                    num_reads,
                    self.seed,
                    region.as_ref(),
                    collapse_method.as_ref(),
                    edge_filter.as_ref(),
                    position_filter.as_ref(),
                    self.only_mapped || position_filter.is_some(),
                    self.suppress_progress,
                )?
            };

            let histograms = if self.histogram {
                Some(
                    read_ids_to_base_mod_calls
                        .get_per_mod_histograms(self.suppress_progress),
                )
            } else {
                None
            };

            let percentiles = read_ids_to_base_mod_calls
                .mle_probs_per_base(self.suppress_progress)
                .into_iter()
                .map(|(canonical_base, mut probs)| {
                    Percentiles::new(&mut probs, &desired_percentiles)
                        .with_context(|| {
                            format!(
                                "failed to calculate threshold for base {}",
                                canonical_base.char()
                            )
                        })
                        .map(|percs| (canonical_base, percs))
                })
                .collect::<AnyhowResult<HashMap<DnaBase, Percentiles>>>()?;

            let sampled_probs =
                SampledProbs::new(histograms, percentiles, self.prefix.clone());

            let mut writer: Box<dyn OutWriter<SampledProbs>> =
                if let Some(p) = &self.out_dir {
                    if !p.exists() {
                        info!("creating directory at {p:?}");
                        std::fs::create_dir_all(p)?;
                    }
                    sampled_probs.check_path(p, self.force)?;
                    Box::new(MultiTableWriter::new(p.clone()))
                } else {
                    Box::new(TsvWriter::new_stdout(None))
                };

            writer.write(sampled_probs)?;

            Ok(())
        })
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct ModSummarize {
    /// Input modBam, can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    in_bam: String,
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
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

    // sampling options
    /// Approximate maximum number of reads to use, especially recommended when
    /// using a large BAM without an index. If an indexed BAM is provided, the
    /// reads will be sampled evenly over the length of the aligned reference.
    /// If a region is passed with the --region option, they will be sampled
    /// over the genomic region. Actual number of reads used may deviate
    /// slightly from this number.
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
    /// No sampling, use all of the reads to calculate the filter thresholds
    /// and generating the summary.
    #[arg(long, group = "sampling_options", default_value_t = false)]
    no_sampling: bool,
    /// Sets a random seed for deterministic running (when using
    /// --sample-frac), the default is non-deterministic, only used when no
    /// BAM index is provided.
    #[arg(short, requires = "sampling_frac", long)]
    seed: Option<u64>,

    // threshold options
    /// Do not perform any filtering, include all base modification calls in
    /// the summary. See filtering.md for details on filtering.
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence base modification calls.
    #[arg(group = "thresholds", short = 'p', long, default_value_t = 0.1)]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per-base. Global filter
    /// threshold can be specified with by a decimal number (e.g. 0.75).
    /// Per-base thresholds can be specified by colon-separated values, for
    /// example C:0.75 specifies a threshold value of 0.75 for cytosine
    /// modification calls. Additional per-base thresholds can be specified
    /// by repeating the option: for example --filter-threshold C:0.75
    /// --filter-threshold A:0.70 or specify a single base option and a
    /// default for all other bases with: --filter-threshold A:0.70
    /// --filter-threshold 0.9 will specify a threshold value of 0.70 for
    /// adenine and 0.9 for all other base modification calls.
    #[arg(
        long,
        group = "thresholds",
        action = clap::ArgAction::Append
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[arg(
    long,
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[arg(long, group = "combine_args", hide_short_help = true)]
    ignore: Option<String>,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,
    /// Only summarize base modification probabilities that are aligned
    /// to the positions in this BED file. (alias: include-positions)
    #[arg(long, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Only use base modification probabilities that are aligned (i.e. ignore
    /// soft-clipped, and inserted bases).
    #[arg(long, default_value_t = false)]
    only_mapped: bool,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// When using regions, interval chunk size in base pairs to process
    /// concurrently. Smaller interval chunk sizes will use less memory but
    /// incur more overhead.
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
}

impl ModSummarize {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let mut reader = get_serial_reader(&self.in_bam)?;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let region = self
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, reader.header()))
            .transpose()?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;

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

        let position_filter = self
            .include_bed
            .as_ref()
            .map(|bed_fp| {
                let targets = get_targets(reader.header(), region.as_ref());
                let chrom_to_tid = targets
                    .iter()
                    .map(|reference_record| {
                        (reference_record.name.as_str(), reference_record.tid)
                    })
                    .collect::<HashMap<&str, u32>>();
                StrandedPositionFilter::from_bed_file(
                    bed_fp,
                    &chrom_to_tid,
                    self.suppress_progress,
                )
            })
            .transpose()?;

        let filter_thresholds = if let Some(raw_thresholds) =
            &self.filter_threshold
        {
            info!("parsing user defined thresholds");
            Some(parse_thresholds(raw_thresholds, per_mod_thresholds.clone())?)
        } else if self.no_filtering {
            info!("not performing filtering");
            Some(MultipleThresholdModCaller::new_passthrough())
        } else {
            None
        };

        let collapse_method =
            if let Some(raw_mod_code_to_ignore) = self.ignore.as_ref() {
                let mod_code_to_ignore =
                    ModCodeRepr::parse(raw_mod_code_to_ignore)?;
                Some(CollapseMethod::ReDistribute(mod_code_to_ignore))
            } else {
                None
            };

        let mod_summary = pool.install(|| {
            let read_ids_to_base_mod_calls = if using_stream(&self.in_bam) {
                reader.set_threads(self.threads)?;
                let record_sampler = RecordSampler::new_from_options(
                    sample_frac,
                    num_reads,
                    self.seed,
                );
                let read_ids_to_base_mod_probs =
                    ReadIdsToBaseModProbs::process_records(
                        reader.records(),
                        !self.suppress_progress,
                        record_sampler,
                        collapse_method.as_ref(),
                        edge_filter.as_ref(),
                        position_filter.as_ref(),
                        self.only_mapped || position_filter.is_some(),
                        false,
                        None,
                        None,
                    )?;
                debug!("sampled {} records", read_ids_to_base_mod_probs.len());
                read_ids_to_base_mod_probs
            } else {
                drop(reader);
                get_sampled_read_ids_to_base_mod_probs::<ReadIdsToBaseModProbs>(
                    &Path::new(&self.in_bam).to_path_buf(),
                    self.threads,
                    self.interval_size,
                    sample_frac,
                    num_reads,
                    self.seed,
                    region.as_ref(),
                    collapse_method.as_ref(),
                    edge_filter.as_ref(),
                    position_filter.as_ref(),
                    self.only_mapped || position_filter.is_some(),
                    self.suppress_progress,
                )?
            };
            let threshold_caller = if let Some(ft) = filter_thresholds {
                // filter thresholds provided, use those
                ft
            } else {
                // calculate the filter thresholds at the requested percentile
                let pct = (self.filter_percentile * 100f32).floor();
                info!("calculating threshold at {pct}(th) percentile");
                calc_thresholds_per_base(
                    &read_ids_to_base_mod_calls,
                    self.filter_percentile,
                    None,
                    per_mod_thresholds,
                    self.suppress_progress,
                )?
            };

            sampled_reads_to_summary(
                read_ids_to_base_mod_calls,
                &threshold_caller,
                region.as_ref(),
                self.suppress_progress,
            )
        })?;

        let mut writer: Box<dyn OutWriter<ModSummary>> = if self.tsv_format {
            Box::new(TsvWriter::new_stdout(None))
        } else {
            Box::new(TableWriter::new())
        };
        writer.write(mod_summary)?;
        Ok(())
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
#[allow(non_camel_case_types)]
enum ModMode {
    #[clap(alias = "ambiguous")]
    explicit,
    implicit,
}

impl ModMode {
    fn to_skip_mode(self) -> SkipMode {
        match self {
            Self::explicit => SkipMode::Explicit,
            Self::implicit => SkipMode::ProbModified,
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct Update {
    /// BAM to update modified base tags in. Can be a path to a file or one of
    /// `-` or `stdin` to specify a stream from standard input.
    in_bam: String,
    /// File to new BAM file to be created or one of `-` or `stdin` to specify
    /// a stream from standard output.
    out_bam: String,
    /// Mode, change mode to this value, options {'explicit', 'implicit'}.
    /// See spec at: https://samtools.github.io/hts-specs/SAMtags.pdf.
    /// 'explicit' ('?') means residues without modification
    /// probabilities will not be assumed canonical or modified. 'implicit'
    /// means residues without explicit modification probabilities are
    /// assumed to be canonical.
    #[arg(short, long, value_enum)]
    mode: Option<ModMode>,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Don't add implicit canonical calls. This flag is important when
    /// converting from one of the implicit modes ( `.` or `""`) to
    /// explicit mode (`?`). By passing this flag, the bases without
    /// associated base modification probabilities will not be assumed to
    /// be canonical. No base modification probability will be written for
    /// these bases, meaning there is no information. The mode will
    /// automatically be set to the explicit mode `?`.
    #[arg(long, default_value_t = false)]
    no_implicit_probs: bool,
    /// Output debug logs to file at this path.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Output SAM format instead of BAM.
    #[arg(long, default_value_t = false)]
    output_sam: bool,
}

fn update_mod_tags(
    mut record: bam::Record,
    no_implicit_calls: bool,
    new_mode: SkipMode,
) -> CliResult<bam::Record> {
    let mod_base_info = ModBaseInfo::new_from_record(&record)?;
    let mm_style = mod_base_info.mm_style;
    let ml_style = mod_base_info.ml_style;

    let mut mm_agg = String::new();
    let mut ml_agg = Vec::new();

    let (converters, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();
    for (base, strand, mut seq_pos_mod_probs) in mod_prob_iter {
        let converter = converters.get(&base).unwrap();
        if no_implicit_calls && new_mode == SkipMode::Explicit {
            seq_pos_mod_probs = seq_pos_mod_probs.remove_implicit_probs();
        } else {
            seq_pos_mod_probs.set_skip_mode(new_mode);
        }
        let (mm, mut ml) =
            format_mm_ml_tag(seq_pos_mod_probs, strand, converter);
        mm_agg.push_str(&mm);
        ml_agg.extend_from_slice(&mut ml);
    }
    record.remove_aux(mm_style.as_bytes()).expect("failed to remove MM tag");
    record.remove_aux(ml_style.as_bytes()).expect("failed to remove ML tag");
    let mm = Aux::String(&mm_agg);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml_agg;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record.push_aux(MM_TAGS[0].as_bytes(), mm).expect("failed to add MM tag");
    record.push_aux(ML_TAGS[0].as_bytes(), ml).expect("failed to add ML tag");

    Ok(record)
}

impl Update {
    fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let threads = self.threads;
        let mut reader = get_serial_reader(&self.in_bam)?;
        reader.set_threads(threads)?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);

        let mut out_bam =
            get_bam_writer(&self.out_bam, &header, self.output_sam)?;
        let spinner = get_ticker();

        spinner.set_message("Updating ModBAM");
        let mut total = 0usize;
        let mut total_failed = 0usize;
        let mut total_skipped = 0usize;

        let to_mode = if let Some(input_mode) = self.mode {
            let skip_mode = input_mode.to_skip_mode();
            if self.no_implicit_probs && skip_mode != SkipMode::Explicit {
                bail!("cannot change to {input_mode:?} and skip implicit probs")
            }
            skip_mode
        } else {
            if self.no_implicit_probs {
                info!("implicit canonical probs will not be present in output");
                info!("setting mode to explicit, `?`");
                SkipMode::Explicit
            } else {
                info!("mode will be set to prob-modified, '.'");
                SkipMode::ProbModified
            }
        };

        for (i, result) in reader.records().enumerate() {
            if let Ok(record) = result {
                let record_name = util::get_query_name_string(&record)
                    .unwrap_or("???".to_owned());
                match update_mod_tags(record, self.no_implicit_probs, to_mode) {
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
                            total = i + 1;
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
#[command(arg_required_else_help = true)]
pub struct CallMods {
    // running args
    /// Input BAM, may be sorted and have associated index available. Can be a
    /// path to a file or one of `-` or `stdin` to specify a stream from
    /// standard input.
    in_bam: String,
    /// Output BAM, can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    out_bam: String,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    // /// Process only the specified region of the BAM when performing
    // transformation. /// Format should be <chrom_name>:<start>-<end> or
    // <chrom_name>. #[arg(long)] todo(arand)
    // region: Option<String>,
    /// Fast fail, stop processing at the first invalid sequence record.
    /// Default behavior is to continue and report failed/skipped records
    /// at the end.
    #[arg(long = "ff", default_value_t = false)]
    fail_fast: bool,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

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
    /// Sample approximately this many reads when estimating the filtering
    /// threshold. If alignments are present reads will be sampled evenly
    /// across aligned genome. If a region is specified, either with the
    /// --region option or the --sample-region option, then reads will be
    /// sampled evenly across the region given. This option is useful for
    /// large BAM files. In practice, 10-50 thousand reads is sufficient to
    /// estimate the model output distribution and determine the filtering
    /// threshold.
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Sample this fraction of the reads when estimating the
    /// filter-percentile. In practice, 50-100 thousand reads is sufficient
    /// to estimate the model output distribution and determine the
    /// filtering threshold. See filtering.md for details on filtering.
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Set a random seed for deterministic running, the default is
    /// non-deterministic, only used when no BAM index is provided.
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Specify a region for sampling reads from when estimating the threshold
    /// probability. If this option is not provided, but --region is
    /// provided, the genomic interval passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    sample_region: Option<String>,
    /// Interval chunk size to process concurrently when estimating the
    /// threshold probability, can be larger than the pileup processing
    /// interval.
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,

    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per primary base. A global
    /// filter threshold can be specified with by a decimal number (e.g.
    /// 0.75). Per-base thresholds can be specified by colon-separated
    /// values, for example C:0.75 specifies a threshold value of 0.75 for
    /// cytosine modification calls. Additional per-base thresholds can be
    /// specified by repeating the option: for example --filter-threshold
    /// C:0.75 --filter-threshold A:0.70 or specify a single base option
    /// and a default for all other bases with: --filter-threshold A:0.70
    /// --filter-threshold 0.9 will specify a threshold value of 0.70 for
    /// adenine and 0.9 for all other base modification calls.
    #[arg(
    long,
    group = "thresholds",
    action = clap::ArgAction::Append,
    alias = "pass_threshold"
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[arg(
    long = "mod-threshold",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Don't filter base modification calls, assign each base modification to
    /// the highest probability prediction.
    #[arg(long, default_value_t = false)]
    no_filtering: bool,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,
    /// Output SAM format instead of BAM.
    #[arg(long, default_value_t = false)]
    output_sam: bool,
}

impl CallMods {
    pub fn run(&self) -> AnyhowResult<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let io_threadpool = tpool::ThreadPool::new(self.threads as u32)?;
        let mut reader = get_serial_reader(&self.in_bam)?;
        reader.set_thread_pool(&io_threadpool)?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);
        let mut bam_writer =
            get_bam_writer(&self.out_bam, &header, self.output_sam)?;
        bam_writer.set_thread_pool(&io_threadpool)?;

        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;

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
            if using_stream(&self.in_bam) {
                bail!(
                    "must specify all thresholds with --filter-threshold and \
                     (optionally) --mod-threshold when using stdin stream"
                )
            }
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .with_context(|| "failed to make threadpool")?;
            pool.install(|| {
                get_threshold_from_options(
                    &Path::new(&self.in_bam).to_path_buf(),
                    self.threads,
                    self.sampling_interval_size,
                    self.sampling_frac,
                    self.num_reads,
                    self.no_filtering,
                    self.filter_percentile,
                    self.seed,
                    sampling_region.as_ref(),
                    per_mod_thresholds,
                    edge_filter.as_ref(),
                    None,
                    None,
                    false,
                    self.suppress_progress,
                )
            })?
        };

        adjust_modbam(
            &mut reader,
            &mut bam_writer,
            &[],
            Some(&caller),
            edge_filter.as_ref(),
            self.fail_fast,
            "Calling Mods, records processed",
            self.suppress_progress,
            false,
        )?;

        Ok(())
    }
}
