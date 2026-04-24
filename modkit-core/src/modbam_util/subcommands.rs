use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::num::ParseFloatError;
use std::ops::AddAssign;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::{Args, Subcommand, ValueEnum};
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressDrawTarget};
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{FetchDefinition, HeaderView, Read};
use rust_htslib::{bam, tpool};

use itertools::Itertools;
use modkit_logging::init_logging;
use rustc_hash::FxHashMap;

use crate::adjust::adjust_modbam;
use crate::command_utils::{
    get_bam_writer, get_motif_lookup_from_parts, get_serial_reader,
    get_threshold_from_options, parse_edge_filter_input, parse_forward_motifs,
    parse_per_mod_thresholds, parse_raw_motifs, parse_thresholds, using_stream,
};
use crate::errs::{MkError, MkResult};
use crate::interval_chunks::{
    ChromCoordinatesFeeder, ReferenceIntervalBatchesFeeder, TotalLength,
    WithPrevEnd,
};
use crate::mod_bam::{
    format_mm_ml_tag, CollapseMethod, ModBaseInfo, SkipMode, ML_TAGS, MM_TAGS,
};
use crate::mod_base_code::{
    BaseAndState, BaseState, DnaBase, ModCodeRepr, ProbHistogram,
};
use crate::modbam_util::check_tags::ModTagViews;
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::ReadIdsToBaseModProbs;
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::reads_sampler::{
    get_sampled_read_ids_to_base_mod_probs, sample_reads_from_interval,
};
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::repair_tags::RepairTags;
use crate::sample_probs::{
    run_extract_probs_workers, AlignedBaseAndModArgmaxProbs,
    AlignedBaseArgmaxProbs, BaseAndModArgmaxProbs, BaseArgmaxProbs,
    ExtractProbsWorker, ProbsExtractor, QualHist, RegionMleProbs,
};
use crate::summarize::{sampled_reads_to_summary, ModSummary};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::calc_thresholds_per_base;
use crate::util::{
    add_modkit_pg_records, filter_reference_records, format_errors_table,
    get_master_progress_bar, get_subroutine_progress_bar, get_targets,
    get_ticker, ReferenceRecord, Region, DEFAULT_NUM_READS,
};
use crate::writers::{
    MultiTableWriter, OutWriter, SampledProbs, TableWriter, TsvWriter,
};

#[derive(Subcommand)]
pub enum EntryModBam {
    #[command(name = "check-tags")]
    CheckTags(EntryCheckTags),
    // SampleReads(EntrySampleReads),
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
}

impl EntryModBam {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryModBam::CheckTags(x) => x.run(),
            EntryModBam::AdjustMods(adjust) => adjust.run(),
            EntryModBam::UpdateTags(update) => update.run(),
            EntryModBam::SampleProbs(sample_mod_base_probs) => {
                sample_mod_base_probs.run()
            }
            EntryModBam::Summary(mod_summarize) => mod_summarize.run(),
            EntryModBam::CallMods(call_mods) => call_mods.run(),
            EntryModBam::Repair(repair) => repair.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryCheckTags {
    /// Input modBam, can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    in_bam: String,
    /// Don't exit 1 when invalid records are found in the input.
    #[arg(long, default_value_t = false)]
    permissive: bool,
    /// Write output tables into this directory. The directory will be created
    /// if it doesn't exist.
    #[clap(help_heading = "IO Options")]
    #[arg(short = 'o', long)]
    out_dir: Option<PathBuf>,
    /// Force overwrite of previous output
    #[clap(help_heading = "IO Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// Prefix output files with this string.
    #[clap(help_heading = "IO Options")]
    #[arg(long)]
    prefix: Option<String>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Perform a linear scan of the modBAM even if the index is found.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    ignore_index: bool,
    /// When using regions, interval chunk size in base pairs to process
    /// concurrently. Smaller interval chunk sizes will use less memory but
    /// incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'i', long, default_value_t = 5_000_000)]
    interval_size: u32,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

    /// Approximate maximum number of reads to use, especially recommended when
    /// using a large BAM without an index. If an indexed BAM is provided, the
    /// reads will be sampled evenly over the length of the aligned reference.
    /// If a region is passed with the --region option, they will be sampled
    /// over the region. Actual number of reads used may deviate
    /// slightly from this number.
    #[clap(help_heading = "Selection Options")]
    #[arg(short = 'n', long)]
    num_reads: Option<usize>,
    /// Check tags on non-primary alignments as well. Keep in mind this
    /// may incur a double-counting of the read with its primary mapping.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, default_value_t = false)]
    allow_non_primary: bool,
    /// Only check alignments that are mapped.
    #[clap(help_heading = "Selection Options")]
    #[arg(
        long = "mapped-only",
        alias = "only-mapped",
        default_value_t = false
    )]
    only_mapped: bool,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,
    /// Use the first N reads without sampling. Shorthand for --ignore-index
    /// --num-reads N.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, conflicts_with_all = ["num_reads", "ignore_index"])]
    head: Option<usize>,
}

impl EntryCheckTags {
    pub fn run(&self) -> anyhow::Result<()> {
        use super::check_tags::output_filenames as ofn;
        let _handle = init_logging(self.log_filepath.as_ref());
        let mut reader = get_serial_reader(&self.in_bam)?;
        if let Some(out_d) = self.out_dir.as_ref() {
            if !out_d.exists() {
                info!("creating directory at {out_d:?}");
                std::fs::create_dir_all(out_d)?;
            }
            for fp in ofn::filenames {
                let out_fn = if let Some(p) = self.prefix.as_ref() {
                    out_d.join(format!("{p}_{fp}"))
                } else {
                    out_d.join(fp)
                };
                if out_fn.exists() && !self.force {
                    bail!("refusing to overwrite {fp:?}, use --force");
                }
            }
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let region = self
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, reader.header()))
            .transpose()?;
        if region.is_some() && using_stream(&self.in_bam) {
            bail!("cannot use --region and stream from stdin");
        }
        let schedule = match (self.num_reads, using_stream(&self.in_bam)) {
            (_, true) | (None, false) => None,
            (Some(num_reads), false) => {
                match bam::IndexedReader::from_path(&self.in_bam) {
                    Ok(_) => {
                        let bam_fp = Path::new(&self.in_bam).to_path_buf();
                        if !bam_fp.exists() {
                            bail!("failed to find ${bam_fp:?}");
                        }
                        Some(SamplingSchedule::from_num_reads(
                            &self.in_bam,
                            num_reads,
                            region.as_ref(),
                            None,
                            !self.only_mapped,
                        )?)
                    }
                    Err(_) => {
                        warn!(
                            "didn't find BAM index, using first (head) \
                             {num_reads} reads"
                        );
                        None
                    }
                }
            }
        };

        let indexed_reader = bam::IndexedReader::from_path(&self.in_bam);
        let has_index = match (indexed_reader, self.region.as_ref()) {
            (Err(e), Some(_)) => {
                bail!("cannot use --region, failed to get BAM index, {e}")
            }
            (Err(_e), None) => false,
            (Ok(_), _) => true,
        };

        let linear_scan = self.ignore_index
            || self.head.is_some()
            || using_stream(&self.in_bam)
            || !has_index;
        let tag_views = pool.install(|| {
            if linear_scan {
                reader.set_threads(self.threads)?;

                let record_sampler = if let Some(&n_reads) =
                    self.num_reads.as_ref().or(self.head.as_ref())
                {
                    info!("checking tags on first {n_reads} reads");
                    RecordSampler::new_num_reads(n_reads)
                } else {
                    RecordSampler::new_passthrough()
                };
                ModTagViews::process_records(
                    reader.records(),
                    !self.suppress_progress,
                    record_sampler,
                    None,
                    None,
                    None,
                    self.only_mapped,
                    self.allow_non_primary,
                    None,
                    None,
                )
            } else {
                let bam_fp = Path::new(&self.in_bam);
                if !bam_fp.exists() {
                    bail!("failed to find modBAM at {bam_fp:?}");
                }

                let reference_records =
                    get_targets(reader.header(), region.as_ref());
                let feeder = ReferenceIntervalBatchesFeeder::new(
                    reference_records,
                    (self.threads as f32 * 1.5f32).floor() as usize,
                    self.interval_size,
                    false,
                    None,
                    None,
                )?;
                pool.install(|| {
                    self.run_check_tags_indexed(
                        bam_fp.to_path_buf(),
                        feeder,
                        schedule,
                    )
                })
            }
        })?;

        tag_views.report(
            self.out_dir.as_ref(),
            self.prefix.as_ref(),
            self.force,
            self.permissive,
        )?;
        Ok(())
    }

    fn run_check_tags_indexed(
        &self,
        bam_fp: PathBuf,
        feeder: ReferenceIntervalBatchesFeeder,
        schedule: Option<SamplingSchedule>,
    ) -> anyhow::Result<ModTagViews> {
        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(ProgressDrawTarget::hidden());
        }
        let prog_length = feeder.total_length() as usize;
        let master_progress = mpb.add(get_master_progress_bar(prog_length));
        master_progress.set_message("genome positions");

        let feeder = feeder.map(|x| x.unwrap()).with_prev_end();
        let (snd, rcv) = crossbeam::channel::bounded(1000);
        let allow_non_primary = self.allow_non_primary;
        let only_mapped = self.only_mapped;
        let bam_fp1 = bam_fp.clone();
        let with_progress = !self.suppress_progress;
        let n_reads = self.num_reads;
        let threads = self.threads;
        std::thread::spawn(move || {
            let mut records_used = 0usize;
            for super_batch in feeder {
                let total_batch_length =
                    super_batch.iter().map(|c| c.total_length()).sum::<u64>();
                let batch_progress =
                    mpb.add(get_subroutine_progress_bar(super_batch.len()));
                batch_progress.set_message("batch progress");
                let results = super_batch
                    .into_par_iter()
                    .progress_with(batch_progress)
                    .map(|batch| {
                        batch
                            .0
                            .into_par_iter()
                            .filter(|cc| {
                                schedule
                                    .as_ref()
                                    .map(|s| s.chrom_has_reads(cc.chrom_tid()))
                                    .unwrap_or(true)
                            })
                            .map(|cc| {
                                let record_sampler = schedule
                                    .as_ref()
                                    .map(|s| {
                                        s.get_record_sampler(
                                            cc.chrom_tid(),
                                            total_batch_length as u32,
                                            cc.start_pos(),
                                            cc.end_pos(),
                                        )
                                    })
                                    .unwrap_or_else(|| {
                                        RecordSampler::new_passthrough()
                                    });
                                sample_reads_from_interval::<ModTagViews>(
                                    &bam_fp1,
                                    cc.chrom_tid(),
                                    cc.start_pos(),
                                    cc.end_pos(),
                                    cc.prev_end(),
                                    record_sampler,
                                    None,
                                    None,
                                    None,
                                    true,
                                    allow_non_primary,
                                    None,
                                )
                            })
                            .collect::<anyhow::Result<Vec<ModTagViews>>>()
                    })
                    .collect::<Vec<anyhow::Result<Vec<ModTagViews>>>>();
                let records_accumulated = results
                    .par_iter()
                    .filter_map(|res| {
                        if let Ok(views) = res {
                            Some(
                                views
                                    .iter()
                                    .map(|x| x.num_reads())
                                    .sum::<usize>(),
                            )
                        } else {
                            None
                        }
                    })
                    .sum::<usize>();

                records_used = records_used.saturating_add(records_accumulated);

                match snd.send(results) {
                    Ok(_) => {}
                    Err(e) => {
                        debug!("failed to send on channel, {e}");
                    }
                }
                master_progress.inc(total_batch_length);
            }
            if !only_mapped {
                debug!("processing unmapped reads");
                let sampler = n_reads
                    .map(|nr| nr.checked_sub(records_used).unwrap_or(0))
                    .map(|n| {
                        debug!("processing {n} unmapped reads");
                        RecordSampler::new_num_reads(n)
                    })
                    .unwrap_or_else(|| RecordSampler::new_passthrough());
                let reader = bam::IndexedReader::from_path(&bam_fp)
                    .and_then(|mut reader| {
                        reader.fetch(FetchDefinition::Unmapped).map(|_| reader)
                    })
                    .and_then(|mut reader| {
                        reader.set_threads(threads).map(|_| reader)
                    });
                match reader {
                    Ok(mut reader) => {
                        let unmapped_tag_views = ModTagViews::process_records(
                            reader.records(),
                            with_progress,
                            sampler,
                            None,
                            None,
                            None,
                            false,
                            true,
                            None,
                            None,
                        );
                        if let Ok(unmapped) = unmapped_tag_views {
                            match snd.send(vec![Ok(vec![unmapped])]) {
                                Ok(_) => {}
                                Err(e) => {
                                    debug!("failed to send on channel, {e}");
                                }
                            }
                        }
                    }
                    Err(e) => {
                        debug!("failed to collect unmapped records, {e}");
                    }
                }
            }
        });

        let tag_views = rcv
            .iter()
            .par_bridge()
            .inspect(|x| {
                let n_errs = x
                    .iter()
                    .filter(|y| y.is_err())
                    // .inspect(|x| {
                    //     dbg!(&x);
                    // })
                    .count();

                if n_errs > 0 {
                    debug!("{n_errs} batches failed");
                }
            })
            .map(|xs| {
                xs.into_iter()
                    .filter_map(|x| x.ok())
                    .flatten()
                    .collect::<Vec<ModTagViews>>()
            })
            .fold(
                || ModTagViews::default(),
                |mut a, b| {
                    b.into_iter().for_each(|x| a.op_mut(x));
                    a
                },
            )
            .reduce(|| ModTagViews::default(), |a, b| a.op(b));
        Ok(tag_views)
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntrySampleReads {
    in_bam: String,
}

impl EntrySampleReads {
    pub fn run(&self) -> anyhow::Result<()> {
        todo!()
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
    #[clap(help_heading = "Output Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,

    /// Modified base code to ignore/remove, see
    /// https://samtools.github.io/hts-specs/SAMtags.pdf for details on
    /// the modified base codes.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, conflicts_with = "convert")]
    ignore: Option<String>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Fast fail, stop processing at the first invalid sequence record.
    /// Default behavior is to continue and report failed/skipped records
    /// at the end.
    #[arg(short, long = "ff", default_value_t = false)]
    fail_fast: bool,
    /// Convert one mod-tag to another, summing the probabilities together if
    /// the retained mod tag is already present.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, action = clap::ArgAction::Append, num_args = 2, conflicts_with_all=["ignore", "filter_probs"])]
    convert: Option<Vec<String>>,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,
    /// Output SAM format instead of BAM.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    output_sam: bool,

    // filtering options
    // sampling args
    /// Filter out the lowest confidence base modification probabilities.
    #[clap(help_heading = "Selection Options")]
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
    #[clap(help_heading = "Sampling Options")]
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
    #[clap(help_heading = "Sampling Options")]
    #[arg(long, requires = "filter_probs", hide_short_help = true)]
    sample_region: Option<String>,
    /// Interval chunk size to process concurrently when estimating the
    /// threshold probability, can be larger than the pileup processing
    /// interval.
    #[clap(help_heading = "Sampling Options")]
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
    #[clap(help_heading = "Filtering Options")]
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
    #[clap(help_heading = "Filtering Options")]
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
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        requires="filter_probs",
        long = "mod-threshold",
        alias = "mod-thresholds",
        action = clap::ArgAction::Append,
        hide_short_help = true,
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Only use base modification probabilities from bases that are aligned
    /// when estimating the filter threshold (i.e. ignore soft-clipped, and
    /// inserted bases).
    #[clap(help_heading = "Selection Options")]
    #[arg(
        long = "mapped-only",
        alias = "only-mapped",
        default_value_t = false,
        hide_short_help = true,
        conflicts_with = "filter_percentile",
        requires = "filter_probs",
        hide_short_help = true
    )]
    only_mapped: bool,

    /// Filter out any base modification call that isn't part of a basecall
    /// sequence motif. This argument can be passed multiple times. Format
    /// is <motif_sequence> <offset>. For example the argument to match CpG
    /// dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
    /// would be `--motif CGCG 2`. Single bases can be used as motifs
    /// to keep only base modification calls for a specific primary base,
    /// for example `--motif C 0`.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, action = clap::ArgAction::Append, num_args = 2)]
    motif: Option<Vec<String>>,
    /// Shorthand for --motif CG 0.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, default_value_t = false)]
    cpg: bool,
    /// Discard base modification calls that match the provided motifs (instead
    /// of keeping them).
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, requires = "motif", default_value_t = false)]
    discard_motifs: bool,

    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
}

impl Adjust {
    pub fn run(&self) -> anyhow::Result<()> {
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

        dbg!(&self.cpg);
        let motifs = parse_forward_motifs(&self.motif, self.cpg)?;
        if let Some(ms) = motifs.as_ref() {
            let patterns = ms.iter().map(|x| x.as_str()).join(",");
            info!("filtering base modification calls to patterns: {patterns}");
        }

        let have_motifs =
            motifs.as_ref().map(|x| !x.is_empty()).unwrap_or(false);

        if edge_filter.is_none()
            && methods.is_empty()
            && !self.filter_probs
            && !have_motifs
        {
            bail!(
                "no edge-filter, ignore, motifs, or convert was provided, no \
                 work to do. Provide --edge-filter, --ignore, --filter-probs, \
                 --motif, or --convert option to use `modkit adjust-mods`"
            )
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
            &motifs,
            self.discard_motifs,
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
    /// Ignore the BAM index (if it exists) and default to a serial scan of the
    /// BAM, otherwise if an index is found, multiple workers will be used to
    /// read the BAM in parallel. Conflicts with --threads since there will
    /// only be one reader.
    #[arg(
        long,
        conflicts_with = "threads",
        default_value_t = false,
        hide_short_help = true
    )]
    pub ignore_index: bool,
    /// Reference sequence in FASTA format. (alias: 'ref')
    #[arg(long = "reference", alias = "ref", short = 'r')]
    reference_fasta: Option<PathBuf>,
    /// Preload the reference sequences, useful when working with many, short
    /// reference sequences such as a transcriptome.
    #[arg(long, default_value_t = false, requires = "reference_fasta")]
    preload_references: bool,
    /// Respect soft masking in the reference FASTA.
    #[arg(
        long,
        short = 'k',
        requires = "reference_fasta",
        default_value_t = false,
        hide_short_help = true
    )]
    mask: bool,
    /// Number of workers to run concurrently. When an aligned, indexed BAM is
    /// found, multiple workers (threads) will be spawned to read disjoint
    /// sections of the BAM in parallel.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Number of IO/decompression threads to use, these are shared across
    /// workers (if an indexed BAM is found) or used for a single reader if
    /// --ignore-index is passed, stdin is used, or an index is not found.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: u32,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Percentiles to calculate, a space separated list of floats (for
    /// example: "10,20,90").
    #[clap(help_heading = "Output Options")]
    #[arg(short, long, default_value_t=String::from("10,50,90"))]
    percentiles: String,
    /// Quantiles to calculate, a space separated list of floats (for example:
    /// "0.1,0.2,0.9").
    #[clap(help_heading = "Output Options")]
    #[arg(short, long, conflicts_with = "percentiles", hide_short_help = true)]
    quantiles: Option<String>,
    /// Tabulate base modification probabilities at motifs (or single bases).
    /// The first argument should be the sequence motif and the second
    /// argument is the 0-based offset to the base to pileup base
    /// modification counts for. For example: --motif CGCG 0 indicates to
    /// collect base mofification probabilities for the first C on the top
    /// strand and the last C (complement to G) on the bottom strand. For
    /// single bases, --motif C 0 or --motif A 0 would restrict to
    /// probabilities on C or A reference bases, respectively.
    /// This argument can be passed multiple times.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, action = clap::ArgAction::Append, num_args = 2, requires = "reference_fasta")]
    motif: Option<Vec<String>>,
    /// Directory to deposit result tables into. Required for model probability
    /// histogram output.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'o', long)]
    out_dir: Option<PathBuf>,
    /// Label to prefix output files with.
    #[clap(help_heading = "Output Options")]
    #[arg(long, requires = "out_dir")]
    prefix: Option<String>,
    /// Overwrite results if present.
    #[clap(help_heading = "Output Options")]
    #[arg(long, requires = "out_dir", default_value_t = false)]
    force: bool,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    edge_filter: Option<String>,

    // probability histogram options
    /// Output histogram of base modification prediction probabilities.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "hist", requires = "out_dir", default_value_t = false)]
    histogram: bool,
    /// Set colors of primary bases in histogram, should be RGB format, e.g.
    /// "#0000FF" is defailt for canonical cytosine
    #[clap(help_heading = "Output Options")]
    #[arg(long="dna-color", requires = "histogram", num_args = 2, action = clap::ArgAction::Append)]
    primary_base_colors: Option<Vec<String>>,
    /// Set colors of modified bases in histogram, should be RGB format, e.g.
    /// "#FF00FF" is default for 5hmC
    #[clap(help_heading = "Output Options")]
    #[arg(long="mod-color", requires = "histogram", num_args = 2, action = clap::ArgAction::Append)]
    mod_base_colors: Option<Vec<String>>,

    // Approximate maximum number of reads to use, especially recommended when
    // using a large BAM without an index. If an indexed BAM is provided, the
    // reads will be sampled evenly over the length of the aligned reference.
    // If a region is passed with the --region option, they will be sampled
    // over the genomic region. Actual number of reads used may deviate
    // slightly from this number when using an indexed BAM.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        conflicts_with = "no_sampling"
    )]
    num_reads: Option<usize>,
    /// Instead of using a defined number of reads, specify a fraction of reads
    /// to sample, for example 0.1 will sample 1/10th of the reads.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        short = 'f',
        long,
        alias = "sample-frac",
        group = "sampling_options",
        conflicts_with = "no_sampling"
    )]
    sampling_frac: Option<f64>,
    /// No sampling, use all of the reads to calculate the filter thresholds.
    #[clap(help_heading = "Sampling Options")]
    #[arg(long, alias = "no_filtering", default_value_t = false)]
    no_sampling: bool,
    /// Random seed for deterministic running, the default is
    /// non-deterministic, only used when no BAM index is provided.
    #[clap(help_heading = "Sampling Options")]
    #[arg(short, conflicts_with = "no_sampling", long)]
    seed: Option<u64>,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    /// Only used when sampling probs from an indexed bam.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
    /// Only sample base modification probabilities that are aligned
    /// to the positions in this BED file. (alias: include-positions)
    #[clap(help_heading = "Selection Options")]
    #[arg(long, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Use non-primary mappings, may cause double-counting of reads with
    /// primary and one or more non-primary mappings. Requires MN tag to be
    /// correct. See documentation for help on proper alignment.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, default_value_t = false)]
    use_non_primary: bool,
}

impl SampleModBaseProbs {
    pub(crate) fn calc_counts_per_chrom(
        interval_size: u32,
        header: &HeaderView,
        num_reads: usize,
        region: Option<&Region>,
    ) -> anyhow::Result<FxHashMap<u32, usize>> {
        let total_reference_length = if let Some(region) = region {
            header
                .tid(region.name.as_bytes())
                .and_then(|tid| header.target_len(tid))
                .ok_or(anyhow!("didn't find {region} in header"))?
        } else {
            (0..header.target_count())
                .map(|tid| header.target_len(tid).unwrap_or(0))
                .sum::<u64>()
        };
        let out = (0..header.target_count())
            .filter(|tid| {
                region
                    .map(|region| {
                        header.tid2name(*tid) == region.name.as_bytes()
                    })
                    .unwrap_or(true)
            })
            .filter_map(|tid| header.target_len(tid).map(|l| (tid, l)))
            .map(|(tid, l)| {
                let frac = l as f64 / total_reference_length as f64;
                // debug!("tid {tid} is {frac} of the total length");
                let n_reads_for_chrom = frac * num_reads as f64;
                // debug!("sampling {n_reads_for_chrom} reads from {tid}");
                let frac_of_chrom_per_interval = if interval_size as u64 > l {
                    1f64
                } else {
                    interval_size as f64 / l as f64
                };
                // debug!(
                //     "each {interval_size} bp interval is \
                //      {frac_of_chrom_per_interval:.3} fraction of {tid} which
                // \      has lengh {l}"
                // );
                let n_reads_per_interval =
                    frac_of_chrom_per_interval * n_reads_for_chrom;
                let fudge_factor = n_reads_per_interval * 0.25;
                let total_reads_per_interval =
                    (n_reads_per_interval + fudge_factor).ceil() as usize;
                // debug!(
                //     "sampling {total_reads_per_interval} total reads per \
                //      interval from {tid} with fudge factor of {fudge_factor}
                // \      ({n_reads_per_interval} calculated)"
                // );
                (tid, total_reads_per_interval)
            })
            .collect();
        Ok(out)
    }

    fn get_region(
        &self,
        header: &HeaderView,
    ) -> anyhow::Result<Option<Region>> {
        if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Ok(Some(Region::parse_str(raw_region, header)?))
        } else {
            Ok(None)
        }
    }

    fn get_position_filter(
        &self,
        reference_records: &[ReferenceRecord],
    ) -> anyhow::Result<Option<StrandedPositionFilter<()>>> {
        self.include_bed
            .as_ref()
            .map(|bed_fp| {
                let chrom_to_tid = reference_records
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
            .transpose()
    }

    fn get_base_mod_probs_from_stream(
        &self,
        mut reader: bam::Reader,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<QualHist> {
        if self.region.is_some() {
            bail!("cannot use --region without indexed, sorted modBAM")
        }
        reader.set_threads(self.io_threads as usize)?;
        let header = reader.header();
        let reference_records = get_targets(header, None);

        let stranded_position_filter =
            self.get_position_filter(&reference_records)?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, false))
            .transpose()?;
        QualHist::from_records(
            reader.records(),
            stranded_position_filter,
            self.num_reads,
            self.sampling_frac,
            self.seed,
            edge_filter.as_ref(),
            self.histogram,
            false,
            &multi_progress,
        )
    }

    fn get_base_mod_probs_from_indexed_hts_file(
        &self,
        reader: bam::IndexedReader,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<QualHist> {
        let bam_fp = Path::new(&self.in_bam).to_path_buf();
        let num_reads = self.num_reads.unwrap_or(DEFAULT_NUM_READS);
        let thread_pool = rust_htslib::tpool::ThreadPool::new(self.io_threads)?;
        let header = reader.header();
        let region = self.get_region(header)?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, false))
            .transpose()?;
        let reference_records = get_targets(header, region.as_ref());
        let stranded_position_filter =
            self.get_position_filter(&reference_records)?;
        // should also adjust start/ends here.
        let reference_records = filter_reference_records(
            reference_records,
            &bam_fp,
            region.as_ref(),
            stranded_position_filter.as_ref(),
            &multi_progress,
        )?;

        let regex_motifs = if let Some(raw_motifs) = self.motif.as_ref() {
            Some(parse_raw_motifs(raw_motifs, false, true)?)
        } else {
            None
        };
        let motif_primary_bases = regex_motifs.as_ref().map(|rms| {
            rms.iter()
                .map(|x| x.motif_info.primary_base)
                .collect::<BTreeSet<DnaBase>>()
        });

        let motif_lookup = get_motif_lookup_from_parts(
            regex_motifs,
            self.reference_fasta.as_ref(),
            false,
            self.mask,
            self.preload_references,
        )?;

        let (chrom_to_counts, rng_sample, sample_frac) =
            if let Some(user_sampling_frac) = self.sampling_frac {
                multi_progress.suspend(|| {
                    info!("sampling {}% of reads", user_sampling_frac * 100f64);
                });
                (None, true, user_sampling_frac)
            } else {
                multi_progress.suspend(|| {
                    info!("attempting to sample at least {num_reads} reads");
                });
                (
                    Some(Arc::new(Self::calc_counts_per_chrom(
                        self.interval_size,
                        header,
                        num_reads,
                        region.as_ref(),
                    )?)),
                    false,
                    0f64,
                )
            };

        let mut workers: Vec<Box<dyn ExtractProbsWorker>> =
            Vec::with_capacity(self.threads);
        let mut motif_bases = [DnaBase::A; 4];
        if let Some(motif_primary_bases) = motif_primary_bases.as_ref() {
            assert!(motif_lookup.is_some());
            for (i, b) in motif_primary_bases.into_iter().enumerate() {
                motif_bases[i] = *b;
            }
        }

        if stranded_position_filter.is_some() || motif_lookup.is_some() {
            if let Some(motif_primary_bases) = motif_primary_bases.as_ref() {
                info!(
                    "scanning over motif bases: {}",
                    motif_primary_bases.iter().join(",")
                );
            }
            if stranded_position_filter.is_some() {
                info!("subsetting to posiitons in BED file");
            }
            if self.histogram {
                info!(
                    "collecting base and modification histograms at aligned \
                     positions"
                );
                for i in 0..self.threads {
                    let w = RegionMleProbs::<
                        AlignedBaseAndModArgmaxProbs,
                        ProbsExtractor,
                    >::new(
                        &bam_fp,
                        self.reference_fasta.as_ref(),
                        self.use_non_primary,
                        motif_bases,
                        edge_filter.as_ref(),
                        &thread_pool,
                        self.seed.unwrap_or(i as u64),
                        rng_sample,
                        sample_frac,
                        chrom_to_counts.clone(),
                    )?;
                    workers.push(Box::new(w));
                }
            } else {
                info!("collecting base-level histograms at aligned positions");
                for i in 0..self.threads {
                    let w = RegionMleProbs::<
                        AlignedBaseArgmaxProbs,
                        ProbsExtractor,
                    >::new(
                        &bam_fp,
                        self.reference_fasta.as_ref(),
                        self.use_non_primary,
                        motif_bases,
                        edge_filter.as_ref(),
                        &thread_pool,
                        self.seed.unwrap_or(i as u64),
                        rng_sample,
                        sample_frac,
                        chrom_to_counts.clone(),
                    )?;
                    workers.push(Box::new(w));
                }
            };
        } else {
            if self.histogram {
                info!(
                    "collecting base and modification histograms, using all \
                     read positions"
                );
                for i in 0..self.threads {
                    let w = RegionMleProbs::<
                        BaseAndModArgmaxProbs,
                        ProbsExtractor,
                    >::new(
                        &bam_fp,
                        self.reference_fasta.as_ref(),
                        self.use_non_primary,
                        motif_bases,
                        edge_filter.as_ref(),
                        &thread_pool,
                        self.seed.unwrap_or(i as u64),
                        rng_sample,
                        sample_frac,
                        chrom_to_counts.clone(),
                    )?;
                    workers.push(Box::new(w));
                }
            } else {
                info!(
                    "collecting base-level histograms, using all read \
                     positions"
                );
                for i in 0..self.threads {
                    let w =
                        RegionMleProbs::<BaseArgmaxProbs, ProbsExtractor>::new(
                            &bam_fp,
                            self.reference_fasta.as_ref(),
                            self.use_non_primary,
                            motif_bases,
                            edge_filter.as_ref(),
                            &thread_pool,
                            self.seed.unwrap_or(i as u64),
                            rng_sample,
                            sample_frac,
                            chrom_to_counts.clone(),
                        )?;
                    workers.push(Box::new(w));
                }
            };
        }
        let feeder = ChromCoordinatesFeeder::new(
            &reference_records,
            self.interval_size,
            motif_lookup,
            false,
            stranded_position_filter,
        )?;
        let qual_hist =
            run_extract_probs_workers(workers, multi_progress.clone(), feeder)?;

        if qual_hist.ok_records < num_reads && !self.no_sampling {
            multi_progress.suspend(|| {
                warn!(
                    "Only sampled {} reads when {} were requested. Has the \
                     input been subsampled already? Consider using \
                     --sample-frac or --no-sampling",
                    qual_hist.ok_records, num_reads
                );
            })
        } else {
            multi_progress.suspend(|| {
                info!("Sampled {} reads", qual_hist.ok_records);
            })
        }
        Ok(qual_hist)
    }

    fn get_base_mod_probs(
        &self,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<QualHist> {
        if using_stream(&self.in_bam) {
            info!(
                "reading BAM from stdin, --threads greater than 1 will not \
                 have an effect (--io-threads will)"
            );
            let reader = bam::Reader::from_stdin()?;
            self.get_base_mod_probs_from_stream(reader, multi_progress.clone())
        } else {
            if self.ignore_index {
                info!("running streaming scan (not splitting into chunks)");
                let reader = bam::Reader::from_path(&self.in_bam)?;
                return self.get_base_mod_probs_from_stream(
                    reader,
                    multi_progress.clone(),
                );
            }
            let Ok(indexed_reader) =
                bam::IndexedReader::from_path(&self.in_bam)
            else {
                let reader = bam::Reader::from_path(&self.in_bam)?;
                return self.get_base_mod_probs_from_stream(
                    reader,
                    multi_progress.clone(),
                );
            };
            info!(
                "collecting modification call probabilities from indexed HTS \
                 file"
            );
            self.get_base_mod_probs_from_indexed_hts_file(
                indexed_reader,
                multi_progress,
            )
        }
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let multi_progress = MultiProgress::new();
        if self.suppress_progress {
            multi_progress
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        if let Some(p) = self.out_dir.as_ref() {
            SampledProbs::check_files(
                p,
                self.prefix.as_ref(),
                self.force,
                self.histogram,
            )?;
        }

        let extra_dna_colors = if let Some(raw_dna_colors) =
            self.primary_base_colors.as_ref()
        {
            if raw_dna_colors.len() % 2 != 0 {
                bail!("invalid number of arguments")
            }
            raw_dna_colors
                .chunks(2)
                .map(|ch| {
                    let dna_base = ch[0]
                        .parse::<char>()
                        .map_err(|e| anyhow!("DNA base should be a char, {e}"))
                        .and_then(|c| DnaBase::parse(c).map_err(|e| e.into()));
                    let color_code = ch[1].to_string();
                    dna_base.map(|b| (b, color_code))
                })
                .collect::<anyhow::Result<HashMap<DnaBase, String>>>()?
        } else {
            HashMap::new()
        };
        let extra_mod_colors =
            if let Some(raw_mod_colors) = self.mod_base_colors.as_ref() {
                if raw_mod_colors.len() % 2 != 0 {
                    bail!("invalid number of arguments")
                }
                raw_mod_colors
                    .chunks(2)
                    .map(|ch| {
                        let mod_code = ModCodeRepr::parse(&ch[0]);
                        let color_code = ch[1].to_string();
                        mod_code.map(|b| (b, color_code))
                    })
                    .collect::<anyhow::Result<HashMap<ModCodeRepr, String>>>()?
            } else {
                HashMap::new()
            };

        let desired_percentiles = self
            .quantiles
            .as_ref()
            .map(|raw_quantiles| {
                parse_percentiles(raw_quantiles).map(|x| {
                    x.into_iter().map(|x| x * 100f32).collect::<Vec<f32>>()
                })
            })
            .unwrap_or_else(|| parse_percentiles(&self.percentiles))?;

        let qual_hist = self.get_base_mod_probs(multi_progress.clone())?;
        let percentiles = qual_hist
            .base_level_percentiles(&desired_percentiles, &multi_progress);

        let histograms = if self.histogram {
            let mut prob_counts = HashMap::new();
            for (i, base_total) in qual_hist.base_totals.iter().enumerate() {
                if *base_total > 0 {
                    let primary_base = DnaBase::try_from(i).unwrap();
                    let counts = qual_hist.hist[i]
                        .iter()
                        .enumerate()
                        .filter(|(_, c)| **c > 0u64)
                        .fold(BTreeMap::new(), |mut agg, (i, c)| {
                            assert!(i <= u8::MAX as usize);
                            agg.insert(i as u8, *c);
                            agg
                        });
                    let base_and_state: BaseAndState =
                        (primary_base, BaseState::Canonical(primary_base));
                    prob_counts.insert(base_and_state, counts);
                }
            }

            for mod_hist in qual_hist.mods_hists.iter() {
                let counts = mod_hist
                    .hist
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| **c > 0u64)
                    .fold(BTreeMap::new(), |mut agg, (i, c)| {
                        assert!(i <= u8::MAX as usize);
                        agg.insert(i as u8, *c);
                        agg
                    });
                let primary_base = mod_hist.dna_base;
                let mod_code = mod_hist.mod_code;
                let base_and_state =
                    (primary_base, BaseState::Modified(mod_code));
                prob_counts.insert(base_and_state, counts);
            }
            Some(ProbHistogram::new(prob_counts))
        } else {
            None
        };

        let sampled_probs = SampledProbs::new(
            histograms,
            percentiles,
            self.prefix.clone(),
            extra_dna_colors,
            extra_mod_colors,
        );
        //
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
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct ModSummarize {
    /// Input modBam, can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    in_bam: String,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Output summary as a tab-separated variables stdout instead of a table.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "tsv", default_value_t = false)]
    tsv_format: bool,
    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

    // sampling options
    /// Approximate maximum number of reads to use, especially recommended when
    /// using a large BAM without an index. If an indexed BAM is provided, the
    /// reads will be sampled evenly over the length of the aligned reference.
    /// If a region is passed with the --region option, they will be sampled
    /// over the region. Actual number of reads used may deviate
    /// slightly from this number.
    #[clap(help_heading = "Sampling Options")]
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
    #[clap(help_heading = "Sampling Options")]
    #[arg(group = "sampling_options", short = 'f', long)]
    sampling_frac: Option<f64>,
    /// No sampling, use all the reads to calculate the filter thresholds
    /// and generating the summary.
    #[clap(help_heading = "Sampling Options")]
    #[arg(long, group = "sampling_options", default_value_t = false)]
    no_sampling: bool,
    /// Sets a random seed for deterministic running (when using
    /// --sample-frac), the default is non-deterministic, only used when no
    /// BAM index is provided.
    #[clap(help_heading = "Sampling Options")]
    #[arg(short, requires = "sampling_frac", long)]
    seed: Option<u64>,

    // threshold options
    /// Do not perform any filtering, include all base modification calls in
    /// the summary. See filtering.md for details on filtering.
    #[clap(help_heading = "Filtering Options")]
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence base modification calls.
    #[clap(help_heading = "Filtering Options")]
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
    #[clap(help_heading = "Filtering Options")]
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
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long = "mod-threshold",
    alias = "mod-thresholds",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, group = "combine_args", hide_short_help = true)]
    ignore: Option<String>,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,
    /// Only summarize base modification probabilities that are aligned
    /// to the positions in this BED file. (alias: include-positions)
    #[clap(help_heading = "Selection Options")]
    #[arg(long, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Only use base modification probabilities that are aligned (i.e. ignore
    /// soft-clipped, and inserted bases).
    #[clap(
        long = "mapped-only",
        alias = "only-mapped",
        help_heading = "Selection Options"
    )]
    #[arg(long, default_value_t = false)]
    mapped_only: bool,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,
    /// When using regions, interval chunk size in base pairs to process
    /// concurrently. Smaller interval chunk sizes will use less memory but
    /// incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'i', long, default_value_t = 1_000_000)]
    interval_size: u32,
}

impl ModSummarize {
    pub fn run(&self) -> anyhow::Result<()> {
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
                        self.mapped_only || position_filter.is_some(),
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
                    self.mapped_only || position_filter.is_some(),
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
    #[clap(help_heading = "Compute Options")]
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
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Output SAM format instead of BAM.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    output_sam: bool,
}

fn update_mod_tags(
    mut record: bam::Record,
    no_implicit_calls: bool,
    new_mode: SkipMode,
) -> MkResult<bam::Record> {
    let mod_base_info = ModBaseInfo::new_from_record(&record)?;
    let mm_style = mod_base_info.mm_style;
    let ml_style = mod_base_info.ml_style;

    let mut mm_agg = String::new();
    let mut ml_agg = Vec::new();

    let (converters, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();
    for (dna_base, strand, mut seq_pos_mod_probs) in mod_prob_iter {
        let converter = converters.get(&dna_base).unwrap();
        if no_implicit_calls && new_mode == SkipMode::Explicit {
            seq_pos_mod_probs = seq_pos_mod_probs.remove_implicit_probs();
        } else {
            seq_pos_mod_probs.set_skip_mode(new_mode);
        }
        let (mm, mut ml) = format_mm_ml_tag(
            seq_pos_mod_probs,
            dna_base,
            &converter.cumulative_counts,
            strand,
        );
        mm_agg.push_str(&mm);
        ml_agg.extend_from_slice(&mut ml);
    }
    record.remove_aux(mm_style.as_bytes())?;
    record.remove_aux(ml_style.as_bytes())?;
    let mm = Aux::String(&mm_agg);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml_agg;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record.push_aux(MM_TAGS[0].as_bytes(), mm)?;
    record.push_aux(ML_TAGS[0].as_bytes(), ml)?;

    Ok(record)
}

impl Update {
    pub fn run(&self) -> anyhow::Result<()> {
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
        let mut error_counts = FxHashMap::<String, usize>::default();

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
                SkipMode::ImplicitUnmodified
            }
        };

        let record_iter = reader
            .records()
            .map(|res| res.map_err(|e| MkError::HtsLibError(e)))
            .map(|res| {
                res.and_then(|record| {
                    update_mod_tags(record, self.no_implicit_probs, to_mode)
                })
            });

        for (i, result) in record_iter.enumerate() {
            match result {
                Ok(record) => {
                    match out_bam
                        .write(&record)
                        .map_err(|e| MkError::HtsLibError(e))
                    {
                        Ok(_) => {
                            spinner.inc(1);
                            total = i + 1;
                        }
                        Err(mk_error) => {
                            error_counts
                                .entry(mk_error.to_string())
                                .or_insert(0usize)
                                .add_assign(1usize);
                        }
                    }
                }
                Err(mk_error) => {
                    error_counts
                        .entry(mk_error.to_string())
                        .or_insert(0usize)
                        .add_assign(1usize);
                }
            }
        }

        spinner.finish_and_clear();

        info!("done, {} records processed", total);

        if !error_counts.is_empty() {
            info!("error/skip counts:");
            let error_table = format_errors_table(&error_counts);
            info!("\n{error_table}");
        }

        Ok(())
    }
}

// todo this command is missing argument groups! It should/could also have
// --ignore
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
    #[arg(long, alias = "log")]
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
    long,
    long = "mod-threshold",
    alias = "mod-thresholds",
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

    /// Filter out any base modification call that isn't part of a basecall
    /// sequence motif This argument can be passed multiple times. Format
    /// is <motif_sequence> <offset>. For example the argument to match CpG
    /// dinucleotides is `--motif CG 0`, or to match CG[5mC]G the argument
    /// would be `--motif CGCG 2`.
    #[arg(long, action = clap::ArgAction::Append, num_args = 2)]
    motif: Option<Vec<String>>,
    /// Shorthand for --motif CG 0.
    #[arg(long, default_value_t = false)]
    cpg: bool,
    /// Discard base modification calls that match the provided motifs (instead
    /// of keeping them).
    #[arg(long, requires = "motif", default_value_t = false)]
    discard_motifs: bool,

    /// Output SAM format instead of BAM.
    #[arg(long, default_value_t = false)]
    output_sam: bool,
}

impl CallMods {
    pub fn run(&self) -> anyhow::Result<()> {
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

        let motifs = parse_forward_motifs(&self.motif, self.cpg)?;
        if let Some(ms) = motifs.as_ref() {
            let patterns = ms.iter().map(|x| x.as_str()).join(",");
            info!("filtering base modification calls to patterns: {patterns}");
        }

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
            &motifs,
            self.discard_motifs,
            "Calling Mods, records processed",
            self.suppress_progress,
            false,
        )?;

        Ok(())
    }
}

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
            Self::implicit => SkipMode::ImplicitUnmodified,
        }
    }
}
