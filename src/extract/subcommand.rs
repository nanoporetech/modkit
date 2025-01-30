use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::bail;
use bio::io::fasta::Reader as FastaReader;
use clap::{Args, Subcommand};
use crossbeam_channel::bounded;
use indicatif::{MultiProgress, ProgressIterator};
use log::{debug, error, info};
use rayon::{ThreadPool, ThreadPoolBuilder};
use rust_htslib::bam::{self, Read};

use crate::command_utils::{
    get_serial_reader, get_threshold_from_options, parse_edge_filter_input,
    parse_per_mod_thresholds, parse_thresholds, using_stream,
};
use crate::extract::args::InputArgs;
use crate::extract::util::ReferencePositionFilter;
use crate::extract::writer::{OutwriterWithMemory, TsvWriterWithContigNames};
use crate::interval_chunks::ReferenceIntervalsFeeder;
use crate::logging::init_logging_smart;
use crate::mod_bam::CollapseMethod;
use crate::mod_base_code::ModCodeRepr;
use crate::read_ids_to_base_mod_probs::{
    ModProfile, PositionModCalls, ReadsBaseModProfile,
};
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::WithRecords;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{get_ticker, Region, KMER_SIZE};
use crate::writers::TsvWriter;

#[derive(Subcommand)]
pub enum ExtractMods {
    /// Transform the probabilities from the MM/ML tags in a modBAM into a
    /// table.
    Full(EntryExtractFull),
    /// Produce a table of read-level base modification calls. This table has,
    /// for each read, one row for each base modification call in that read
    /// using the same thresholding algorithm as in pileup, or summary (see
    /// online documentation for details on thresholds).
    Calls(EntryExtractCalls),
}

impl ExtractMods {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            ExtractMods::Full(x) => x.run(),
            ExtractMods::Calls(x) => x.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryExtractFull {
    #[clap(flatten)]
    input_args: InputArgs,
    /// Path to reference FASTA to extract reference context information from.
    /// Required for motif selection.
    #[arg(long, alias = "ref")]
    pub reference: Option<PathBuf>,
}

impl EntryExtractFull {
    fn using_stdin(&self) -> bool {
        using_stream(&self.input_args.in_bam)
    }

    fn load_regions(
        &self,
        name_to_tid: &HashMap<&str, u32>,
        region: Option<&Region>,
        contigs: &HashMap<String, Vec<u8>>,
        master_progress_bar: &MultiProgress,
        thread_pool: &ThreadPool,
    ) -> anyhow::Result<(
        Option<ReferenceIntervalsFeeder>,
        ReferencePositionFilter,
    )> {
        super::util::load_regions(
            &self.input_args,
            self.using_stdin(),
            name_to_tid,
            region,
            contigs,
            master_progress_bar,
            thread_pool,
        )
    }

    pub(crate) fn run(&self) -> anyhow::Result<()> {
        let stream_out = using_stream(self.input_args.out_path.as_str());
        let _handle = init_logging_smart(
            self.input_args.log_filepath.as_ref(),
            stream_out,
        );
        if self.input_args.out_threads == 0 {
            bail!("output threads must be >= 1")
        }

        if self.input_args.kmer_size > KMER_SIZE {
            bail!("kmer size must be less than or equal to {KMER_SIZE}")
        }

        let multi_prog = MultiProgress::new();
        if self.input_args.suppress_progress || stream_out {
            multi_prog.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let pool = ThreadPoolBuilder::new()
            .num_threads(self.input_args.threads)
            .build()?;

        let collapse_method = match &self.input_args.ignore {
            Some(raw_mod_code) => {
                let mod_code = ModCodeRepr::parse(raw_mod_code)?;
                Some(CollapseMethod::ReDistribute(mod_code))
            }
            None => None,
        };
        let edge_filter = self
            .input_args
            .edge_filter
            .as_ref()
            .map(|raw| {
                parse_edge_filter_input(raw, self.input_args.invert_edge_filter)
            })
            .transpose()?;

        let mut reader = get_serial_reader(&self.input_args.in_bam)?;
        let header = reader.header().to_owned();

        let queue_size = self.input_args.queue_size;
        let (snd, rcv) = bounded(queue_size);

        let tid_to_name = (0..header.target_count())
            .filter_map(|tid| {
                match String::from_utf8(header.tid2name(tid).to_vec()) {
                    Ok(contig) => Some((tid, contig)),
                    Err(e) => {
                        error!(
                            "failed to parse contig {tid}, {}",
                            e.to_string()
                        );
                        None
                    }
                }
            })
            .collect::<HashMap<u32, String>>();
        let name_to_tid = tid_to_name
            .iter()
            .map(|(tid, name)| (name.as_str(), *tid))
            .collect::<HashMap<&str, u32>>();

        let chrom_to_seq = match self.reference.as_ref() {
            Some(fp) => {
                let reader = FastaReader::from_file(fp)?;
                let pb = multi_prog.add(get_ticker());
                pb.set_message("parsing FASTA records");
                reader
                    .records()
                    .progress_with(pb)
                    .filter_map(|r| r.ok())
                    .filter(|record| name_to_tid.get(record.id()).is_some())
                    .map(|record| {
                        (record.id().to_owned(), record.seq().to_vec())
                    })
                    .collect::<HashMap<String, Vec<u8>>>()
            }
            None => HashMap::new(),
        };

        let region = self
            .input_args
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, &header))
            .transpose()?;

        let (references_and_intervals, reference_position_filter) = self
            .load_regions(
                &name_to_tid,
                region.as_ref(),
                &chrom_to_seq,
                &multi_prog,
                &pool,
            )?;

        // allowed to use the sampling schedule if there is an index, if
        // asked for num_reads with no index, scan first N reads
        let schedule = match (self.input_args.num_reads, self.using_stdin()) {
            (_, true) | (None, false) => None,
            (Some(num_reads), false) => {
                match bam::IndexedReader::from_path(&self.input_args.in_bam) {
                    Ok(_) => Some(SamplingSchedule::from_num_reads(
                        &self.input_args.in_bam,
                        num_reads,
                        region.as_ref(),
                        reference_position_filter.include_pos.as_ref(),
                        reference_position_filter.include_unmapped_reads,
                    )?),
                    Err(_) => {
                        debug!(
                            "cannot use sampling schedule without index, \
                             keeping first {num_reads} reads"
                        );
                        None
                    }
                }
            }
        };

        let n_failed = multi_prog.add(get_ticker());
        n_failed.set_message("~records failed");
        let n_skipped = multi_prog.add(get_ticker());
        n_skipped.set_message("~records skipped");
        let n_used = multi_prog.add(get_ticker());
        n_used.set_message("~records used");
        let n_rows = multi_prog.add(get_ticker());
        n_rows.set_message("rows written");

        reader.set_threads(self.input_args.threads)?;
        let n_reads = self.input_args.num_reads;
        let threads = self.input_args.threads;
        let mapped_only = self.input_args.mapped_only;
        let in_bam = self.input_args.in_bam.clone();
        let kmer_size = self.input_args.kmer_size;
        let allow_non_primary = self.input_args.allow_non_primary;
        let remove_inferred = self.input_args.ignore_implicit;

        pool.spawn(move || {
            super::util::run_extract_reads(
                reader,
                in_bam,
                references_and_intervals,
                schedule,
                collapse_method,
                edge_filter,
                allow_non_primary,
                kmer_size,
                remove_inferred,
                reference_position_filter,
                snd,
                queue_size,
                n_reads,
                threads,
                mapped_only,
                multi_prog,
            );
        });

        let output_header = if self.input_args.no_headers {
            None
        } else {
            Some(ModProfile::header())
        };
        let mut writer: Box<dyn OutwriterWithMemory<ReadsBaseModProfile>> =
            match self.input_args.out_path.as_str() {
                "stdout" | "-" => {
                    let tsv_writer = TsvWriter::new_stdout(output_header);
                    let writer = TsvWriterWithContigNames::new(
                        tsv_writer,
                        tid_to_name,
                        chrom_to_seq,
                    )?;
                    Box::new(writer)
                }
                _ => {
                    if self.input_args.bgzf {
                        let tsv_writer = TsvWriter::new_gzip(
                            &self.input_args.out_path,
                            self.input_args.force,
                            self.input_args.out_threads,
                            output_header,
                        )?;
                        let writer = TsvWriterWithContigNames::new(
                            tsv_writer,
                            tid_to_name,
                            chrom_to_seq,
                        )?;
                        Box::new(writer)
                    } else {
                        let tsv_writer = TsvWriter::new_file(
                            &self.input_args.out_path,
                            self.input_args.force,
                            output_header,
                        )?;
                        let writer = TsvWriterWithContigNames::new(
                            tsv_writer,
                            tid_to_name,
                            chrom_to_seq,
                        )?;
                        Box::new(writer)
                    }
                }
            };

        for result in rcv {
            match result {
                Ok(mod_profile) => {
                    n_used.inc(mod_profile.num_reads() as u64);
                    n_failed.inc(mod_profile.num_fails as u64);
                    n_skipped.inc(mod_profile.num_skips as u64);
                    match writer.write(mod_profile, kmer_size) {
                        Ok(n) => n_rows.inc(n),
                        Err(e) => {
                            error!("failed to write {}", e.to_string());
                        }
                    }
                }
                Err(e) => {
                    debug!(
                        "failed to calculate read-level mod probs, {}",
                        e.to_string()
                    );
                }
            }
        }

        n_failed.finish_and_clear();
        n_skipped.finish_and_clear();
        n_used.finish_and_clear();
        n_rows.finish_and_clear();
        info!(
            "processed {} reads, {} rows, skipped ~{} reads, failed ~{} reads",
            writer.num_reads(),
            n_rows.position(),
            n_skipped.position(),
            n_failed.position()
        );
        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryExtractCalls {
    #[clap(flatten)]
    input_args: InputArgs,
    /// Path to reference FASTA to extract reference context information from.
    /// If no reference is provided, `ref_kmer` column will be "." in the
    /// output. (alias: ref)
    #[arg(long, alias = "ref")]
    pub reference: Option<PathBuf>,

    /// Only output base modification calls that pass the minimum confidence
    /// threshold. (alias: pass)
    #[clap(help_heading = "Selection Options")]
    #[arg(long, alias = "pass", default_value_t = false)]
    pass_only: bool,
    // sampling and filtering
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
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        long,
        alias = "mod-threshold",
        action = clap::ArgAction::Append,
        hide_short_help = true
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Don't estimate the pass threshold, all calls will "pass".
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        conflicts_with_all = ["mod_thresholds", "filter_threshold", "pass_only"],
        long,
        default_value_t = false,
        hide_short_help = true
    )]
    no_filtering: bool,
    /// Interval chunk size in base pairs to process concurrently when
    /// estimating the threshold probability.
    #[clap(help_heading = "Sampling Options")]
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,
    /// Sample this fraction of the reads when estimating the pass-threshold.
    /// In practice, 10-100 thousand reads is sufficient to estimate the model
    /// output distribution and determine the filtering threshold. See
    /// filtering.md for details on filtering.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Sample this many reads when estimating the filtering threshold. If a
    /// sorted, indexed modBAM is provided reads will be sampled evenly
    /// across aligned genome. If a region is specified, with the --region,
    /// then reads will be sampled evenly across the region given.
    /// This option is useful for large BAM files. In practice, 10-50 thousand
    /// reads is sufficient to estimate the model output distribution and
    /// determine the filtering threshold.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    sample_num_reads: usize,
    /// Set a random seed for deterministic running, the default is
    /// non-deterministic when using `sampling_frac`. When using `num_reads`
    /// the output is still deterministic.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
}

impl EntryExtractCalls {
    fn using_stdin(&self) -> bool {
        using_stream(&self.input_args.in_bam)
    }

    fn run(&self) -> anyhow::Result<()> {
        let stream_out = using_stream(self.input_args.out_path.as_str());

        let _handle = init_logging_smart(
            self.input_args.log_filepath.as_ref(),
            stream_out,
        );

        if self.input_args.out_threads == 0 {
            bail!("output threads must be >= 1")
        }

        if self.input_args.kmer_size > KMER_SIZE {
            bail!("kmer size must be less than or equal to {KMER_SIZE}")
        }

        let multi_prog = MultiProgress::new();
        if self.input_args.suppress_progress || stream_out {
            multi_prog.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let pool = ThreadPoolBuilder::new()
            .num_threads(self.input_args.threads)
            .build()?;

        let collapse_method = match &self.input_args.ignore {
            Some(raw_mod_code) => {
                let mod_code = ModCodeRepr::parse(raw_mod_code)?;
                Some(CollapseMethod::ReDistribute(mod_code))
            }
            None => None,
        };
        let edge_filter = self
            .input_args
            .edge_filter
            .as_ref()
            .map(|raw| {
                parse_edge_filter_input(raw, self.input_args.invert_edge_filter)
            })
            .transpose()?;

        let mut reader = get_serial_reader(&self.input_args.in_bam)?;
        let header = reader.header().to_owned();

        let tid_to_name = (0..header.target_count())
            .filter_map(|tid| {
                match String::from_utf8(header.tid2name(tid).to_vec()) {
                    Ok(contig) => Some((tid, contig)),
                    Err(e) => {
                        error!(
                            "failed to parse contig {tid}, {}",
                            e.to_string()
                        );
                        None
                    }
                }
            })
            .collect::<HashMap<u32, String>>();

        let name_to_tid = tid_to_name
            .iter()
            .map(|(tid, name)| (name.as_str(), *tid))
            .collect::<HashMap<&str, u32>>();

        let chrom_to_seq = match self.reference.as_ref() {
            Some(fp) => {
                let reader = FastaReader::from_file(fp)?;
                let pb = multi_prog.add(get_ticker());
                pb.set_message("parsing FASTA records");
                reader
                    .records()
                    .progress_with(pb)
                    .filter_map(|r| r.ok())
                    .filter(|record| name_to_tid.get(record.id()).is_some())
                    .map(|record| {
                        (record.id().to_owned(), record.seq().to_vec())
                    })
                    .collect::<HashMap<String, Vec<u8>>>()
            }
            None => HashMap::new(),
        };

        let region = self
            .input_args
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, &header))
            .transpose()?;

        let per_mod_thresholds = self
            .mod_thresholds
            .as_ref()
            .map(|raw_per_mod_thresholds| {
                parse_per_mod_thresholds(raw_per_mod_thresholds)
            })
            .transpose()?;

        let (references_and_intervals, reference_position_filter) =
            super::util::load_regions(
                &self.input_args,
                self.using_stdin(),
                &name_to_tid,
                region.as_ref(),
                &chrom_to_seq,
                &multi_prog,
                &pool,
            )?;

        let caller = if !self.no_filtering {
            // stdin input and want a threshold, not allowed
            if self.using_stdin() && self.filter_threshold.is_none() {
                bail!(
                    "\
                        cannot use stdin and estimate a filter threshold, set \
                     the threshold on the command line with \
                     --filter-threshold and/or --mod-threshold (or set \
                     --no-filtering)."
                )
            }
            if let Some(raw_threshold) = &self.filter_threshold {
                parse_thresholds(raw_threshold, per_mod_thresholds)?
            } else {
                let in_bam = Path::new(&self.input_args.in_bam).to_path_buf();
                if !in_bam.exists() {
                    bail!(
                        "failed to find input modBAM file at {}",
                        self.input_args.in_bam
                    );
                }
                pool.install(|| {
                    get_threshold_from_options(
                        &in_bam,
                        self.input_args.threads,
                        self.sampling_interval_size,
                        self.sampling_frac,
                        self.sample_num_reads,
                        false,
                        self.filter_percentile,
                        self.seed,
                        region.as_ref(),
                        per_mod_thresholds,
                        edge_filter.as_ref(),
                        collapse_method.as_ref(),
                        reference_position_filter.include_pos.as_ref(),
                        reference_position_filter.only_mapped_positions(),
                        self.input_args.suppress_progress,
                    )
                })?
            }
        } else {
            MultipleThresholdModCaller::new_passthrough()
        };

        let output_header = if self.input_args.no_headers {
            None
        } else {
            Some(PositionModCalls::header())
        };
        let mut writer: Box<dyn OutwriterWithMemory<ReadsBaseModProfile>> =
            match self.input_args.out_path.as_str() {
                "stdout" | "-" => {
                    let tsv_writer = TsvWriter::new_stdout(output_header);
                    let writer = TsvWriterWithContigNames::new_with_caller(
                        tsv_writer,
                        tid_to_name,
                        chrom_to_seq,
                        caller,
                        self.pass_only,
                    )?;
                    Box::new(writer)
                }
                _ => {
                    if self.input_args.bgzf {
                        let tsv_writer = TsvWriter::new_gzip(
                            &self.input_args.out_path,
                            self.input_args.force,
                            self.input_args.out_threads,
                            output_header,
                        )?;
                        let writer = TsvWriterWithContigNames::new_with_caller(
                            tsv_writer,
                            tid_to_name,
                            chrom_to_seq,
                            caller,
                            self.pass_only,
                        )?;
                        Box::new(writer)
                    } else {
                        let tsv_writer = TsvWriter::new_file(
                            &self.input_args.out_path,
                            self.input_args.force,
                            output_header,
                        )?;
                        let writer = TsvWriterWithContigNames::new_with_caller(
                            tsv_writer,
                            tid_to_name,
                            chrom_to_seq,
                            caller,
                            self.pass_only,
                        )?;
                        Box::new(writer)
                    }
                }
            };

        let schedule = match (self.input_args.num_reads, self.using_stdin()) {
            (_, true) | (None, false) => None,
            (Some(num_reads), false) => {
                match bam::IndexedReader::from_path(&self.input_args.in_bam) {
                    Ok(_) => Some(SamplingSchedule::from_num_reads(
                        &self.input_args.in_bam,
                        num_reads,
                        region.as_ref(),
                        reference_position_filter.include_pos.as_ref(),
                        reference_position_filter.include_unmapped_reads,
                    )?),
                    Err(_) => {
                        debug!(
                            "cannot use sampling schedule without index, \
                             keeping first {num_reads} reads"
                        );
                        None
                    }
                }
            }
        };

        let queue_size = self.input_args.queue_size;
        let (snd, rcv) = bounded(queue_size);

        let n_failed = multi_prog.add(get_ticker());
        n_failed.set_message("~records failed");
        let n_skipped = multi_prog.add(get_ticker());
        n_skipped.set_message("~records skipped");
        let n_used = multi_prog.add(get_ticker());
        n_used.set_message("~records used");
        let n_rows = multi_prog.add(get_ticker());
        n_rows.set_message("rows written");

        reader.set_threads(self.input_args.threads)?;
        let n_reads = self.input_args.num_reads;
        let threads = self.input_args.threads;
        let mapped_only = self.input_args.mapped_only;
        let in_bam = self.input_args.in_bam.clone();
        let kmer_size = self.input_args.kmer_size;
        let allow_non_primary = self.input_args.allow_non_primary;
        let remove_inferred = self.input_args.ignore_implicit;

        pool.spawn(move || {
            super::util::run_extract_reads(
                reader,
                in_bam,
                references_and_intervals,
                schedule,
                collapse_method,
                edge_filter,
                allow_non_primary,
                kmer_size,
                remove_inferred,
                reference_position_filter,
                snd,
                queue_size,
                n_reads,
                threads,
                mapped_only,
                multi_prog,
            );
        });

        for result in rcv {
            match result {
                Ok(mod_profile) => {
                    n_used.inc(mod_profile.num_reads() as u64);
                    n_failed.inc(mod_profile.num_fails as u64);
                    n_skipped.inc(mod_profile.num_skips as u64);
                    match writer.write(mod_profile, 0) {
                        Ok(n) => n_rows.inc(n),
                        Err(e) => {
                            error!("failed to write {}", e.to_string());
                        }
                    }
                }
                Err(e) => {
                    debug!(
                        "failed to calculate read-level mod probs, {}",
                        e.to_string()
                    );
                }
            }
        }

        n_failed.finish_and_clear();
        n_skipped.finish_and_clear();
        n_used.finish_and_clear();
        n_rows.finish_and_clear();
        info!(
            "processed {} reads, {} rows, skipped ~{} reads, failed ~{} reads",
            writer.num_reads(),
            n_rows.position(),
            n_skipped.position(),
            n_failed.position()
        );
        Ok(())
    }
}
