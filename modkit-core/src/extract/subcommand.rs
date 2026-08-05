use itertools::Itertools;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context};
use clap::{Args, Subcommand};
use common_macros::hash_map;
use crossbeam_channel::bounded;
use indicatif::{MultiProgress, ProgressIterator};
use log::{debug, error, info};
use rayon::{ThreadPool, ThreadPoolBuilder};
use rust_htslib::bam::{self, Read};

use modkit_logging::{init_logging, init_logging_smart};

use crate::command_utils::{
    get_serial_reader, parse_edge_filter_input, parse_per_base_thresholds,
    parse_per_mod_thresholds, parse_raw_thresholds_string_with_default,
    parse_sampling_fraction, parse_thresholds, using_stream,
};
use crate::extract::args::InputArgs;
use crate::extract::util::{
    calc_per_base_thresholds_from_stream, ReadModStatsProcessor,
    ReadModsStatsRecord, ReferencePositionFilter,
};
use crate::extract::writer::{
    CanWriteReadModStatsRecords, OutwriterWithMemory, ReadModStatsWriter,
    TsvWriterWithContigNames,
};
use crate::fasta::HtsLibFastaRecords;
use crate::interval_chunks::ReferenceIntervalBatchesFeeder;
use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCodeRepr, ModifiedBasesOptions};
use crate::motifs::motif_bed::MotifPositionLookup;
use crate::position_filter::StrandedPositionFilter;
use crate::read_ids_to_base_mod_probs::{
    ModProfile, PositionModCalls, ReadsBaseModProfile,
};
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::WithRecords;
use crate::sample_probs::calc_per_base_thresholds_from_indexed_hts_file;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{format_errors_table, get_ticker, Region, KMER_SIZE};
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
    /// Produce a table where modification counts are summarized on the read
    /// level. This table will have one record per valid read and count the
    /// number of modified and unmodified bases for each base modification
    /// requested.
    ReadStats(EntryReadStats),
}

impl ExtractMods {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            ExtractMods::Full(x) => x.run(),
            ExtractMods::Calls(x) => x.run(),
            ExtractMods::ReadStats(x) => x.run(),
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
        Option<ReferenceIntervalBatchesFeeder>,
        ReferencePositionFilter,
        Option<MotifPositionLookup>,
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
                let reader = HtsLibFastaRecords::from_file(fp)?;
                let pb = multi_prog.add(get_ticker());
                pb.set_message("parsing FASTA records");
                reader
                    .into_iter()
                    .progress_with(pb)
                    .filter_map(|r| r.ok())
                    .filter(|(record_id, _record)| {
                        name_to_tid.get(record_id.as_str()).is_some()
                    })
                    .map(|(record_id, record)| {
                        (record_id, record.as_bytes().to_vec())
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

        let (
            references_and_intervals,
            reference_position_filter,
            motif_position_lookup,
        ) = self.load_regions(
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

        let with_motifs = self.input_args.motif.is_some();
        let output_header = if self.input_args.no_headers {
            None
        } else {
            Some(ModProfile::header(with_motifs))
        };
        let mut writer: Box<dyn OutwriterWithMemory<ReadsBaseModProfile>> =
            match self.input_args.out_path.as_str() {
                "stdout" | "-" => {
                    let tsv_writer = TsvWriter::new_stdout(output_header);
                    let writer = TsvWriterWithContigNames::new(
                        tsv_writer,
                        tid_to_name,
                        chrom_to_seq,
                        with_motifs,
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
                            with_motifs,
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
                            with_motifs,
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
                    match writer
                        .write(mod_profile, motif_position_lookup.as_ref())
                    {
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
    /// Preload the reference sequences, useful when working with many, short
    /// reference sequences such as a transcriptome.
    #[arg(long, requires = "reference", default_value_t = false)]
    preload_references: bool,

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
        long = "mod-threshold",
        alias = "mod-thresholds",
        action = clap::ArgAction::Append,
        hide_short_help = true
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Maximum threhsold value to use, if the estimated threshold value for
    /// any primary sequence base is greater than this value, use this value
    /// instead. To set a default value for all primary bases, use
    /// --max-threshold 0.8, for example. Or use per-base threshold values like
    /// --max-threshold C:0.8.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long,
    group = "thresholds",
    alias = "max-threshold",
    action = clap::ArgAction::Append,
    )]
    max_filter_threshold: Option<Vec<String>>,
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
    /// Must be a finite value in the inclusive range [0, 1].
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true,
        value_parser = parse_sampling_fraction
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
    fn calc_per_base_thersholds(
        &self,
        region: Option<&Region>,
        per_mod_thresholds: Option<&HashMap<ModCodeRepr, f32>>,
        edge_filter: Option<&EdgeFilter>,
        position_filter: Option<StrandedPositionFilter<()>>,
        io_threadpool: &rust_htslib::tpool::ThreadPool,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<MultipleThresholdModCaller> {
        let max_thresholds_per_base = self
            .max_filter_threshold
            .as_ref()
            .map(|raw_thresholds| {
                parse_raw_thresholds_string_with_default(raw_thresholds)
            })
            .transpose()?;
        let per_base_thresholds = if bam::IndexedReader::from_path(
            &self.input_args.in_bam,
        )
        .is_err()
        {
            calc_per_base_thresholds_from_stream(
                &Path::new(&self.input_args.in_bam).to_path_buf(),
                self.sample_num_reads,
                false,
                position_filter,
                edge_filter,
                self.filter_percentile,
                self.input_args.io_threads as usize,
                max_thresholds_per_base,
                &multi_progress,
            )?
        } else {
            calc_per_base_thresholds_from_indexed_hts_file(
                &Path::new(&self.input_args.in_bam).to_path_buf(),
                self.reference.as_ref(),
                self.filter_percentile,
                false,
                self.preload_references,
                self.input_args.mask,
                self.sampling_frac,
                self.sample_num_reads,
                max_thresholds_per_base,
                position_filter,
                self.input_args.motif.as_ref(),
                self.input_args.cpg,
                self.input_args.threads,
                self.sampling_interval_size,
                self.seed,
                region,
                edge_filter,
                io_threadpool,
                multi_progress,
            )?
        };
        let per_base_thresholds = hash_map! {
            DnaBase::A => per_base_thresholds[DnaBase::A as usize],
            DnaBase::C => per_base_thresholds[DnaBase::C as usize],
            DnaBase::G => per_base_thresholds[DnaBase::G as usize],
            DnaBase::T => per_base_thresholds[DnaBase::T as usize],
        };
        Ok(MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds.cloned().unwrap_or_else(HashMap::new),
            0f32,
        ))
    }

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
                let reader = HtsLibFastaRecords::from_file(fp)?;
                let pb = multi_prog.add(get_ticker());
                pb.set_message("parsing FASTA records");
                reader
                    .into_iter()
                    .progress_with(pb)
                    .filter_map(|r| r.ok())
                    .filter(|(read_id, _record)| {
                        name_to_tid.get(read_id.as_str()).is_some()
                    })
                    .map(|(record_id, record)| {
                        (record_id, record.as_bytes().to_vec())
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

        let (
            references_and_intervals,
            reference_position_filter,
            motif_position_lookup,
        ) = super::util::load_regions(
            &self.input_args,
            self.using_stdin(),
            &name_to_tid,
            region.as_ref(),
            &chrom_to_seq,
            &multi_prog,
            &pool,
        )?;

        let io_threadpool =
            rust_htslib::tpool::ThreadPool::new(self.input_args.io_threads)?;

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
                self.calc_per_base_thersholds(
                    region.as_ref(),
                    per_mod_thresholds.as_ref(),
                    edge_filter.as_ref(),
                    reference_position_filter.include_pos.clone(),
                    &io_threadpool,
                    multi_prog.clone(),
                )?
            }
        } else {
            MultipleThresholdModCaller::new_passthrough()
        };
        // TODO(arand) once I refactor extract, I'll want to keep this around.
        drop(io_threadpool);

        let with_motifs = self.input_args.motif.is_some();
        let output_header = if self.input_args.no_headers {
            None
        } else {
            Some(PositionModCalls::header(with_motifs))
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
                        with_motifs,
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
                            with_motifs,
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
                            with_motifs,
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
                    match writer
                        .write(mod_profile, motif_position_lookup.as_ref())
                    {
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
pub struct EntryReadStats {
    /// Input modBAM. Can be a path to a file or "-"/"stdin" for standard input
    in_bam: String,
    /// Path to file for output, default will be to stdout
    out_path: String,
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
    #[arg(long)]
    filter_threshold: Vec<String>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    #[arg(long)]
    mod_thresholds: Option<Vec<String>>,
    /// Specify which modcodes to tabulate, should be <primary_base>:<mod_code>
    #[arg(
        long,
        alias = "mod-code",
        value_parser = clap::value_parser!(ModifiedBasesOptions),
        num_args = 1..,
        value_delimiter = ' '
    )]
    mod_codes: Vec<ModifiedBasesOptions>,
    /// Path to file to write run log.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Queue size, maximum number of reads to be held in memory waiting for
    /// records to be written
    #[arg(long, default_value_t = 10_000usize)]
    queue_size: usize,
    #[arg(long, default_value_t = 4usize)]
    io_threads: usize,
    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    pub suppress_progress: bool,
}

impl EntryReadStats {
    fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());
        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        if self.mod_codes.is_empty() {
            bail!("need to provide at least 1 base:mod-code pair")
        }
        let per_base_thresholds = self.parse_thresholds()?;
        let per_mod_thresholds = self.parse_per_mod_thresholds()?;
        let mut header = vec![
            "read_id".to_string(),
            "chrom".to_string(),
            "aln_start".to_string(),
        ];
        let bases_to_codes = self
            .mod_codes
            .iter()
            .map(|x| (x.primary_base, x.mod_code))
            .collect::<Vec<(DnaBase, ModCodeRepr)>>();

        let grouped = bases_to_codes
            .iter()
            .fold(HashMap::new(), |mut agg, (base, code)| {
                agg.entry(*base).or_insert(HashSet::new()).insert(*code);
                agg
            })
            .into_iter()
            .map(|(b, codes)| {
                (b, codes.into_iter().sorted().collect::<Vec<_>>())
            })
            .collect::<HashMap<DnaBase, Vec<ModCodeRepr>>>();

        let mut reader = match self.in_bam.as_str() {
            "-" | "stdin" => bam::Reader::from_stdin()?,
            p @ _ => {
                let fp = Path::new(&self.in_bam);
                if !fp.exists() {
                    bail!("didn't find BAM at {p}");
                }
                bam::Reader::from_path(fp)?
            }
        };
        let contigs = reader
            .header()
            .target_names()
            .iter()
            .map(|raw_name| {
                raw_name.iter().map(|b| *b as char).collect::<String>()
            })
            .collect::<Vec<String>>();

        reader.set_threads(self.io_threads)?;

        let mut mod_codes = Vec::new();
        let mut mod_code_idxs = Vec::new();
        let mut mods_per_base = [0u8; 4];
        for (base, codes) in grouped.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            header.push(format!("unmodified_{base}"));
            if self.has_filtering() {
                header.push(format!("fail_unmodified_{base}"));
            }
            mods_per_base[*base as usize] =
                codes.len().try_into().context("too many mods for {base}")?;
            for code in codes {
                mod_codes.push(*code);
                let idx = mod_code_idxs.len();
                mod_code_idxs.push(idx);
                // mod_idxs_per_base[*base as usize].push(idx);
                header.push(format!("modified_{code}"));
                if self.has_filtering() {
                    header.push(format!("fail_modified_{code}"));
                }
            }
            header.push(format!("other_modified_{base}"));
            if self.has_filtering() {
                header.push(format!("fail_other_modified_{base}"));
            }
        }
        header.push("read_length".to_string());

        let (start_idx_per_base, _) = mods_per_base.iter().fold(
            (Vec::new(), 0usize),
            |(mut agg, cum_sum), next| {
                agg.push(cum_sum);
                (agg, cum_sum.saturating_add(*next as usize))
            },
        );
        assert_eq!(start_idx_per_base.len(), 4usize);
        let start_idx_per_base = [
            start_idx_per_base[0],
            start_idx_per_base[1],
            start_idx_per_base[2],
            start_idx_per_base[3],
        ];

        let header = header.join(",");
        let header = format!("{header}\n");

        let (snd_results, rcv_results) = crossbeam::channel::unbounded();
        let (snd_mem, rcv_mem) = crossbeam::channel::unbounded();
        for _ in 0..self.queue_size {
            match snd_mem.send(ReadModsStatsRecord::<16>::empty()) {
                Ok(_) => {}
                Err(_) => {
                    bail!("failed to initiate memory allocation")
                }
            }
        }

        let mut read_processor = ReadModStatsProcessor::<16>::new(
            rcv_mem,
            snd_results,
            reader,
            mods_per_base,
            &mod_codes,
            &mod_code_idxs,
            per_base_thresholds,
            per_mod_thresholds,
            mpb.clone(),
        );

        let mut writer: Box<
            dyn CanWriteReadModStatsRecords<ReadModsStatsRecord<16>>,
        > = {
            match self.out_path.as_str() {
                "-" | "stdout" => Box::new(ReadModStatsWriter::new_stdout(
                    self.has_filtering(),
                    &header,
                    mods_per_base,
                    start_idx_per_base,
                    snd_mem,
                    contigs,
                    mpb.clone(),
                )?),
                _ => {
                    let fp = Path::new(self.out_path.as_str());
                    Box::new(ReadModStatsWriter::new_file(
                        self.has_filtering(),
                        fp,
                        &header,
                        mods_per_base,
                        start_idx_per_base,
                        snd_mem,
                        contigs,
                        mpb.clone(),
                    )?)
                }
            }
        };

        let reader_handle =
            std::thread::spawn(move || read_processor.process_bam());

        let mut errs = FxHashMap::default();
        let done_records = mpb.add(get_ticker());
        done_records.set_message("records processed successfully");

        for result in rcv_results {
            match result {
                Ok(record) => {
                    writer.write(record)?;
                    done_records.inc(1);
                }
                Err(e) => {
                    let c = errs.entry(e.to_string()).or_insert(0usize);
                    *c = c.saturating_add(1);
                }
            }
        }

        match reader_handle.join() {
            Ok(_) => {}
            Err(e) => {
                bail!("reader thread paniced, {e:?}")
            }
        };
        if !errs.is_empty() {
            let table = format_errors_table(&errs);
            mpb.suspend(|| {
                info!("errors:\n{table}");
            })
        }
        Ok(())
    }

    fn parse_thresholds(&self) -> anyhow::Result<[f32; 4]> {
        if self.filter_threshold.is_empty() {
            info!("no thresholds provided, no filtering will be done");
            return Ok([0f32; 4]);
        }
        let (default_threshold, per_base_thresholds) =
            parse_per_base_thresholds(&self.filter_threshold)?;
        let bases_with_thresholds = per_base_thresholds.keys().join(",");
        let bases_without_thresholds =
            [DnaBase::A, DnaBase::C, DnaBase::G, DnaBase::T]
                .into_iter()
                .filter(|x| !per_base_thresholds.contains_key(x))
                .collect::<Vec<_>>();

        match (default_threshold, per_base_thresholds.is_empty()) {
            (Some(x), true) => {
                info!("using threshold value {x} for all primary bases");
                Ok([x; 4])
            }
            (Some(x), false) => {
                if bases_with_thresholds.is_empty() {
                    info!(
                        "using provided thresholds for {bases_with_thresholds}"
                    );
                } else {
                    info!(
                        "using provided thresholds for \
                         {bases_with_thresholds} and {x} for {}",
                        bases_without_thresholds.iter().join(",")
                    );
                }
                let mut thresholds = [x; 4];
                for (base, thresh) in per_base_thresholds {
                    thresholds[base as usize] = thresh;
                }
                Ok(thresholds)
            }
            (None, true) => {
                info!("no thresholds provided, no filtering will be done");
                Ok([0f32; 4])
            }
            (None, false) => {
                if bases_with_thresholds.is_empty() {
                    info!(
                        "using provided thresholds for {bases_with_thresholds}"
                    );
                } else {
                    info!(
                        "using provided thresholds for \
                         {bases_with_thresholds} and 0 for {}",
                        bases_without_thresholds.iter().join(",")
                    );
                }
                let mut thresholds = [0f32; 4];
                for (base, thresh) in per_base_thresholds {
                    thresholds[base as usize] = thresh;
                }
                Ok(thresholds)
            }
        }
    }

    fn parse_per_mod_thresholds(
        &self,
    ) -> anyhow::Result<Vec<(ModCodeRepr, f32)>> {
        if let Some(raw_mod_thresholds) = &self.mod_thresholds.as_ref() {
            let per_mod_thresholds =
                parse_per_mod_thresholds(raw_mod_thresholds)?;
            Ok(per_mod_thresholds.into_iter().collect())
        } else {
            Ok(Vec::new())
        }
    }

    fn has_filtering(&self) -> bool {
        !(self.filter_threshold.is_empty() && self.mod_thresholds.is_none())
    }
}
