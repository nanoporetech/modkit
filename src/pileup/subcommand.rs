use crate::command_utils::{
    get_threshold_from_options, parse_per_mod_thresholds, parse_thresholds,
};
use crate::interval_chunks::IntervalChunks;
use crate::logging::init_logging;
use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::mod_base_code::ModCode;
use crate::motif_bed::{MotifLocations, RegexMotif};
use crate::pileup::{process_region, ModBasePileup, PileupNumericOptions};
use crate::position_filter::StrandedPositionFilter;
use crate::util::{
    get_master_progress_bar, get_subroutine_progress_bar, get_targets,
    get_ticker, parse_partition_tags, Region,
};
use crate::writers::{
    BedGraphWriter, BedMethylWriter, OutWriter, PartitioningBedMethylWriter,
};
use anyhow::{anyhow, bail, Context};

use crate::reads_sampler::sampling_schedule::IdxStats;
use clap::{Args, ValueEnum};
use crossbeam_channel::bounded;
use indicatif::{MultiProgress, ParallelProgressIterator};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::io::BufWriter;
use std::path::PathBuf;

#[derive(Args)]
pub struct ModBamPileup {
    // running args
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,
    /// Output file (or directory with --bedgraph option) to write results into.
    /// Specify "-" or "stdout" to direct output to stdout.
    out_bed: String,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Process only the specified region of the BAM when performing pileup.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// Maximum number of records to use when calculating puleup. This argument is
    /// passed to the pileup engine. If you have high depth data, consider
    /// increasing this value substantially. Must be less than 2147483647 or
    /// an error will be raised.
    #[arg(long, default_value_t = 8000, hide_short_help = true)]
    max_depth: u32,

    // processing args
    /// Number of threads to use while processing chunks concurrently.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

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
    /// Interval chunk size in base pairs to process concurrently when estimating the threshold
    /// probability, can be larger than the pileup processing interval.
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,
    /// BED file that will restrict threshold estimation and pileup results to
    /// positions overlapping intervals in the file. (alias: include-positions)
    #[arg(long, hide_short_help = true, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Include unmapped base modifications when estimating the pass threshold
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = false,
        conflicts_with = "include_bed"
    )]
    include_unmapped: bool,

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
    #[arg(long = "ref", alias = "reference", short = 'r')]
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
    /// Discard base modification calls that are this many bases from the start or the end
    /// of the read. For example, a value of 10 will require that the base modification is
    /// at least the 11th base or 11 bases from the end.
    #[arg(long, hide_short_help = true)]
    edge_filter: Option<usize>,

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
    #[arg(long)]
    prefix: Option<String>,
    /// Partition output into multiple bedMethyl files based on tag-value pairs. The output
    /// will be multiple bedMethyl files with the format
    /// `<prefix>_<tag_value_1>_<tag_value_2>_<tag_value_n>.bed` prefix is optional and set
    /// with the `--prefix` flag.
    #[arg(long)]
    partition_tag: Option<Vec<String>>,
}

impl ModBamPileup {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        // do this first so we fail when the file isn't readable
        let header = bam::IndexedReader::from_path(&self.in_bam)
            .map(|reader| reader.header().to_owned())?;

        // options parsing below
        let region = self
            .region
            .as_ref()
            .map(|raw_region| {
                info!("parsing region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let sampling_region = self
            .sample_region
            .as_ref()
            .map(|raw_region| {
                info!("parsing sample region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|trim_num| EdgeFilter::new(*trim_num, *trim_num));
        let per_mod_thresholds = self
            .mod_thresholds
            .as_ref()
            .map(|raw_per_mod_thresholds| {
                parse_per_mod_thresholds(raw_per_mod_thresholds)
            })
            .transpose()?;
        let partition_tags = self
            .partition_tag
            .as_ref()
            .map(|raw_tags| parse_partition_tags(raw_tags))
            .transpose()?;
        let tids = get_targets(&header, region.as_ref());
        let position_filter = self
            .include_bed
            .as_ref()
            .map(|bed_fp| {
                let chrom_to_tid = tids
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
        // use the path here instead of passing the reader directly to avoid potentially
        // changing mutable internal state of the reader.
        IdxStats::new_from_path(
            &self.in_bam,
            region.as_ref(),
            position_filter.as_ref(),
        )
        .and_then(|index_stats| {
            if index_stats.mapped_read_count > 0 {
                Ok(())
            } else {
                Err(anyhow!("did not find any mapped reads, perform alignment first or use \
                modkit extract and/or modkit summary to inspect unaligned modBAMs"))
            }
        })?;

        if self.filter_percentile > 1.0 {
            bail!("filter percentile must be <= 1.0")
        }
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
        let mut writer: Box<dyn OutWriter<ModBasePileup>> =
            match (self.bedgraph, partition_tags.is_some()) {
                (true, _) => Box::new(BedGraphWriter::new(
                    &out_fp_str,
                    self.prefix.as_ref(),
                    partition_tags.is_some(),
                )?),
                (false, true) => Box::new(PartitioningBedMethylWriter::new(
                    &self.out_bed,
                    self.only_tabs,
                    self.prefix.as_ref(),
                )?),
                (false, false) => match out_fp_str.as_str() {
                    "stdout" | "-" => {
                        let writer = BufWriter::new(std::io::stdout());
                        Box::new(BedMethylWriter::new(writer, !self.only_tabs))
                    }
                    _ => {
                        let fh = std::fs::File::create(out_fp_str)
                            .context("failed to make output file")?;
                        let writer = BufWriter::new(fh);
                        Box::new(BedMethylWriter::new(writer, !self.only_tabs))
                    }
                },
            };

        // start the actual work here
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .with_context(|| "failed to make threadpool")?;
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
                        edge_filter.as_ref(),
                        threshold_collapse_method.as_ref(),
                        position_filter.as_ref(),
                        !self.include_unmapped,
                        self.suppress_progress,
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
                    position_filter.as_ref(),
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
        if self.suppress_progress {
            master_progress
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        let tid_progress =
            master_progress.add(get_master_progress_bar(tids.len()));
        tid_progress.set_message("contigs");
        let write_progress = master_progress.add(get_ticker());
        write_progress.set_message("rows written");
        let skipped_reads = master_progress.add(get_ticker());
        skipped_reads.set_message("~records skipped");
        let processed_reads = master_progress.add(get_ticker());
        processed_reads.set_message("~records processed");

        let force_allow = self.force_allow_implicit;
        let max_depth = self.max_depth;

        std::thread::spawn(move || {
            pool.install(|| {
                for target in tids {
                    let intervals = IntervalChunks::new(
                        target.start,
                        target.length,
                        interval_size,
                        target.tid,
                        motif_locations.as_ref(),
                    )
                    .filter(|(start, end)| {
                        position_filter
                            .as_ref()
                            .map(|pf| {
                                pf.overlaps_not_stranded(
                                    target.tid,
                                    *start as u64,
                                    *end as u64,
                                )
                            })
                            .unwrap_or(true)
                    })
                    .collect::<Vec<(u32, u32)>>();
                    let n_intervals = intervals.len();
                    let interval_progress = master_progress
                        .add(get_subroutine_progress_bar(n_intervals));
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
                                        max_depth,
                                        motif_locations.as_ref(),
                                        edge_filter.as_ref(),
                                        partition_tags.as_ref(),
                                        position_filter.as_ref(),
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
                    processed_reads
                        .inc(mod_base_pileup.processed_records as u64);
                    skipped_reads.inc(mod_base_pileup.skipped_records as u64);
                    let rows_written = writer.write(mod_base_pileup)?;
                    write_progress.inc(rows_written);
                }
                Err(message) => {
                    debug!("> unexpected error {message}");
                }
            }
        }
        let rows_processed = write_progress.position();
        let n_skipped_reads = skipped_reads.position();
        let n_skipped_message = if n_skipped_reads == 0 {
            format!("zero reads")
        } else {
            format!("~{n_skipped_reads} reads")
        };
        let n_processed_reads = processed_reads.position();
        write_progress.finish_and_clear();
        processed_reads.finish_and_clear();
        skipped_reads.finish_and_clear();
        info!("Done, processed {rows_processed} rows. Processed ~{n_processed_reads} reads and \
            skipped {n_skipped_message}.");
        Ok(())
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[allow(non_camel_case_types)]
enum Presets {
    traditional,
}
