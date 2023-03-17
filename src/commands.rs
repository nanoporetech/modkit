use std::collections::{HashMap, HashSet};
use std::io::BufWriter;
use std::num::ParseFloatError;
use std::path::PathBuf;
use std::thread;

use anyhow::{Context, Result as AnyhowResult};
use clap::{Args, Subcommand, ValueEnum};
use crossbeam_channel::bounded;
use indicatif::{
    MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle,
};
use log::{debug, info};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{HeaderView, Read};

use crate::errs::{InputError, RunError};
use crate::interval_chunks::IntervalChunks;
use crate::logging::init_logging;
use crate::mod_bam::{
    collapse_mod_probs, format_mm_ml_tag, CollapseMethod, ModBaseInfo,
    SkipMode, ML_TAGS, MM_TAGS,
};
use crate::mod_base_code::ModCode;
use crate::mod_pileup::{process_region, ModBasePileup, PileupNumericOptions};
use crate::motif_bed::{motif_bed, MotifLocations, RegexMotif};
use crate::summarize::summarize_modbam;
use crate::thresholds::{
    calc_threshold_from_bam, sample_modbase_probs, Percentiles,
};
use crate::util;
use crate::util::{
    add_modkit_pg_records, get_spinner, record_is_secondary, ReferenceRecord,
    Region,
};
use crate::writers::{BedGraphWriter, BedMethylWriter, OutWriter, TsvWriter};

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
    /// silent '.' to '?' or '.'.
    UpdateTags(Update),
    /// Calculate an estimate of the distribution of mod-base prediction
    /// probabilities distribution.
    SampleProbs(SampleModBaseProbs),
    /// Summarize the mod tags present in a BAM and get basic statistics
    Summary(ModSummarize),
    /// Create BED file with all locations of a motif
    MotifBed(MotifBed),
}

impl Commands {
    pub fn run(&self) -> Result<(), String> {
        match self {
            Self::AdjustMods(x) => x.run(),
            Self::Pileup(x) => x.run(),
            Self::SampleProbs(x) => x.run(),
            Self::Summary(x) => x.run(),
            Self::MotifBed(x) => x.run(),
            Self::UpdateTags(x) => x.run(),
        }
    }
}

fn check_collapse_method(raw_method: &str) -> Result<String, String> {
    match raw_method {
        "norm" => Ok(raw_method.to_owned()),
        "dist" => Ok(raw_method.to_owned()),
        _ => Err(format!("unknown method {raw_method}")),
    }
}

#[derive(Args)]
pub struct Adjust {
    /// BAM file to collapse mod call from.
    in_bam: PathBuf,
    /// File path to new BAM file to be created.
    out_bam: PathBuf,
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
    /// Method to use to collapse mod calls, 'norm', 'dist'. A full description
    /// of the methods can be found in collapse.md.
    #[arg(
        long,
        default_value_t = String::from("dist"),
        value_parser = check_collapse_method,
        requires = "ignore",
    )]
    method: String,
    /// Convert one mod-tag to another, summing the probabilities together if
    /// the retained mod tag is already present.
    #[arg(group = "prob_args", long, action = clap::ArgAction::Append, num_args = 2)]
    convert: Option<Vec<char>>,
    /// Output debug logs to file at this path.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
}

type CliResult<T> = Result<T, RunError>;

fn record_is_valid(record: &bam::Record) -> Result<(), RunError> {
    if record_is_secondary(&record) {
        return Err(RunError::new_skipped("not primary"));
    }
    if record.seq_len() == 0 {
        return Err(RunError::new_failed("seq is zero length"));
    }
    Ok(())
}

fn adjust_mod_probs(
    mut record: bam::Record,
    methods: &[CollapseMethod],
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
        for method in methods {
            seq_pos_mod_probs = collapse_mod_probs(seq_pos_mod_probs, method);
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
        .push_aux(mm_style.as_bytes(), mm)
        .expect("failed to add MM tag");
    record
        .push_aux(ml_style.as_bytes(), ml)
        .expect("failed to add ML tag");

    Ok(record)
}

impl Adjust {
    pub fn run(&self) -> Result<(), String> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let fp = &self.in_bam;
        let out_fp = &self.out_bam;
        let mut reader =
            bam::Reader::from_path(fp).map_err(|e| e.to_string())?;
        let threads = self.threads;
        reader.set_threads(threads).map_err(|e| e.to_string())?;
        let mut header = bam::Header::from_template(reader.header());
        add_modkit_pg_records(&mut header);
        let mut out_bam =
            bam::Writer::from_path(out_fp, &header, bam::Format::Bam)
                .map_err(|e| e.to_string())?;

        let fail_fast = self.fail_fast;

        let methods = if let Some(convert) = &self.convert {
            let mut conversions = HashMap::new();
            for chunk in convert.chunks(2) {
                assert_eq!(chunk.len(), 2);
                let from = ModCode::parse_raw_mod_code(chunk[0])?;
                let to = ModCode::parse_raw_mod_code(chunk[1])?;
                let froms = conversions.entry(to).or_insert(HashSet::new());
                froms.insert(from);
            }
            for (to_code, from_codes) in conversions.iter() {
                info!(
                    "Converting {} to {}",
                    from_codes.iter().map(|c| c.char()).collect::<String>(),
                    to_code.char()
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
            let mod_code_to_remove = ModCode::parse_raw_mod_code(self.ignore)?;
            info!(
                "{}",
                format!(
                    "Removing mod base {} from {}, new bam {}",
                    mod_code_to_remove.char(),
                    fp.to_str().unwrap_or("???"),
                    out_fp.to_str().unwrap_or("???")
                )
            );
            let method =
                CollapseMethod::parse_str(&self.method, mod_code_to_remove)?;
            vec![method]
        };

        let spinner = get_spinner();
        spinner.set_message("Adjusting ModBAM");
        let mut total = 0usize;
        let mut total_failed = 0usize;
        let mut total_skipped = 0usize;
        for (i, result) in reader.records().enumerate() {
            if let Ok(record) = result {
                let record_name = util::get_query_name_string(&record)
                    .unwrap_or("???".to_owned());
                match adjust_mod_probs(record, &methods) {
                    Err(RunError::BadInput(InputError(err)))
                    | Err(RunError::Failed(err)) => {
                        if fail_fast {
                            return Err(err.to_string());
                        } else {
                            debug!("read {} failed, {}", record_name, err);
                            total_failed += 1;
                        }
                    }
                    Err(RunError::Skipped(_reason)) => {
                        total_skipped += 1;
                    }
                    Ok(record) => {
                        if let Err(err) = out_bam.write(&record) {
                            if fail_fast {
                                return Err(format!(
                                    "failed to write {}",
                                    err.to_string()
                                ));
                            } else {
                                debug!("failed to write {}", err);
                                total_failed += 1;
                            }
                        } else {
                            spinner.inc(1);
                            total = i;
                        }
                    }
                }
            } else {
                if fail_fast {
                    let err = result.err().unwrap().to_string();
                    return Err(err);
                }
                total_failed += 1;
            }
        }
        spinner.finish_and_clear();

        info!(
            "done, {} records processed, {} failed, {} skipped",
            total + 1,
            total_failed,
            total_skipped
        );
        Ok(())
    }
}

#[derive(Args)]
pub struct ModBamPileup {
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,

    /// Output file to write results into.
    out_bed: PathBuf,

    /// Number of threads to use while processing chunks concurrently.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[arg(long)]
    log_filepath: Option<PathBuf>,

    /// Interval chunk size to process concurrently. Smaller interval chunk
    /// sizes will use less memory but incur more overhead.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,

    /// Sample this fraction of the reads when estimating the
    /// `filter-percentile`. In practice, 50-100 thousand reads is sufficient to
    /// estimate the model output distribution and determine the filtering
    /// threshold.
    #[arg(short = 'f', long, default_value_t = 0.1, hide_short_help = true)]
    sampling_frac: f64,

    /// Set a random seed for deterministic running, the default is non-deterministic.
    #[arg(long, hide_short_help = true)]
    seed: Option<u64>,

    /// Do not perform any filtering, include all mod base calls in output. See
    /// filtering.md for details on filtering.
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,

    /// Filter out mod-calls where the probability of the predicted
    /// variant is below this percentile. For example, 0.1 will filter
    /// out the 10% lowest confidence modification calls.
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,

    /// Filter threshold, drop calls below this probability.
    #[arg(group = "thresholds", long, hide_short_help = true)]
    filter_threshold: Option<f32>,

    /// Collapse _in_situ_ by redistributing base modification probability
    /// equally across other options. For example, if collapsing 'h', with 'm'
    /// and canonical options, half of the probability of 'h' will be added to
    /// both 'm' and 'C'. A full description of the methods can be found in
    /// collapse.md.
    #[arg(long, group = "combine_args", hide = true)]
    collapse: Option<char>,

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

    /// Process only the specified region of the BAM when performing pileup.
    /// Format should be <chrom_name>:<start>-<end>.
    #[arg(long)]
    region: Option<String>,

    /// Only output counts at CpG motifs. Requires a reference sequence to be
    /// provided.
    #[arg(long, requires = "reference_fasta", default_value_t = false)]
    cpg: bool,

    /// Reference sequence in FASTA format. Required for CpG motif filtering.
    #[arg(long = "ref")]
    reference_fasta: Option<PathBuf>,

    /// Prepare data for comparison to whole genome bisulfite sequencing runs.
    /// This setting will ignore modification calls other than 5mC, by removing
    /// them using the redistribute method described in collapse.md.
    #[arg(
        long,
        requires = "cpg",
        conflicts_with = "combine_mods",
        default_value_t = false
    )]
    bisulfite: bool,

    /// Combine mod calls, all counts of modified bases are summed together. See
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
}

impl ModBamPileup {
    fn get_targets(
        &self,
        header: &HeaderView,
        region: Option<Region>,
    ) -> Vec<ReferenceRecord> {
        (0..header.target_count())
            .filter_map(|tid| {
                let chrom_name =
                    String::from_utf8(header.tid2name(tid).to_vec()).unwrap_or("???".to_owned());
                if let Some(region) = &region {
                    if chrom_name == region.name {
                        Some(ReferenceRecord::new(tid, region.start, region.length(), chrom_name))
                    } else {
                        None
                    }
                } else {
                    match header.target_len(tid) {
                        Some(size) => Some(ReferenceRecord::new(tid, 0, size as u32, chrom_name)),
                        None => {
                            debug!("> no size information for {chrom_name} (tid: {tid})");
                            None
                        }
                    }
                }

            })
            .collect::<Vec<ReferenceRecord>>()
    }

    fn run(&self) -> AnyhowResult<(), String> {
        let _handle = init_logging(self.log_filepath.as_ref());
        // do this first so we fail when the file isn't readable
        let header = bam::IndexedReader::from_path(&self.in_bam)
            .map_err(|e| e.to_string())
            .map(|reader| reader.header().to_owned())?;

        let pileup_options = if self.bisulfite {
            if self.combine_mods {
                PileupNumericOptions::Combine
            } else {
                PileupNumericOptions::Collapse(CollapseMethod::ReDistribute(
                    ModCode::h,
                ))
            }
        } else {
            match (self.combine_mods, &self.collapse) {
                (false, None) => PileupNumericOptions::Passthrough,
                (true, _) => PileupNumericOptions::Combine,
                (_, Some(raw_mod_code)) => {
                    let mod_code = ModCode::parse_raw_mod_code(*raw_mod_code)?;
                    let method = CollapseMethod::ReDistribute(mod_code);
                    PileupNumericOptions::Collapse(method)
                }
            }
        };

        let threshold = if self.no_filtering {
            0f32
        } else if let Some(user_threshold) = self.filter_threshold {
            info!(
                "Using user-defined threshold probability: {}",
                user_threshold
            );
            user_threshold
        } else {
            info!(
                "Determining filter threshold probability using sampling \
                frequency {}",
                self.sampling_frac
            );
            // todo need to calc threshold based on collapsed probs
            calc_threshold_from_bam(
                &self.in_bam,
                self.threads,
                self.sampling_frac,
                self.filter_percentile,
                self.seed,
            )?
        };

        info!("Using filter threshold {}", threshold);

        let region = if let Some(raw_region) = &self.region {
            info!("parsing region {raw_region}");
            Some(Region::parse_str(raw_region)?)
        } else {
            None
        };

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .with_context(|| "failed to make threadpool")
            .map_err(|e| e.to_string())?;

        let tids = self.get_targets(&header, region);
        let (motif_locations, tids) = if self.cpg {
            let fasta_fp = self
                .reference_fasta
                .as_ref()
                .ok_or("reference fasta is required for CpG")?;
            let regex_motif = RegexMotif::parse_string("CG", 0).unwrap();
            debug!("filtering output to only CpG motifs");
            if self.combine_strands {
                debug!("combining + and - strand counts");
            }
            let names_to_tid = tids
                .iter()
                .map(|target| (target.name.as_str(), target.tid))
                .collect::<HashMap<&str, u32>>();
            let motif_locations = pool.install(|| {
                MotifLocations::from_fasta(fasta_fp, regex_motif, &names_to_tid)
                    .map_err(|e| e.to_string())
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
        let combine_strands = self.combine_strands;

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
                                        threshold,
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

        let out_fp_str = self.out_bed.clone();
        let mut writer: Box<dyn OutWriter<ModBasePileup>> = if self.bedgraph {
            Box::new(BedGraphWriter::new(out_fp_str)?)
        } else {
            let out_fp = std::fs::File::create(out_fp_str)
                .context("failed to make output file")
                .map_err(|e| e.to_string())?;
            Box::new(BedMethylWriter::new(
                BufWriter::new(out_fp),
                self.only_tabs,
            ))
        };

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
    /// Input BAM with modified base tags.
    in_bam: PathBuf,
    /// Fraction of reads to sample, for example 0.1 will sample
    /// 1/10th of the reads.
    #[arg(short = 'f', long, default_value_t = 0.1)]
    sampling_frac: f64,
    /// Number of threads to use reading BAM.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Random seed for deterministic running, the default is non-deterministic.
    #[arg(short, long)]
    seed: Option<u64>,
    /// Percentiles to calculate, a space separated list of floats.
    #[arg(short, long, default_value_t=String::from("0.1,0.5,0.9"))]
    percentiles: String,
}

impl SampleModBaseProbs {
    fn run(&self) -> AnyhowResult<(), String> {
        let mut bam = bam::Reader::from_path(&self.in_bam).unwrap();
        bam.set_threads(self.threads).unwrap();

        let mut probs =
            sample_modbase_probs(&mut bam, self.seed, self.sampling_frac)
                .map_err(|e| e.to_string())?;
        let desired_percentiles = parse_percentiles(&self.percentiles)
            .with_context(|| {
                format!("failed to parse percentiles: {}", &self.percentiles)
            })
            .map_err(|e| e.to_string())?;
        println!(
            "{}",
            Percentiles::new(&mut probs, &desired_percentiles)?.report()
        );
        Ok(())
    }
}

#[derive(Args)]
pub struct ModSummarize {
    /// Input ModBam file.
    in_bam: PathBuf,
    /// Number of threads to use reading BAM.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
}

impl ModSummarize {
    pub fn run(&self) -> AnyhowResult<(), String> {
        let mod_summary = summarize_modbam(&self.in_bam, self.threads)
            .map_err(|e| e.to_string())?;
        let mut writer = TsvWriter::new();
        writer.write(mod_summary).map_err(|e| e.to_string())?;
        Ok(())
    }
}

#[derive(Args)]
pub struct MotifBed {
    /// Input FASTA file
    fasta: PathBuf,
    /// Motif to search for within FASTA
    motif: String,
    /// Offset within motif.
    offset: usize,
}

impl MotifBed {
    fn run(&self) -> Result<(), String> {
        motif_bed(&self.fasta, &self.motif, self.offset);
        Ok(())
    }
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
