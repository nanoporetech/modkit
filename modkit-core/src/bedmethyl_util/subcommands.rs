use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};
use bigtools::{
    beddata::BedParserStreamingIterator, BigWigWrite, InputSortType,
};
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressDrawTarget};
use itertools::Itertools;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use rust_htslib::tbx::{Read as TbxRead, Reader as TbxReader};
use rustc_hash::{FxHashMap, FxHashSet};

use modkit_logging::init_logging;

use crate::bedmethyl_util::BedMethylStream;
use crate::command_utils::calculate_chunk_size;
use crate::dmr::bedmethyl::BedMethylLine;
use crate::dmr::isoform::{
    parse_gtf, transcript_pos0_to_genomic0, GtfId, GtfTranscript,
};
use crate::errs::MkError;
use crate::interval_chunks::{
    ChromCoordinates, ReferenceIntervalBatchesFeeder, TotalLength,
};
use crate::mod_base_code::ModCodeRepr;
use crate::position_filter::Iv;
use crate::tabix::{HtsTabixHandler, ParseBedLine};
use crate::util::{
    create_out_directory, get_guage, get_subroutine_progress_bar, get_ticker,
    read_sequence_lengths_file, ReferenceRecord, StrandRule,
};
use crate::writers::bedmethyl_header;

#[derive(Subcommand)]
pub enum EntryBedMethyl {
    /// Perform an outer join on two or more bedMethyl files, summing their
    /// counts for records that overlap
    #[command(name = "merge")]
    MergeBedMethyl(EntryMergeBedMethyl),
    /// Make a BigWig track from a bedMethyl file or stream.
    /// For details on the BigWig format see https://doi.org/10.1093/bioinformatics/btq351.
    #[command(name = "tobigwig")]
    ToBigWig(EntryToBigWig),
    /// Map a transcriptome-aligned bedMethyl file to genome-coordinates based
    /// on a GTF
    #[command(name = "map-to-genome")]
    MapToGenome(EntryMapToGenome),
}

impl EntryBedMethyl {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryBedMethyl::MergeBedMethyl(x) => x.run(),
            EntryBedMethyl::ToBigWig(x) => x.run(),
            EntryBedMethyl::MapToGenome(x) => x.run(),
        }
    }
}

/// Value of `--min-samples`: either an explicit count or every input ("all").
#[derive(Clone, Debug)]
enum MinSamples {
    All,
    AtLeast(usize),
}

fn parse_min_samples(s: &str) -> Result<MinSamples, String> {
    if s.eq_ignore_ascii_case("all") {
        return Ok(MinSamples::All);
    }
    // A value of 1 is rejected: it would keep every position (an outer join),
    // which is the behaviour when the option is omitted, so it does nothing.
    match s.parse::<usize>() {
        Ok(n) if n > 1 => Ok(MinSamples::AtLeast(n)),
        Ok(_) => Err("must be an integer greater than 1, or \"all\"".to_string()),
        Err(_) => Err(format!(
            "invalid value '{s}': expected an integer greater than 1 or \"all\""
        )),
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryMergeBedMethyl {
    /// Input bedMethyl table(s). Should be bgzip-compressed and have an
    /// associated Tabix index. The tabix index will be assumed to be
    /// $this_file.tbi.
    #[arg(num_args(2..))]
    in_bedmethyl: Vec<PathBuf>,
    /// Specify the output file to write the results table.
    #[clap(help_heading = "Output Options")]
    #[arg(long, short = 'o', alias = "out")]
    out_bed: String,
    /// TSV of genome sizes, should be <chrom>\t<size_in_bp>
    #[arg(long, short = 'g')]
    genome_sizes: PathBuf,

    /// Force overwrite the output file.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    force: bool,
    /// Output a header with the bedMethyl.
    #[clap(help_heading = "Output Options")]
    #[arg(
        long = "header",
        alias = "with-header",
        alias = "include_header",
        default_value_t = false
    )]
    with_header: bool,

    /// Output bedMethyl where the delimiter of columns past column 10 are
    /// space-delimited instead of tab-delimited. This option can be useful
    /// for some browsers and parsers that don't expect the extra columns
    /// of the bedMethyl format.
    #[clap(help_heading = "Output Options")]
    #[arg(
        long = "mixed-delim",
        alias = "mixed-delimiters",
        default_value_t = false,
        hide_short_help = true
    )]
    mixed_delimiters: bool,

    /// Chunk size for how many start..end regions for each chromosome to read.
    /// Larger values will lead to faster merging at the expense of memory
    /// usage, while smaller values will be slower with lower memory usage.
    /// This option will only impact large bedmethyl files.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, hide_short_help = true)]
    chunk_size: Option<usize>,

    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,

    /// Specify a file to write debug logs to.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of batches (of size chunk size) allowed to be in a pre-written
    /// state at once. Increasing this number will increase memory usage.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 30)]
    queue_size: usize,
    /// Number of tabix/bgzf threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 2)]
    io_threads: usize,

    /// Only output a position present in at least this many input bedMethyl
    /// files. Accepts an integer greater than 1, or "all" to require the
    /// position in every input (an inner join across replicates).
    #[clap(help_heading = "Filtering Options")]
    #[arg(long, value_parser = parse_min_samples)]
    min_samples: Option<MinSamples>,
    /// Minimum valid coverage for an input's record to count towards a position,
    /// for both the --min-samples tally and the summed counts.
    #[clap(help_heading = "Filtering Options")]
    #[arg(long)]
    min_sample_coverage: Option<u64>,
}

type BedMethylChunk = Vec<BedMethylLine>;

fn merge_data(
    readers: &[HtsTabixHandler<BedMethylLine>],
    chrom_coordinates: ChromCoordinates,
    tid_to_name: &FxHashMap<u32, String>,
    io_threads: usize,
    min_samples: usize,
    min_sample_coverage: u64,
) -> anyhow::Result<BedMethylChunk> {
    type Key = (u64, ModCodeRepr, StrandRule);
    // this is safe because of how we constructed this
    let contig = tid_to_name.get(&chrom_coordinates.chrom_tid).unwrap();
    let range = (chrom_coordinates.start_pos as u64)
        ..(chrom_coordinates.end_pos as u64);
    // value is the merged record plus a tally of how many inputs contributed to
    // it (each input has at most one record per key), used for the --min-samples
    // (inner-join) filter below.
    let mut merged_data = FxHashMap::<Key, (BedMethylLine, usize)>::default();

    // rationale:
    // iterate over every possible contig
    // iterate over the set chunks of start..end regions for that particular
    // contig obtain lines from each bedmethyl for the
    // contig:start..end insert or update (crude join) existing
    // lines in a hashmap write the hashmap to a new bedmethyl
    // recreate hashmap and repeat process for next contig/regions
    for index in readers.iter() {
        let lines = index.read_bedmethyl(&contig, &range, io_threads)?;

        for line in lines {
            let line = line?;

            // an input only contributes to a position when its record has at
            // least the requested valid coverage
            if line.valid_coverage < min_sample_coverage {
                continue;
            }

            merged_data
                .entry((line.start(), line.raw_mod_code, line.strand))
                // modify the methyl data if an entry is found
                .and_modify(|(methyl, n_samples)| {
                    methyl.count_methylated += line.count_methylated;
                    methyl.valid_coverage += line.valid_coverage;
                    methyl.count_canonical += line.count_canonical;
                    methyl.count_other += line.count_other;
                    methyl.count_delete += line.count_delete;
                    methyl.count_fail += line.count_fail;
                    methyl.count_diff += line.count_diff;
                    methyl.count_nocall += line.count_nocall;
                    *n_samples += 1;
                })
                .or_insert((line, 1));
        }
    }

    // get just the bedmethyllines for writing, keeping only positions present in
    // at least --min-samples inputs (default 1 == outer join, unchanged)
    let merged_data = merged_data
        .into_values()
        .filter(|(_, n_samples)| *n_samples >= min_samples)
        .map(|(methyl, _)| methyl)
        .sorted_by(|a, b| {
            debug_assert_eq!(a.chrom, b.chrom);
            match a.start().cmp(&b.start()) {
                Ordering::Equal => match a.strand.cmp(&b.strand) {
                    Ordering::Equal => a.raw_mod_code.cmp(&b.raw_mod_code),
                    o @ _ => o,
                },
                o @ _ => o,
            }
        })
        .collect::<BedMethylChunk>();

    Ok(merged_data)
}

impl EntryMergeBedMethyl {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

        // Resolve --min-samples ("all" -> number of inputs; omitted -> 1 = outer join).
        let n_inputs = self.in_bedmethyl.len();
        let min_samples: usize = match &self.min_samples {
            None => 1,
            Some(MinSamples::All) => n_inputs,
            Some(MinSamples::AtLeast(n)) => *n,
        };
        if min_samples > n_inputs {
            bail!(
                "--min-samples ({}) is greater than the number of input \
                 bedMethyl files ({})",
                min_samples,
                n_inputs
            );
        }
        let min_sample_coverage: u64 = self.min_sample_coverage.unwrap_or(0);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        // set up the writer here so we fail before doing any work (if there are
        // problems).
        let out_fp_str = self.out_bed.clone();
        let mut writer: Box<BufWriter<dyn Write>> = match out_fp_str.as_str() {
            "stdout" | "-" => {
                let writer = BufWriter::new(std::io::stdout());
                Box::new(writer)
            }
            _ => {
                create_out_directory(&out_fp_str)?;
                let fh = if self.force {
                    File::create(out_fp_str)?
                } else {
                    File::create_new(out_fp_str)?
                };
                let writer = BufWriter::new(fh);
                Box::new(writer)
            }
        };
        if self.with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        }

        let readers = self
            .in_bedmethyl
            .iter()
            .filter_map(|bedmethyl| {
                let index: HtsTabixHandler<BedMethylLine> =
                    match HtsTabixHandler::from_path(&bedmethyl) {
                        Ok(reader) => reader,
                        Err(_) => {
                            return None;
                        }
                    };

                Some(index)
            })
            .collect::<Vec<HtsTabixHandler<BedMethylLine>>>();

        // get set of contigs from all files
        // done this way in case one file has a set of contigs that the other
        // bedmethyl files do not have
        let tabix_contigs: HashSet<String> =
            readers.iter().flat_map(|handler| handler.get_contigs()).collect();
        // get the chrom sizes so we know how large to make the chunks
        let reference_records = read_sequence_lengths_file(&self.genome_sizes)
            .map(|sizes| {
                sizes
                    .into_iter()
                    .enumerate()
                    .filter(|(_, (c, _))| {
                        let in_tabix = tabix_contigs.contains(c);
                        if !in_tabix {
                            debug!(
                                "{c} is not present in tabix headers and will \
                                 be skipped"
                            );
                        }
                        in_tabix
                    })
                    .map(|(tid, (name, length))| {
                        ReferenceRecord::new(tid as u32, 0, length as u32, name)
                    })
                    .collect::<Vec<ReferenceRecord>>()
            })
            .with_context(|| {
                format!(
                    "failed to read genome sizes at {:?}",
                    &self.genome_sizes
                )
            })?;
        let tid_to_name = reference_records
            .iter()
            .map(|rr| (rr.tid, rr.name.to_owned()))
            .collect::<FxHashMap<u32, String>>();

        let chunk_size = calculate_chunk_size(
            self.chunk_size,
            self.interval_size,
            self.threads,
        );

        let feeder = ReferenceIntervalBatchesFeeder::new(
            reference_records,
            chunk_size,
            self.interval_size,
            false,
            None,
            None,
        )?;

        let mpb = MultiProgress::new();
        let contig_progress =
            mpb.add(get_subroutine_progress_bar(feeder.total_length()));
        contig_progress.set_message("genome positions");
        let rows_written = mpb.add(get_ticker());
        rows_written.set_message("rows written");
        let errored_batches = mpb.add(get_ticker());
        errored_batches.set_message("batch errors");

        let (snd, rcv) = crossbeam::channel::bounded(self.queue_size);
        let gauge = mpb.add(get_guage(self.queue_size));
        gauge.set_message("enqueued batches");
        gauge.set_position(snd.len() as u64);

        let io_threads = self.io_threads;
        pool.spawn(move || {
            feeder
                .into_iter()
                .filter_map(|r| match r {
                    Ok(mcc) => Some(mcc),
                    Err(e) => {
                        debug!("failed to retrieve batch, {e}");
                        None
                    }
                })
                .map(|m_cc| {
                    let batch_len = m_cc.total_length();
                    let results = m_cc
                        .into_par_iter()
                        .map(|batch| {
                            batch.0.into_par_iter().map(|chrom_coordinates| {
                                merge_data(
                                    &readers,
                                    chrom_coordinates,
                                    &tid_to_name,
                                    io_threads,
                                    min_samples,
                                    min_sample_coverage,
                                )
                            })
                        })
                        .flatten()
                        .collect::<Vec<anyhow::Result<BedMethylChunk>>>();
                    contig_progress.inc(batch_len);
                    results
                })
                .for_each(|batch| match snd.send(batch) {
                    Ok(_) => {
                        gauge.set_position(snd.len() as u64);
                    }
                    Err(e) => {
                        error!("failed to send on channel, {e}");
                    }
                })
        });

        for batch_result in rcv {
            for result in batch_result {
                match result {
                    Ok(lines) => {
                        let rows = pool.install(move || {
                            lines
                                .into_par_iter()
                                .map(|bml| bml.to_line())
                                .collect::<Vec<String>>()
                        });
                        for row in rows {
                            writer.write(row.as_bytes())?;
                            rows_written.inc(1);
                        }
                    }
                    Err(e) => {
                        debug!("{e}");
                        errored_batches.inc(1);
                    }
                }
            }
        }

        Ok(())
    }
}

#[derive(Debug, clap::Args)]
#[group(required = true, multiple = false)]
pub struct ChromSizes {
    /// A chromosome sizes file. Each line should be a chromosome and its
    /// size in bases, separated by whitespace. A fasta index (.fai) works as
    /// well. Use instead of the bam header.
    #[clap(short = 'g', long = "sizes")]
    chromsizes: Option<PathBuf>,
    /// modBAM from which the pileup was generated. Chromosome sizes are
    /// gathered from the header. Use instead of the chromosome sizes file.
    #[clap(short = 'b', long = "header")]
    input_bam: Option<PathBuf>,
}

impl ChromSizes {
    fn get_sequence_lengths(&self) -> anyhow::Result<HashMap<String, u32>> {
        match (self.chromsizes.as_ref(), self.input_bam.as_ref()) {
            (Some(fp), _) => read_sequence_lengths_file(fp).map(|sizes| {
                sizes
                    .into_iter()
                    .map(|(ch, sz)| (ch, sz as u32))
                    .collect::<HashMap<String, u32>>()
            }),
            (_, Some(fp)) => {
                let reader = bam::Reader::from_path(fp)?;
                let header = reader.header();
                (0..header.target_count())
                    .map(|tid| {
                        header
                            .target_len(tid)
                            .ok_or(anyhow!(
                                "header missing length for tid: {tid}"
                            ))
                            .map(|l| {
                                let contig_name = header
                                    .tid2name(tid)
                                    .iter()
                                    .map(|&b| b as char)
                                    .collect::<String>();
                                (contig_name, l as u32)
                            })
                    })
                    .collect()
            }
            // should be disallowed by Clap
            _ => bail!("chromsizes or header must be provided"),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryToBigWig {
    /// Input bedmethyl, uncompressed, "-" or "stdin" indicates an input
    /// stream.
    in_bedmethyl: String,
    /// Output bigWig filename.
    out_fp: PathBuf,
    #[clap(flatten)]
    chrom_sizes: ChromSizes,
    /// Make a bigWig track where the values are the percent of bases with this
    /// modification, use multiple comma-separated codes to combine counts. For
    /// example --mod-code m makes a track of the 5mC percentages and
    /// --mod-codes h,m will make a track of the combined counts from 5hmC
    /// and 5mC. Combining counts for different primary bases will cause an
    /// error (e.g. --mod-codes a,h will error).
    #[arg(
        short = 'm',
        long,
        value_delimiter = ',',
        required = true,
        alias = "mod-code"
    )]
    mod_codes: Vec<String>,

    /// Report the percentages on the negative strand as negative values. The
    /// data range will be [-100, 100].
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    negative_strand_values: bool,

    /// Set the number of threads to use. This tool will typically use ~225%
    /// CPU on a HDD. SDDs may be higher. (IO bound)
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 6)]
    pub nthreads: usize,

    /// Set the maximum of zooms to create.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'z', long, default_value_t = 10, hide_short_help = true)]
    pub nzooms: u32,

    /// Set the zoom resolutions to use (overrides the --nzooms argument).
    #[clap(help_heading = "Output Options")]
    #[arg(long, value_delimiter = ',', num_args = 1.., hide_short_help = true)]
    pub zooms: Option<Vec<u32>>,

    /// Don't use compression.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'u', long, default_value_t = false, hide_short_help = true)]
    pub uncompressed: bool,

    /// Number of items to bundle in r-tree.
    #[arg(long, default_value_t = 256)]
    #[clap(help_heading = "Output Options", hide_short_help = true)]
    pub block_size: u32,

    /// Number of data points bundled at lowest level.
    #[clap(help_heading = "Output Options", hide_short_help = true)]
    #[arg(long, default_value_t = 1024)]
    pub items_per_slot: u32,

    /// Do not create temporary files for intermediate data.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = false)]
    pub inmemory: bool,

    /// If input bedMethyl has sorting of the same scheme as `sort`, this
    /// option may speed up conversion.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    pub force_chromosome_ordering: bool,

    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,

    /// Hide the progress bar
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
}

impl EntryToBigWig {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());
        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(ProgressDrawTarget::hidden());
        }
        let counter = mpb.add(get_ticker());
        counter.set_message("records processed");

        let include_codes = self
            .mod_codes
            .iter()
            .map(|raw| ModCodeRepr::parse(raw))
            .collect::<anyhow::Result<FxHashSet<ModCodeRepr>>>()?;

        if include_codes.is_empty() {
            bail!("must provide at least one modification code to use")
        }

        let chrom_sizes = self.chrom_sizes.get_sequence_lengths()?;
        mpb.suspend(|| {
            info!("loaded {} chromosomes", chrom_sizes.len());
        });

        let mut outb = BigWigWrite::create_file(&self.out_fp, chrom_sizes)?;
        outb.options.max_zooms = self.nzooms;
        outb.options.manual_zoom_sizes = self.zooms.clone();
        outb.options.compress = !self.uncompressed;
        outb.options.input_sort_type = InputSortType::ALL;
        outb.options.block_size = self.block_size;
        outb.options.inmemory = self.inmemory;

        let in_stream: Box<dyn BufRead> = match self.in_bedmethyl.as_str() {
            "-" | "stdin" => Box::new(BufReader::new(std::io::stdin().lock())),
            p @ _ => {
                let fp = Path::new(p);
                Box::new(BufReader::new(File::open(fp)?))
            }
        };

        let in_stream = match BedMethylStream::new(
            in_stream,
            include_codes,
            self.negative_strand_values,
            counter.clone(),
        ) {
            Ok(x) => x,
            Err(MkError::EmptyInput) => {
                mpb.suspend(|| {
                    warn!("empty input, producing empty bigwig");
                });
                return Ok(());
            }
            Err(e) => bail!("{e}"),
        };
        let vals = BedParserStreamingIterator::new(
            in_stream,
            !self.force_chromosome_ordering,
        );

        let rt = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(self.nthreads)
            .build()?;

        outb.write(vals, rt)?;
        let message = format!("finished, wrote {} records", counter.position());

        if self.suppress_progress {
            debug!("{message}");
        } else {
            info!("{message}");
        }
        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryMapToGenome {
    input_bedmethyl: PathBuf,
    output_bedmethyl: String,
    #[arg(long, short = 't', alias = "tx")]
    transcript_id: String,
    #[arg(long)]
    transcript_version: Option<u32>,
    #[arg(long)]
    gtf: PathBuf,
    #[arg(long, default_value_t = false)]
    header: bool,
}

impl EntryMapToGenome {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(None);

        let multi_progress = MultiProgress::new();
        let tx_models = parse_gtf(&self.gtf, &multi_progress)?
            .into_iter()
            .collect::<FxHashMap<_, _>>();

        let tx_id = if let Some(tx_version) = self.transcript_version {
            GtfTranscript::new(self.transcript_id.clone(), tx_version)
        } else {
            GtfTranscript::from_str(&self.transcript_id)
                .context("failed to parse transcript id")?
        };
        let Some(tm) = tx_models.get(&tx_id) else {
            bail!("didn't find {} in GTF", &self.transcript_id)
        };

        let processed_records = multi_progress.add(get_ticker());
        processed_records.set_message("processed records");
        let errored = multi_progress.add(get_ticker());
        errored.set_message("errored records");

        let mut reader = TbxReader::from_path(&self.input_bedmethyl)
            .context("failed to get tabix reader")?;
        let Some(tid) = reader
            .seqnames()
            .iter()
            .find(|x| x.starts_with(&self.transcript_id))
            .and_then(|contig_name| {
                multi_progress.suspend(|| {
                    info!(
                        "found contig {contig_name} for transcript-id {}",
                        self.transcript_id
                    );
                });
                reader.tid(contig_name.as_str()).ok()
            })
        else {
            bail!("didn't find {} in tabix header", self.transcript_id)
        };

        let mut writer: Box<dyn Write> = match self.output_bedmethyl.as_str() {
            "-" | "stdout" => Box::new(BufWriter::new(stdout().lock())),
            p @ _ => Box::new(BufWriter::new(File::create(p)?)),
        };
        if self.header {
            writer.write(bedmethyl_header().as_bytes())?;
        }

        reader.fetch(tid, 0, tm.transcript_len)?;
        for res in reader.records() {
            let Ok(mut bml) = res
                .map_err(|e| MkError::HtsLibError(e))
                .and_then(|bs| {
                    String::from_utf8(bs)
                        .map_err(|e| MkError::InvalidBedMethyl(e.to_string()))
                })
                .and_then(|raw| BedMethylLine::parse(&raw))
            else {
                errored.inc(1);
                continue;
            };
            let Ok(genome_start) =
                transcript_pos0_to_genomic0(&tm, bml.start())
            else {
                errored.inc(1);
                continue;
            };
            let genome_stop = match bml.strand {
                StrandRule::Positive | StrandRule::Both => {
                    genome_start.saturating_add(1)
                }
                StrandRule::Negative => genome_start.saturating_sub(1),
            };

            bml.chrom = tm.chrom.clone();
            bml.interval =
                Iv { start: genome_start, stop: genome_stop, val: () };
            writer.write(bml.to_line().as_bytes())?;
            processed_records.inc(1);
        }

        multi_progress.suspend(|| {
            info!(
                "done, processed {} records, {} failed",
                processed_records.position(),
                errored.position()
            );
        });

        Ok(())
    }
}
