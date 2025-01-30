use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};

use anyhow::bail;
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressDrawTarget, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use bigtools::{
    beddata::BedParserStreamingIterator, BigWigWrite, InputSortType,
};

use crate::bedmethyl_util::BedMethylStream;
use crate::command_utils::calculate_chunk_size;
use crate::dmr::bedmethyl::BedMethylLine;
use crate::interval_chunks::{ChromCoordinates, ReferenceIntervalsFeeder};
use crate::logging::init_logging;
use crate::mod_base_code::ModCodeRepr;
use crate::tabix::HtsTabixHandler;
use crate::util::{
    create_out_directory, get_subroutine_progress_bar, get_ticker,
    read_sequence_lengths_file, ReferenceRecord, StrandRule,
};
use crate::writers::BedMethylWriter;
use crate::writers::PileupWriter;

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
}

impl EntryBedMethyl {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryBedMethyl::MergeBedMethyl(x) => x.run(),
            EntryBedMethyl::ToBigWig(x) => x.run(),
        }
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
    #[arg(long, short = 'o', alias = "out")]
    out_bed: String,
    /// TSV of genome sizes, should be <chrom>\t<size_in_bp>
    #[arg(long, short = 'g')]
    genome_sizes: PathBuf,

    /// Force overwrite the output file.
    #[arg(long, default_value_t = false)]
    force: bool,
    /// Output a header with the bedMethyl.
    #[arg(long, alias = "include_header", default_value_t = false)]
    with_header: bool,

    /// Output bedMethyl where the delimiter of columns past column 10 are
    /// space-delimited instead of tab-delimited. This option can be useful
    /// for some browsers and parsers that don't expect the extra columns
    /// of the bedMethyl format.
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
    #[arg(long, hide_short_help = true)]
    chunk_size: Option<usize>,

    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,

    /// Specify a file to write debug logs to.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of tabix/bgzf threads to use.
    #[arg(long, default_value_t = 2)]
    io_threads: usize,
}

type BedMethylChunk = Vec<BedMethylLine>;

fn merge_data(
    readers: &[HtsTabixHandler<BedMethylLine>],
    chrom_coordinates: ChromCoordinates,
    tid_to_name: &FxHashMap<u32, String>,
    io_threads: usize,
) -> anyhow::Result<BedMethylChunk> {
    type Key = (u64, ModCodeRepr, StrandRule);
    // this is safe because of how we constructed this
    let contig = tid_to_name.get(&chrom_coordinates.chrom_tid).unwrap();
    let range = (chrom_coordinates.start_pos as u64)
        ..(chrom_coordinates.end_pos as u64);
    let mut merged_data = FxHashMap::<Key, BedMethylLine>::default();

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

            merged_data
                .entry((line.start(), line.raw_mod_code, line.strand))
                // modify the methyl data if an entry is found
                .and_modify(|methyl| {
                    methyl.count_methylated += line.count_methylated;
                    methyl.valid_coverage += line.valid_coverage;
                    methyl.count_canonical += line.count_canonical;
                    methyl.count_other += line.count_other;
                    methyl.count_delete += line.count_delete;
                    methyl.count_fail += line.count_fail;
                    methyl.count_diff += line.count_diff;
                    methyl.count_nocall += line.count_nocall;
                })
                .or_insert(line);
        }
    }

    // get just the bedmethyllines for writing
    let merged_data = merged_data
        .into_values()
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

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        // setup the writer here so we fail before doing any work (if there are
        // problems).
        let out_fp_str = self.out_bed.clone();
        let mut writer: Box<dyn PileupWriter<Vec<BedMethylLine>>> =
            match out_fp_str.as_str() {
                "stdout" | "-" => {
                    let writer = BufWriter::new(std::io::stdout());
                    Box::new(BedMethylWriter::new(
                        writer,
                        self.mixed_delimiters,
                        self.with_header,
                    )?)
                }
                _ => {
                    create_out_directory(&out_fp_str)?;
                    let fh = if self.force {
                        File::create(out_fp_str)?
                    } else {
                        File::create_new(out_fp_str)?
                    };
                    let writer = BufWriter::new(fh);
                    Box::new(BedMethylWriter::new(
                        writer,
                        self.mixed_delimiters,
                        self.with_header,
                    )?)
                }
            };

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

        let n_contigs = reference_records.len();
        let feeder = ReferenceIntervalsFeeder::new(
            reference_records,
            chunk_size,
            self.interval_size,
            false,
            None,
            None,
        )?;

        let mpb = indicatif::MultiProgress::new();
        let contig_progress = mpb.add(get_subroutine_progress_bar(n_contigs));
        contig_progress.set_message("contigs processed");
        let rows_written = mpb.add(get_ticker());
        rows_written.set_message("merging contigs");
        let errored_batches = mpb.add(get_ticker());
        errored_batches.set_message("batch errors");

        let (snd, rcv) = crossbeam::channel::bounded(1000);
        let io_threads = self.io_threads;
        pool.install(move || {
            feeder
                .into_iter()
                .progress_with(contig_progress)
                .filter_map(|r| match r {
                    Ok(mcc) => Some(mcc),
                    Err(e) => {
                        debug!("failed to retrieve batch, {e}");
                        None
                    }
                })
                .map(|m_cc| {
                    m_cc.into_par_iter()
                        .map(|batch| {
                            batch.0.into_par_iter().map(|chrom_coordinates| {
                                merge_data(
                                    &readers,
                                    chrom_coordinates,
                                    &tid_to_name,
                                    io_threads,
                                )
                            })
                        })
                        .flatten()
                        .collect::<Vec<anyhow::Result<BedMethylChunk>>>()
                })
                .for_each(|batch| match snd.send(batch) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("failed to send on channel, {e}");
                    }
                })
        });

        for batch_result in rcv {
            for result in batch_result {
                match result {
                    Ok(lines) => {
                        rows_written.inc(lines.len() as u64);
                        writer.write(lines, &[])?;
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

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryToBigWig {
    /// Input bedmethyl, uncompressed, "-" or "stdin" indicates an input
    /// stream.
    in_bedmethyl: String,
    /// Output bigWig to make.
    out_fp: PathBuf,
    /// A chromosome sizes file. Each line should be have a chromosome and its
    /// size in bases, separated by whitespace. A fasta index (.fai) works as
    /// well.
    #[arg(long = "sizes", short = 'g')]
    chromsizes: PathBuf,

    /// Make a bigWig track where the values are the percent of bases with this
    /// modification, use multiple comma-separated codes to combine counts. For
    /// example --mod-code m makes a track of the 5mC percentages and
    /// --mod-codes h,m will make a track of the combined counts from 5hmC
    /// and 5mC. Combining counts for different primary bases will cause an
    /// error (e.g. --mod-codes a,h).
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
    #[arg(long, default_value_t = false)]
    negative_strand_values: bool,

    /// Set the number of threads to use. This tool will typically use ~225%
    /// CPU on a HDD. SDDs may be higher. (IO bound)
    #[arg(short = 't', long, default_value_t = 6)]
    pub nthreads: usize,

    /// Set the maximum of zooms to create.
    #[arg(short = 'z', long, default_value_t = 10)]
    pub nzooms: u32,

    /// Set the zoom resolutions to use (overrides the --nzooms argument).
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub zooms: Option<Vec<u32>>,

    /// Don't use compression.
    #[arg(short = 'u', long, default_value_t = false)]
    pub uncompressed: bool,

    /// Number of items to bundle in r-tree.
    #[arg(long, default_value_t = 256)]
    #[arg()]
    pub block_size: u32,

    /// Number of data points bundled at lowest level.
    #[arg(long)]
    #[arg(default_value_t = 1024)]
    pub items_per_slot: u32,

    /// Do not create temporary files for intermediate data.
    #[arg(long, default_value_t = false)]
    pub inmemory: bool,

    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,

    /// Hide the progress bar
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

        let chrom_sizes =
            read_sequence_lengths_file(&self.chromsizes).map(|sizes| {
                sizes
                    .into_iter()
                    .map(|(ch, sz)| (ch, sz as u32))
                    .collect::<HashMap<String, u32>>()
            })?;

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

        let in_stream = BedMethylStream::new(
            in_stream,
            include_codes,
            self.negative_strand_values,
            counter.clone(),
        )?;
        let vals = BedParserStreamingIterator::new(in_stream, false);

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
