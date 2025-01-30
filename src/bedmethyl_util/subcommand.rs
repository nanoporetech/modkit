use std::cmp::Ordering;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use anyhow::bail;
use clap::{Args, Subcommand};
use indicatif::ProgressIterator;
use itertools::Itertools;
use log::{debug, error, info};
use log_once::debug_once;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::command_utils::calculate_chunk_size;
use crate::dmr::bedmethyl::BedMethylLine;
use crate::interval_chunks::{ChromCoordinates, ReferenceIntervalsFeeder};
use crate::logging::init_logging;
use crate::mod_base_code::ModCodeRepr;
use crate::tabix::HtsTabixHandler;
use crate::util::{
    create_out_directory, get_subroutine_progress_bar, get_ticker,
    load_sequence_lengths_file, ReferenceRecord, StrandRule,
};
use crate::writers::BedMethylWriter;
use crate::writers::PileupWriter;

#[derive(Subcommand)]
pub enum EntryBedMethyl {
    /// Perform an outer join on two or more bedMethyl files, summing their
    /// counts for records that overlap
    #[command(name = "merge")]
    MergeBedMethyl(EntryMergeBedMethyl),
}

impl EntryBedMethyl {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryBedMethyl::MergeBedMethyl(x) => x.run(),
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
        let mut parse_fails = 0u64;
        match index.read_bedmethyl(&contig, &range, io_threads) {
            Ok(lines) => {
                let lines = lines.into_iter().filter_map(|r| match r {
                    Ok(record) => Some(record),
                    Err(_) => {
                        parse_fails += 1;
                        None
                    }
                });
                for line in lines {
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
            Err(e) => {
                debug!(
                    "error fetching {}:{}-{} from {:?}, {e}",
                    contig, range.start, range.end, index.indexed_fp
                );
            }
        }
        if parse_fails > 0 {
            debug_once!(
                "failed to parse {parse_fails} bedmethyl records from {:?}",
                index.indexed_fp
            );
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
                            info!(
                                "failed to initialize reader for {bedmethyl:?}"
                            );
                            return None;
                        }
                    };

                Some(index)
            })
            .collect::<Vec<HtsTabixHandler<BedMethylLine>>>();

        if readers.is_empty() {
            bail!(
                "failed to get any bedMethyl readers, are the files \
                 tabix-indexed and bgzip-compressed?"
            );
        }

        // get set of contigs from all files
        // done this way in case one file has a set of contigs that the other
        // bedmethyl files do not have
        let tabix_contigs: HashSet<String> =
            readers.iter().flat_map(|handler| handler.get_contigs()).collect();
        info!(
            "Collected {} contigs from {} readers",
            tabix_contigs.len(),
            readers.len()
        );

        let reference_records = load_sequence_lengths_file(&self.genome_sizes)
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
