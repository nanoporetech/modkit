use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};
use clap::Args;
use indicatif::ParallelProgressIterator;
use log::{debug, error};
use noodles::bgzf;
use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::model::{llk_ratio, ModificationCounts};
use crate::dmr::{BedMethylLine, DmrInterval};
use crate::logging::init_logging;
use crate::position_filter::Iv;
use crate::util::get_master_progress_bar;

fn parse_roi_bed<P: AsRef<Path>>(fp: P) -> anyhow::Result<Vec<DmrInterval>> {
    let intervals = BufReader::new(File::open(fp)?)
        .lines()
        .filter_map(|r| match r {
            Ok(l) => Some(l),
            Err(e) => {
                error!(
                    "error fetching line from regions BED, {}",
                    e.to_string()
                );
                None
            }
        })
        .map(|line| DmrInterval::parse_str(&line))
        .collect::<anyhow::Result<Vec<DmrInterval>>>()?;
    if intervals.is_empty() {
        bail!("didn't parse any regions")
    } else {
        Ok(intervals)
    }
}

fn aggregate_counts(
    bm_lines: &[BedMethylLine],
) -> (HashMap<char, usize>, usize) {
    let grouped_by_position: FxHashMap<u64, Vec<&BedMethylLine>> = bm_lines
        .iter()
        .fold(FxHashMap::default(), |mut acc, bm_line| {
            acc.entry(bm_line.start())
                .or_insert(Vec::new())
                .push(bm_line);
            acc
        });
    let (counts_per_code, total) = grouped_by_position.into_iter().fold(
        (HashMap::new(), 0),
        |(mut acc, mut total_so_far), (_pos, grouped)| {
            let valid_covs = grouped
                .iter()
                .map(|bml| bml.valid_coverage)
                .collect::<FxHashSet<u64>>();
            let chroms = grouped
                .iter()
                .map(|bml| &bml.chrom)
                .collect::<FxHashSet<&String>>();
            let valid_coverage = grouped[0].valid_coverage as usize;
            assert_eq!(valid_covs.len(), 1);
            assert_eq!(
                chroms.len(),
                1,
                "should only get 1 chrom, got {} {:?}, {:?}",
                chroms.len(),
                &chroms,
                &grouped
            );
            for x in grouped {
                *acc.entry(x.raw_mod_code).or_insert(0) +=
                    x.count_methylated as usize;
            }
            total_so_far += valid_coverage;
            (acc, total_so_far)
        },
    );
    (counts_per_code, total)
}

struct AggregatedCounts {
    mod_code_counts: HashMap<char, usize>,
    total: usize,
}

fn get_mod_counts_for_condition(
    reader: &mut bgzf::Reader<File>,
    chunks: &[IndexChunk],
    interval: &Iv,
) -> anyhow::Result<(HashMap<char, usize>, usize)> {
    let mut bedmethyl_lines = Vec::new();
    // todo could be a fold instead
    for chunk in chunks {
        reader.seek(chunk.start())?;
        let mut lines = Vec::new();
        'readloop: loop {
            let mut buf = String::new();
            let _byts = reader.read_line(&mut buf).unwrap();
            assert!(buf.starts_with("chr"), "illegal line? {}", buf);
            lines.push(buf);
            let cur_pos = reader.virtual_position();
            if cur_pos >= chunk.end() {
                break 'readloop;
            }
        }
        bedmethyl_lines.extend(
            lines
                .into_iter()
                .filter_map(|l| BedMethylLine::parse(l.as_str()).ok())
                .filter(|bml| interval.overlap(bml.start(), bml.stop())),
        );
    }

    Ok(aggregate_counts(&bedmethyl_lines))
}

fn get_modification_counts(
    control_bedmethyl: &PathBuf,
    exp_bedmethyl: &PathBuf,
    control_chunks: &[IndexChunk],
    exp_chunks: &[IndexChunk],
    dmr_interval: DmrInterval,
) -> anyhow::Result<(ModificationCounts, f64)> {
    let mut control_reader =
        File::open(control_bedmethyl).map(bgzf::Reader::new)?;
    let mut exp_bed_reader =
        File::open(exp_bedmethyl).map(bgzf::Reader::new)?;
    if control_chunks.len() != 1 {
        debug!("more than 1 control chunk?");
    }
    if exp_chunks.len() != 1 {
        debug!("more than 1 control chunk?");
    }
    let (control_methyl_counts, control_total) = get_mod_counts_for_condition(
        &mut control_reader,
        &control_chunks,
        &dmr_interval.interval,
    )?;
    let (exp_methyl_counts, exp_total) = get_mod_counts_for_condition(
        &mut exp_bed_reader,
        &exp_chunks,
        &dmr_interval.interval,
    )?;

    let modification_counts = ModificationCounts::new(
        dmr_interval.start(),
        dmr_interval.stop(),
        control_methyl_counts,
        control_total,
        exp_methyl_counts,
        exp_total,
        dmr_interval,
    );
    let stat = llk_ratio(&modification_counts)?;
    Ok((modification_counts, stat))
}

#[derive(Args)]
pub struct BedMethylDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    control_bed_methyl: PathBuf,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    exp_bed_methyl: PathBuf,
    /// Regions over which to compare methylation levels. Should be tab-separated (spaces
    /// allowed in the "name" column). Requires chrom, chromStart, chromEnd, and Name columns.
    /// Strand is currently ignored.
    #[arg(long, short = 'r')]
    regions_bed: PathBuf,
    /// File to write logs to, it's recommended to use this option.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use, WARNING: currently this will open a file per thread.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
}

impl BedMethylDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let control_index = noodles::tabix::read(format!(
            "{}.tbi",
            self.control_bed_methyl.to_str().unwrap()
        ))
        .context("failed to load control index")?;
        let exp_index = noodles::tabix::read(format!(
            "{}.tbi",
            self.exp_bed_methyl.to_str().unwrap()
        ))
        .context("failed to load control index")?;
        let regions_of_interest = parse_roi_bed(&self.regions_bed)?;

        let control_chr_lookup = control_index
            .header()
            .ok_or_else(|| anyhow!("failed to get control tabix header"))?
            .reference_sequence_names()
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        let exp_chr_lookup = exp_index
            .header()
            .ok_or_else(|| anyhow!("failed to get exp tabix header"))?
            .reference_sequence_names()
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        let _pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()?;

        let pb = get_master_progress_bar(regions_of_interest.len());
        let mut modification_counts_agg = regions_of_interest
                .into_par_iter()
                .progress_with(pb)
                .filter_map(|dmr_interval| {
                    match (control_chr_lookup.get(&dmr_interval.chrom), exp_chr_lookup.get(&dmr_interval.chrom)) {
                        (Some(control_chr_id), Some(exp_chr_id)) => {
                            Some((*control_chr_id, *exp_chr_id, dmr_interval))
                        },
                        (None, _) => {
                            debug!("didn't find chrom id for {} in control tabix header", &dmr_interval.chrom);
                            None
                        },
                        (_, None) => {
                            debug!("didn't find chrom id for {} in experimental tabix header", &dmr_interval.chrom);
                            None
                        }
                    }
                })
                .filter_map(|(control_chr_id, exp_chr_id, dmr_interval)| {
                    let control_chunks = dmr_interval
                        .get_index_chunks(&control_index, control_chr_id);
                    let exp_chunks =
                        dmr_interval.get_index_chunks(&exp_index, exp_chr_id);
                    match (control_chunks, exp_chunks) {
                        (Ok(control_chunks), Ok(exp_chunks)) => {
                            Some((control_chunks, exp_chunks, dmr_interval))
                        },
                        (Err(e), _) => {
                            debug!("failed to index into control bedMethyl for chrom id {}, {}", control_chr_id, e.to_string());
                            None
                        },
                        (_, Err(e)) => {
                            debug!("failed to index into experiment bedMethyl for chrom id {}, {}",exp_chr_id, e.to_string());
                            None
                        }
                    }
                })
                .filter_map(|(control_chunks,exp_chunks, dmr_interval)| {
                    match get_modification_counts(
                        &self.control_bed_methyl,
                        &self.exp_bed_methyl,
                        &control_chunks,
                        &exp_chunks,
                        dmr_interval.clone()
                    ) {
                        Ok(modification_counts) => Some(modification_counts),
                        Err(e) => {
                            debug!("failed to get modification counts for interval {:?}, {}", dmr_interval, e.to_string());
                            None
                        }
                    }
                }).collect::<Vec<(ModificationCounts, f64)>>();

        modification_counts_agg
            .sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
        let mut writer = BufWriter::new(std::io::stdout());
        for (counts, stat) in modification_counts_agg {
            writer.write(counts.to_row(stat)?.as_bytes()).unwrap();
        }

        Ok(())
    }
}
