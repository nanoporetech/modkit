use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};
use clap::Args;
use log::error;
use noodles::bgzf;
use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::model::{llk_ratio, ModificationCounts};
use crate::dmr::{BedMethylLine, DmrInterval};
use crate::logging::init_logging;
use crate::position_filter::Iv;

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

fn get_mod_counts(
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
#[derive(Args)]
pub struct BedMethylDmr {
    control_bed_methyl: PathBuf,
    exp_bed_methyl: PathBuf,
    #[arg(long, short = 'r')]
    regions_bed: PathBuf,
    #[arg(long)]
    log_filepath: Option<PathBuf>,
}

impl BedMethylDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let mut control_reader =
            File::open(&self.control_bed_methyl).map(bgzf::Reader::new)?;
        let mut exp_bed_reader =
            File::open(&self.exp_bed_methyl).map(bgzf::Reader::new)?;
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

        let mut modification_counts_agg = Vec::new();
        for dmr_interval in regions_of_interest {
            let control_chr_id =
                *control_chr_lookup.get(&dmr_interval.chrom).unwrap();
            let exp_chr_id = *exp_chr_lookup.get(&dmr_interval.chrom).unwrap();

            let control_chunks = dmr_interval
                .get_index_chunks(&control_index, control_chr_id)?;
            assert_eq!(control_chunks.len(), 1);
            let (control_methyl_counts, control_total) = get_mod_counts(
                &mut control_reader,
                &control_chunks,
                &dmr_interval.interval,
            )?;

            let exp_chunks =
                dmr_interval.get_index_chunks(&exp_index, exp_chr_id)?;
            let (exp_methyl_counts, exp_total) = get_mod_counts(
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
            modification_counts_agg.push((modification_counts, stat));
        }

        modification_counts_agg
            .sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
        let mut writer = BufWriter::new(std::io::stdout());
        for (counts, stat) in modification_counts_agg {
            writer.write(counts.to_row(stat)?.as_bytes()).unwrap();
        }

        Ok(())
    }
}
