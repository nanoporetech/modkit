use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};
use bio::io::fasta::Reader as FastaReader;
use clap::Args;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use log::{debug, error};
use noodles::bgzf;
use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::model::{AggregatedCounts, ModificationCounts};
use crate::dmr::{BedMethylLine, DmrInterval};
use crate::logging::init_logging;
use crate::motif_bed::{find_motif_hits, RegexMotif};
use crate::position_filter::{GenomeLapper, Iv, StrandedPositionFilter};
use crate::util::{get_master_progress_bar, get_ticker, Strand};

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
    chrom_id: u32,
    position_filter: &StrandedPositionFilter,
) -> anyhow::Result<AggregatedCounts> {
    let grouped_by_position: FxHashMap<u64, Vec<&BedMethylLine>> = bm_lines
        .iter()
        .filter(|bm_line| match bm_line.strand {
            '+' => position_filter.contains(
                chrom_id as i32,
                bm_line.start(),
                Strand::Positive,
            ),
            '-' => position_filter.contains(
                chrom_id as i32,
                bm_line.start(),
                Strand::Negative,
            ),
            '.' => position_filter.overlaps_not_stranded(
                chrom_id,
                bm_line.start(),
                bm_line.stop(),
            ),
            _ => {
                debug!(
                    "encountered illegal strand in bedmethyl {}",
                    bm_line.strand
                );
                false
            }
        })
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
    AggregatedCounts::try_new(counts_per_code, total)
}

fn get_mod_counts_for_condition(
    reader: &mut bgzf::Reader<File>,
    chunks: &[IndexChunk],
    interval: &Iv,
    chrom_id: u32,
    position_filter: &StrandedPositionFilter,
) -> anyhow::Result<AggregatedCounts> {
    let mut bedmethyl_lines = Vec::new();
    // todo could be a fold instead
    for chunk in chunks {
        reader.seek(chunk.start())?;
        let mut lines = Vec::new();
        'readloop: loop {
            let mut buf = String::new();
            let _byts = reader.read_line(&mut buf).unwrap();
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

    aggregate_counts(&bedmethyl_lines, chrom_id, position_filter)
}

fn get_modification_counts(
    control_bedmethyl: &PathBuf,
    exp_bedmethyl: &PathBuf,
    control_chunks: &[IndexChunk],
    exp_chunks: &[IndexChunk],
    dmr_interval: DmrInterval,
    position_filter: &StrandedPositionFilter,
    chrom_id: u32,
) -> anyhow::Result<ModificationCounts> {
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
    let control_counts = get_mod_counts_for_condition(
        &mut control_reader,
        &control_chunks,
        &dmr_interval.interval,
        chrom_id,
        &position_filter,
    )?;
    let experimental_counts = get_mod_counts_for_condition(
        &mut exp_bed_reader,
        &exp_chunks,
        &dmr_interval.interval,
        chrom_id,
        &position_filter,
    )?;

    ModificationCounts::new(
        dmr_interval.start(),
        dmr_interval.stop(),
        control_counts,
        experimental_counts,
        dmr_interval,
    )
}

#[derive(Args)]
pub struct BedMethylDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    control_bed_methyl: PathBuf,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    exp_bed_methyl: PathBuf,
    /// Path to file to direct output, '-' will direct to stdout.
    out_path: String,
    /// Regions over which to compare methylation levels. Should be tab-separated (spaces
    /// allowed in the "name" column). Requires chrom, chromStart, chromEnd, and Name columns.
    /// Strand is currently ignored.
    #[arg(long, short = 'r')]
    regions_bed: PathBuf,
    /// Path to reference fasta for the pileup.
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Bases to use to calculate DMR, may be multiple. For example, to calculate
    /// differentially methylated regions using only cytosine modifications use --base C.
    #[arg(short = 'b', long = "base")]
    modified_bases: Vec<char>,
    /// File to write logs to, it's recommended to use this option.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use, WARNING: currently this will open a file per thread.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
}

impl BedMethylDmr {
    fn get_stranded_position_filter<'a>(
        &self,
        name_to_id: HashMap<String, usize>,
        multi_pb: &MultiProgress,
        modified_bases: &[RegexMotif],
    ) -> anyhow::Result<StrandedPositionFilter> {
        let fasta_reader = FastaReader::from_file(&self.reference_fasta)?;
        let reader_pb = multi_pb.add(get_ticker());
        reader_pb.set_message("sequences read");

        let (snd, rcv) = crossbeam_channel::unbounded();
        let mask = self.mask;

        std::thread::spawn(move || {
            fasta_reader
                .records()
                .progress_with(reader_pb)
                .filter_map(|res| res.ok())
                .filter_map(|record| {
                    name_to_id.get(record.id()).map(|tid| (record, *tid))
                })
                .filter_map(|(record, tid)| {
                    String::from_utf8(record.seq().to_vec())
                        .map(|s| if mask { s } else { s.to_ascii_uppercase() })
                        .ok()
                        .map(|s| (s, tid))
                })
                .for_each(|(sequence, tid)| match snd.send((sequence, tid)) {
                    Ok(_) => {}
                    Err(e) => {
                        error!(
                            "failed to send sequence on channel, {}",
                            e.to_string()
                        );
                    }
                });
        });
        let (pos_positions, neg_positions) = rcv
            .iter()
            .par_bridge()
            .fold(
                || (FxHashMap::default(), FxHashMap::default()),
                |(mut pos_agg, mut neg_agg), (sequence, tid)| {
                    let (pos_strand_intervals, neg_strand_intervals) =
                        modified_bases
                            .iter()
                            .flat_map(|base| {
                                find_motif_hits(sequence.as_str(), base)
                            })
                            .fold(
                                (Vec::<Iv>::new(), Vec::<Iv>::new()),
                                |(mut pos, mut neg), (position, strand)| {
                                    match strand {
                                        Strand::Positive => pos.push(Iv {
                                            start: position as u64,
                                            stop: (position + 1) as u64,
                                            val: (),
                                        }),
                                        Strand::Negative => neg.push(Iv {
                                            start: position as u64,
                                            stop: (position + 1) as u64,
                                            val: (),
                                        }),
                                    }
                                    (pos, neg)
                                },
                            );
                    debug!("found {} positive-strand and {} negative-strand positions on \
                        {tid}", pos_strand_intervals.len(), neg_strand_intervals.len());
                    let pos_lp = GenomeLapper::new(pos_strand_intervals);
                    let neg_lp = GenomeLapper::new(neg_strand_intervals);
                    pos_agg.insert(tid as u32, pos_lp);
                    neg_agg.insert(tid as u32, neg_lp);
                    (pos_agg, neg_agg)
                },
            )
            .reduce(
                || (FxHashMap::default(), FxHashMap::default()),
                |(a_pos, a_neg), (b_pos, b_neg)| {
                    let pos = a_pos
                        .into_iter()
                        .chain(b_pos.into_iter())
                        .collect::<FxHashMap<u32, GenomeLapper>>();
                    let neg = a_neg
                        .into_iter()
                        .chain(b_neg.into_iter())
                        .collect::<FxHashMap<u32, GenomeLapper>>();
                    (pos, neg)
                },
            );

        Ok(StrandedPositionFilter {
            pos_positions,
            neg_positions,
        })
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if self.modified_bases.is_empty() {
            bail!("need to specify at least 1 modified base")
        }
        for b in self.modified_bases.iter() {
            match *b {
                'A' | 'C' | 'G' | 'T' => {}
                _ => bail!("modified base needs to be A, C, G, or T."),
            }
        }
        let _pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()?;
        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        // initial checks
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

        let mut writer: Box<dyn Write> = {
            match self.out_path.as_str() {
                "-" => Box::new(BufWriter::new(std::io::stdout())),
                _ => {
                    let fp = Path::new(&self.out_path);
                    if fp.exists() && !self.force {
                        bail!(
                            "refusing to overwrite existing file {}",
                            &self.out_path
                        )
                    } else {
                        let fh = File::create(fp)?;
                        Box::new(BufWriter::new(fh))
                    }
                }
            }
        };

        let regions_of_interest = parse_roi_bed(&self.regions_bed)?;

        let control_contig_lookup = control_index
            .header()
            .ok_or_else(|| anyhow!("failed to get control tabix header"))?
            .reference_sequence_names()
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        let exp_contig_lookup = exp_index
            .header()
            .ok_or_else(|| anyhow!("failed to get exp tabix header"))?
            .reference_sequence_names()
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        let motifs = self
            .modified_bases
            .iter()
            .map(|c| RegexMotif::parse_string(&format!("{c}"), 0))
            .collect::<anyhow::Result<Vec<RegexMotif>>>()?;

        let position_filter = self.get_stranded_position_filter(
            control_contig_lookup.clone(),
            &mpb,
            &motifs,
        )?;

        let pb = mpb.add(get_master_progress_bar(regions_of_interest.len()));
        let mut modification_counts_agg = regions_of_interest
                .into_par_iter()
                .progress_with(pb)
                .filter_map(|dmr_interval| {
                    match (control_contig_lookup.get(&dmr_interval.chrom), exp_contig_lookup.get(&dmr_interval.chrom)) {
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
                            Some((control_chr_id, control_chunks, exp_chunks, dmr_interval))
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
                .filter_map(|(control_chrom_id, control_chunks,exp_chunks, dmr_interval)| {
                    match get_modification_counts(
                        &self.control_bed_methyl,
                        &self.exp_bed_methyl,
                        &control_chunks,
                        &exp_chunks,
                        dmr_interval.clone(),
                            &position_filter,
                        control_chrom_id as u32,
                    ) {
                        Ok(modification_counts) => Some(modification_counts),
                        Err(e) => {
                            debug!("failed to get modification counts for interval {:?}, {}", dmr_interval, e.to_string());
                            None
                        }
                    }
                }).collect::<Vec<ModificationCounts>>();

        modification_counts_agg
            .sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        // let mut writer = BufWriter::new(std::io::stdout());
        for modification_counts in modification_counts_agg {
            writer
                .write(modification_counts.to_row()?.as_bytes())
                .unwrap();
        }

        Ok(())
    }
}

#[cfg(test)]
mod dmr_unit_tests {
    use crate::dmr::subcommand::parse_roi_bed;
    use crate::dmr::DmrInterval;
    use crate::position_filter::Iv;
    use std::path::Path;

    #[test]
    fn test_roi_parsing() {
        let fp = "tests/resources/sim_cpg_regions.bed";
        let rois = parse_roi_bed(fp).unwrap();
        let expected = [
            DmrInterval {
                interval: Iv {
                    start: 10172120,
                    stop: 10172545,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r1".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10217487,
                    stop: 10218336,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r2".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10034963,
                    stop: 10035266,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r3".to_string(),
            },
        ]
        .to_vec();
        assert_eq!(rois, expected);
    }
}
