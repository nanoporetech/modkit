use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use bio::io::fasta::Reader as FastaReader;
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressIterator};
use log::{debug, error, info};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::dmr::model::ModificationCounts;
use crate::dmr::pairwise::get_modification_counts;
use crate::dmr::util::{parse_roi_bed, DmrIntervalIter};
use crate::logging::init_logging;
use crate::motif_bed::{find_motif_hits, RegexMotif};
use crate::position_filter::{GenomeLapper, Iv, StrandedPositionFilter};
use crate::util::{get_master_progress_bar, get_ticker, Strand};

#[derive(Subcommand)]
pub enum BedMethylDmr {
    Pair(PairwiseDmr),
    Multi(MultiSampleDmr),
}

impl BedMethylDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Self::Pair(x) => x.run(),
            Self::Multi(x) => x.run(),
        }
    }
}

#[derive(Args)]
pub struct PairwiseDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    #[arg(short = 'a')]
    control_bed_methyl: PathBuf,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample. There should be
    /// a tabix index with the same name and .tbi next to this file.
    #[arg(short = 'b')]
    exp_bed_methyl: PathBuf,
    /// Path to file to direct output, optional, no argument will direct output to stdout.
    #[arg(short = 'o', long)]
    out_path: Option<String>,
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
    #[arg(short, alias = "base")]
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

impl PairwiseDmr {
    fn get_stranded_position_filter(
        &self,
        name_to_id: Arc<HashMap<String, usize>>,
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
        for fp in [&self.control_bed_methyl, &self.exp_bed_methyl] {
            if !fp.exists() {
                bail!(
                    "input file {} not found",
                    &self
                        .control_bed_methyl
                        .to_str()
                        .unwrap_or("UTF-8-decode failure")
                )
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
        .context("failed to load 'a' index")?;
        let exp_index = noodles::tabix::read(format!(
            "{}.tbi",
            self.exp_bed_methyl.to_str().unwrap()
        ))
        .context("failed to load 'b' index")?;

        let mut writer: Box<dyn Write> = {
            match self.out_path.as_ref() {
                None => Box::new(BufWriter::new(std::io::stdout())),
                Some(fp) => {
                    let p = Path::new(fp);
                    if p.exists() && !self.force {
                        bail!("refusing to overwrite existing file {}", fp)
                    } else {
                        let fh = File::create(p)?;
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
        let control_contig_lookup = Arc::new(control_contig_lookup);

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

        let chunk_size = (self.threads as f32 * 1.5f32).floor() as usize;
        info!("processing {chunk_size} DMR intervals concurrently");
        let control_fn = self
            .control_bed_methyl
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));
        let exp_fn = self
            .exp_bed_methyl
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));

        let pb = mpb.add(get_master_progress_bar(regions_of_interest.len()));
        pb.set_message("regions processed");
        let failures = mpb.add(get_ticker());
        failures.set_message("regions failed to process");
        let successes = mpb.add(get_ticker());
        successes.set_message("regions processed successfully");

        let dmr_interval_iter = DmrIntervalIter::new(
            control_fn,
            exp_fn,
            control_contig_lookup.clone(),
            exp_contig_lookup,
            control_index,
            exp_index,
            regions_of_interest.into_iter().collect(),
            chunk_size,
            failures.clone(),
        );

        let (snd, rcv) = crossbeam_channel::bounded(1000);

        let control_bedmethyl_fp = self.control_bed_methyl.clone();
        let exp_bedmethyl_fp = self.exp_bed_methyl.clone();
        std::thread::spawn(move || {
            for chunks in dmr_interval_iter {
                let mut result: Vec<anyhow::Result<ModificationCounts>> =
                    vec![];
                let (res, _) = rayon::join(
                    || {
                        chunks
                            .into_par_iter()
                            .map(|dmr_chunk| {
                                get_modification_counts(
                                    &control_bedmethyl_fp,
                                    &exp_bedmethyl_fp,
                                    &dmr_chunk.control_chunks,
                                    &dmr_chunk.exp_chunks,
                                    dmr_chunk.dmr_interval,
                                    &position_filter,
                                    dmr_chunk.chrom_id,
                                )
                            })
                            .collect::<Vec<_>>()
                    },
                    || {
                        result.into_iter().for_each(|counts| {
                            pb.inc(1);
                            match snd.send(counts) {
                                Ok(_) => {}
                                Err(e) => error!(
                                    "failed to send results, {}",
                                    e.to_string()
                                ),
                            }
                        })
                    },
                );
                result = res;
                result.into_iter().for_each(|counts| {
                    pb.inc(1);
                    match snd.send(counts) {
                        Ok(_) => {}
                        Err(e) => {
                            error!("failed to send results, {}", e.to_string())
                        }
                    }
                });
            }
            pb.finish_and_clear();
        });

        for result in rcv {
            match result {
                Ok(counts) => {
                    writer.write(counts.to_row()?.as_bytes())?;
                    successes.inc(1);
                }
                Err(e) => {
                    debug!("unexpected error, {}", e.to_string());
                }
            }
        }

        info!(
            "{} regions processed successfully and {} regions failed",
            successes.position(),
            failures.position()
        );

        Ok(())
    }
}

#[derive(Args)]
pub struct MultiSampleDmr {}

impl MultiSampleDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        unimplemented!()
    }
}
