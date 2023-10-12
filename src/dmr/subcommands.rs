use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use bio::io::fasta::Reader as FastaReader;
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use noodles::csi::Index as CsiIndex;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::multi_sample::{
    get_reference_modified_base_positions, n_choose_2, DmrSample,
};
use crate::dmr::pairwise::run_pairwise_dmr;
use crate::dmr::util::{parse_roi_bed, DmrIntervalIter};
use crate::logging::init_logging;
use crate::mod_base_code::DnaBase;
use crate::position_filter::{GenomeLapper, Iv, StrandedPositionFilter};
use crate::util::{
    get_master_progress_bar, get_subroutine_progress_bar, get_ticker,
};

#[derive(Subcommand)]
pub enum BedMethylDmr {
    /// Compare regions in a pair of samples (for example, tumor and normal or
    /// control and experiment). A sample is input as a bgzip pileup bedMethyl
    /// (produced by pileup, for example) that has an associated tabix index.
    /// Output is a BED file with the score column indicating the magnitude of
    /// the difference in methylation between the two samples. See the online
    /// documentation for additional details.
    Pair(PairwiseDmr),
    /// Compare regions between all pairs of samples (for example a trio sample
    /// set or haplotyped trio sample set). As with `pair` all inputs must be bgzip
    /// compressed bedMethyl files with associated tabix indices. Each sample
    /// must be assigned a name. Output is a directory of BED files with the score column
    /// indicating the magnitude of the difference in methylation between the
    /// two samples indicated in the file name. See the online documentation for
    /// additional details.
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
    /// a tabix index with the same name and .tbi next to this file or the --index-a option
    /// must be provided.
    #[arg(short = 'a')]
    control_bed_methyl: PathBuf,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample. There should be
    /// a tabix index with the same name and .tbi next to this file or the --index-b option
    /// must be provided.
    #[arg(short = 'b')]
    exp_bed_methyl: PathBuf,
    /// Path to file to direct output, optional, no argument will direct output to stdout.
    #[arg(short = 'o', long)]
    out_path: Option<String>,
    /// Regions BED file over which to compare methylation levels. Should be tab-separated (spaces
    /// allowed in the "name" column). Requires chrom, chromStart and chromEnd. The Name column is
    /// optional. Strand is currently ignored.
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
    /// Number of threads to use.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// Path to tabix index associated with -a (--control-bed-methyl) bedMethyl file.
    #[arg(long)]
    index_a: Option<PathBuf>,
    /// Path to tabix index associated with -b (--exp-bed-methyl) bedMethyl file.
    #[arg(long)]
    index_b: Option<PathBuf>,
}

impl PairwiseDmr {
    fn get_stranded_position_filter(
        &self,
        name_to_id: Arc<HashMap<String, usize>>,
        multi_pb: &MultiProgress,
        modified_bases: &[DnaBase],
    ) -> anyhow::Result<StrandedPositionFilter> {
        let fasta_reader = FastaReader::from_file(&self.reference_fasta)?;
        let reader_pb = multi_pb.add(get_ticker());
        reader_pb.set_message("sequences read");
        let positions_pb = multi_pb.add(get_ticker());
        positions_pb.set_message("positions found");

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

        let pos_bases = modified_bases
            .iter()
            .map(|b| b.char())
            .collect::<FxHashSet<char>>();
        let neg_bases = modified_bases
            .iter()
            .map(|b| b.complement().char())
            .collect::<FxHashSet<char>>();
        let (pos_positions, neg_positions) = rcv
            .iter()
            .par_bridge()
            .fold(
                || (FxHashMap::default(), FxHashMap::default()),
                |(mut pos_agg, mut neg_agg), (sequence, tid)| {
                    let (pos_strand_intervals, neg_strand_intervals) = sequence.chars().enumerate()
                        .fold((Vec::new(), Vec::new()), |(mut pos_agg, mut neg_agg), (pos, base)| {
                            if pos_bases.contains(&base) {
                                pos_agg.push(Iv {
                                    start: pos as u64,
                                    stop: (pos + 1) as u64,
                                    val: ()
                                })
                            } else if neg_bases.contains(&base) {
                                neg_agg.push(Iv {
                                    start: pos as u64,
                                    stop: (pos + 1) as u64,
                                    val: ()
                                })
                            };
                            (pos_agg, neg_agg)
                        });
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

    fn load_index(
        bedmethyl_path: &PathBuf,
        specified_index: Option<&PathBuf>,
    ) -> anyhow::Result<(CsiIndex, PathBuf)> {
        if let Some(index_fp) = specified_index {
            info!(
                "loading user-specified index at {}",
                index_fp.to_string_lossy()
            );
            noodles::tabix::read(index_fp)
                .with_context(|| {
                    format!("failed to read index at {:?}", index_fp)
                })
                .map(|idx| (idx, index_fp.clone()))
        } else {
            let bedmethyl_fp = bedmethyl_path.to_str().ok_or_else(|| {
                anyhow!("could not format control (a) bedmethyl filepath")
            })?;
            let index_fp = format!("{}.tbi", bedmethyl_fp);
            info!(
                "looking for index associated with {} at {}",
                bedmethyl_fp, &index_fp
            );
            let index_path = Path::new(&index_fp).to_path_buf();

            noodles::tabix::read(&index_fp)
                .context("foo")
                .map(|idx| (idx, index_path))
        }
    }
    fn check_modified_bases(&self) -> anyhow::Result<()> {
        Self::validate_modified_bases(&self.modified_bases)
    }

    fn validate_modified_bases(bases: &[char]) -> anyhow::Result<()> {
        if bases.is_empty() {
            bail!("need to specify at least 1 modified base")
        }
        for b in bases.iter() {
            match *b {
                'A' | 'C' | 'G' | 'T' => {}
                _ => bail!("modified base needs to be A, C, G, or T."),
            }
        }

        Ok(())
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        self.check_modified_bases()?;
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
        let (control_index, _) =
            Self::load_index(&self.control_bed_methyl, self.index_a.as_ref())?;
        let (exp_index, _) =
            Self::load_index(&self.exp_bed_methyl, self.index_b.as_ref())?;

        let writer: Box<dyn Write> = {
            match self.out_path.as_ref() {
                None => Box::new(BufWriter::new(std::io::stdout())),
                Some(fp) => {
                    let p = Path::new(fp);
                    if let Some(parent) = p.parent() {
                        if !parent.exists() {
                            info!(
                                "creating output directory {}",
                                parent.to_str().unwrap_or("failed to parse")
                            );
                            std::fs::create_dir_all(parent)?;
                        }
                    }
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
        info!("loaded {} regions", regions_of_interest.len());

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
            .map(|c| DnaBase::parse(*c))
            .collect::<anyhow::Result<Vec<DnaBase>>>()?;

        let position_filter = self.get_stranded_position_filter(
            control_contig_lookup.clone(),
            &mpb,
            &motifs,
        )?;

        let chunk_size = (self.threads as f32 * 1.5f32).floor() as usize;
        info!("processing {chunk_size} regions concurrently");

        let pb = mpb.add(get_master_progress_bar(regions_of_interest.len()));
        pb.set_message("regions processed");
        let failures = mpb.add(get_ticker());
        failures.set_message("regions failed to process");

        let dmr_interval_iter = DmrIntervalIter::new(
            &self.control_bed_methyl,
            &self.exp_bed_methyl,
            control_contig_lookup.clone(),
            exp_contig_lookup,
            control_index,
            exp_index,
            regions_of_interest.into_iter().collect(),
            chunk_size,
            failures.clone(),
        );

        let success_count = run_pairwise_dmr(
            &self.control_bed_methyl,
            &self.exp_bed_methyl,
            dmr_interval_iter,
            position_filter,
            writer,
            pb,
        )?;

        info!(
            "{} regions processed successfully and {} regions failed",
            success_count,
            failures.position()
        );

        Ok(())
    }
}

#[derive(Args)]
pub struct MultiSampleDmr {
    /// Two or more named samples to compare. Two arguments are required <path> <name>.
    #[arg(short = 's', long = "sample", num_args = 2)]
    samples: Vec<String>,
    /// Optional, paths to tabix indices associated with named samples. Two arguments
    /// are required <path> <name> where <name> corresponds to the name of the sample
    /// given to the -s/--sample argument.
    #[arg(short = 'i', long = "index", num_args = 2)]
    indices: Vec<String>,
    /// Regions BED file over which to compare methylation levels. Should be tab-separated (spaces
    /// allowed in the "name" column). Requires chrom, chromStart and chromEnd. The Name column is
    /// optional. Strand is currently ignored.
    #[arg(long, short = 'r')]
    regions_bed: PathBuf,
    /// Directory to place output DMR results in BED format.
    #[arg(short = 'o', long)]
    out_dir: PathBuf,
    /// Prefix files in directory with this label
    #[arg(short = 'p', long)]
    prefix: Option<String>,
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
    /// Number of threads to use.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
}

impl MultiSampleDmr {
    fn get_writer(
        &self,
        a_name: &str,
        b_name: &str,
    ) -> anyhow::Result<Box<BufWriter<File>>> {
        let fp = if let Some(p) = self.prefix.as_ref() {
            self.out_dir
                .join(format!("{}_{}_{}.bed", p, a_name, b_name))
        } else {
            self.out_dir.join(format!("{}_{}.bed", a_name, b_name))
        };
        if fp.exists() && !self.force {
            bail!(
                "refusing to overwrite {:?}",
                fp.to_str().unwrap_or("failed decode")
            )
        } else {
            let fh = File::create(fp)?;
            Ok(Box::new(BufWriter::new(fh)))
        }
    }

    fn get_stranded_position_filter(
        &self,
        pos_positions: &HashMap<String, Vec<u64>>,
        neg_positions: &HashMap<String, Vec<u64>>,
        name_to_tid: Arc<HashMap<String, usize>>,
    ) -> anyhow::Result<StrandedPositionFilter> {
        let pos_lps = pos_positions
            .iter()
            .filter_map(|(name, positions)| {
                name_to_tid.get(name).map(|tid| {
                    let ivs = positions
                        .iter()
                        .map(|p| Iv {
                            start: *p,
                            stop: *p + 1,
                            val: (),
                        })
                        .collect::<Vec<Iv>>();
                    let lp = GenomeLapper::new(ivs);
                    (*tid as u32, lp)
                })
            })
            .collect::<FxHashMap<u32, GenomeLapper>>();

        let neg_lps = neg_positions
            .iter()
            .filter_map(|(name, positions)| {
                name_to_tid.get(name).map(|&tid| {
                    let ivs = positions
                        .iter()
                        .map(|p| Iv {
                            start: *p,
                            stop: *p + 1,
                            val: (),
                        })
                        .collect::<Vec<Iv>>();
                    let lp = GenomeLapper::new(ivs);
                    (tid as u32, lp)
                })
            })
            .collect::<FxHashMap<u32, GenomeLapper>>();

        Ok(StrandedPositionFilter {
            pos_positions: pos_lps,
            neg_positions: neg_lps,
        })
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let _pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()?;

        PairwiseDmr::validate_modified_bases(&self.modified_bases)?;
        let indices = self.indices
            .chunks(2)
            .filter_map(|raw| {
                if raw.len() != 2 {
                    error!("illegal index pair {:?}, should be length 2 of the form <path> <name>", raw);
                    None
                } else {
                    let fp = Path::new(raw[0].as_str()).to_path_buf();
                    let name = raw[1].to_string();
                    if fp.exists() {
                        Some((name, fp))
                    } else {
                        error!("index for {name} at {} not found", &raw[0]);
                        None
                    }
                }
            })
            .collect::<HashMap<String, PathBuf>>();

        let samples = self.samples
            .chunks(2)
            .filter_map(|raw| {
                if raw.len() != 2 {
                    error!("illegal sample pair {:?}, should be length 2 of the form <path> <name>", raw);
                    None
                } else {
                    let fp = Path::new(raw[0].as_str()).to_path_buf();
                    let name = raw[1].to_string();
                    if fp.exists() {
                        let specified_index = indices.get(&name);
                        if let Ok((_, index_fp)) = PairwiseDmr::load_index(&fp, specified_index) {
                            Some(DmrSample::new(fp, index_fp, name))
                        } else {
                            error!("failed to load tabix index for {name}");
                            None
                        }
                    } else {
                        error!("bedMethyl for {name} at {} not found", &raw[0]);
                        None
                    }
                }
            }).collect::<Vec<DmrSample>>();

        if samples.len() < 2 {
            bail!("failed to collect at least 2 samples");
        }
        let motifs = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<anyhow::Result<Vec<DnaBase>>>()
            .context("failed to parse modified base")?;

        let regions_of_interest = parse_roi_bed(&self.regions_bed)?;
        info!("loaded {} regions", regions_of_interest.len());

        let chunk_size = (self.threads as f32 * 1.5f32).floor() as usize;
        info!("processing {chunk_size} regions concurrently");

        let mpb = MultiProgress::new();
        let sample_pb =
            mpb.add(get_master_progress_bar(n_choose_2(samples.len())?));
        let n_regions = regions_of_interest.len();
        let (positive_positions, negative_positions) =
            get_reference_modified_base_positions(
                &self.reference_fasta,
                &motifs,
                self.mask,
                mpb.clone(),
            )?;

        for pair in samples
            .iter()
            .combinations(2)
            .progress_with(sample_pb.clone())
        {
            let a = pair[0];
            let b = pair[1];
            sample_pb
                .set_message(format!("comparing {} and {}", a.name, b.name));
            let pb = mpb.add(get_subroutine_progress_bar(n_regions));
            pb.set_message("regions processed");
            let failures = mpb.add(get_ticker());
            failures.set_message("regions failed to process");

            let (a_index, _) =
                PairwiseDmr::load_index(&a.bedmethyl_fp, Some(&a.index))?;
            let (b_index, _) =
                PairwiseDmr::load_index(&b.bedmethyl_fp, Some(&b.index))?;
            let control_contig_lookup = a_index
                .header()
                .ok_or_else(|| {
                    anyhow!("failed to get {} tabix header", &a.name)
                })?
                .reference_sequence_names()
                .iter()
                .enumerate()
                .map(|(idx, r)| (r.to_owned(), idx))
                .collect::<HashMap<String, usize>>();
            let control_contig_lookup = Arc::new(control_contig_lookup);

            let exp_contig_lookup = b_index
                .header()
                .ok_or_else(|| {
                    anyhow!("failed to get {} tabix header", b.name)
                })?
                .reference_sequence_names()
                .iter()
                .enumerate()
                .map(|(idx, r)| (r.to_owned(), idx))
                .collect::<HashMap<String, usize>>();

            let position_filter = self.get_stranded_position_filter(
                &positive_positions,
                &negative_positions,
                control_contig_lookup.clone(),
            )?;

            let dmr_interval_iter = DmrIntervalIter::new(
                &a.bedmethyl_fp,
                &b.bedmethyl_fp,
                control_contig_lookup.clone(),
                exp_contig_lookup,
                a_index,
                b_index,
                regions_of_interest.clone().into_iter().collect(),
                chunk_size,
                failures.clone(),
            );

            let writer = self.get_writer(&a.name, &b.name)?;
            let success_count = run_pairwise_dmr(
                &a.bedmethyl_fp,
                &b.bedmethyl_fp,
                dmr_interval_iter,
                position_filter,
                writer,
                pb,
            )?;
            debug!(
                "{} regions processed successfully and {} regions failed for pair {} {}",
                success_count,
                failures.position(),
                &a.name, &b.name,
            );
        }

        Ok(())
    }
}
