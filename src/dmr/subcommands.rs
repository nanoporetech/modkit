use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use rustc_hash::FxHashMap;

use crate::dmr::pairwise::run_pairwise_dmr;
use crate::dmr::single_site::SingleSiteDmrAnalysis;
use crate::dmr::tabix::{IndexHandler, MultiSampleIndex};
use crate::dmr::util::{parse_roi_bed, HandleMissing, RoiIter};
use crate::genome_positions::GenomePositions;
use crate::logging::init_logging;
use crate::mod_base_code::{DnaBase, ModCodeRepr, MOD_CODE_TO_DNA_BASE};
use crate::util::{
    create_out_directory, get_master_progress_bar, get_subroutine_progress_bar,
    get_ticker,
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
    /// set or haplotyped trio sample set). As with `pair` all inputs must be
    /// bgzip compressed bedMethyl files with associated tabix indices.
    /// Each sample must be assigned a name. Output is a directory of BED
    /// files with the score column indicating the magnitude of the
    /// difference in methylation between the two samples indicated in the
    /// file name. See the online documentation for additional details.
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
#[command(arg_required_else_help = true)]
pub struct PairwiseDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There
    /// should be a tabix index with the same name and .tbi next to this
    /// file or the --index-a option must be provided.
    #[arg(short = 'a')]
    control_bed_methyl: Vec<PathBuf>,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample.
    /// There should be a tabix index with the same name and .tbi next to
    /// this file or the --index-b option must be provided.
    #[arg(short = 'b')]
    exp_bed_methyl: Vec<PathBuf>,
    /// Path to file to direct output, optional, no argument will direct output
    /// to stdout.
    #[arg(short = 'o', long)]
    out_path: Option<String>,
    /// Include header in output
    #[arg(long, default_value_t = false)]
    header: bool,
    /// BED file of regions over which to compare methylation levels. Should be
    /// tab-separated (spaces allowed in the "name" column). Requires
    /// chrom, chromStart and chromEnd. The Name column is optional. Strand
    /// is currently ignored. When omitted, methylation levels are compared at
    /// each site.
    #[arg(long, short = 'r', alias = "regions")]
    regions_bed: Option<PathBuf>,
    /// Path to reference fasta for used in the pileup/alignment.
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Run segmentation, output segmented differentially methylated regions to
    /// this file.
    #[arg(long = "segment", conflicts_with = "regions_bed")]
    segmentation_fp: Option<PathBuf>,

    /// Maximum number of base pairs between modified bases for them to be
    /// segmented together.
    #[arg(long, requires = "segmentation_fp", default_value_t = 5000)]
    max_gap_size: u64,
    /// Prior probability of a differentially methylated position
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.1,
        hide_short_help = true
    )]
    dmr_prior: f64,
    /// Maximum probability of continuing a differentially methylated block,
    /// decay will be dynamic based on proximity to the next position.
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.9,
        hide_short_help = true
    )]
    diff_stay: f64,
    /// Significance factor, effective p-value necessary to favor the
    /// "Different" state.
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.01,
        hide_short_help = true
    )]
    significance_factor: f64,
    /// Use logarithmic decay for "Different" stay probability
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = false,
        hide_short_help = true
    )]
    log_transition_decay: bool,
    /// After this many base pairs, the transition probability will become the
    /// prior probability of encountering a differentially modified
    /// position.
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 500,
        hide_short_help = true
    )]
    decay_distance: u32,
    /// Preset HMM segmentation parameters for higher propensity to switch from
    /// "Same" to "Different" state. Results will be shorter segments, but
    /// potentially higher sensitivity.
    #[arg(
        long,
        requires = "segmentation_fp",
        conflicts_with_all=["log_transition_decay", "significance_factor", "diff_stay", "dmr_prior"],
        default_value_t=false
    )]
    fine_grained: bool,
    /// Bases to use to calculate DMR, may be multiple. For example, to
    /// calculate differentially methylated regions using only cytosine
    /// modifications use --base C.
    #[arg(short, long="base", alias = "modified-bases", action=clap::ArgAction::Append)]
    modified_bases: Vec<char>,
    /// Extra assignments of modification codes to their respective primary
    /// bases. In general, modkit dmr will use the SAM specification to
    /// know which modification codes are appropriate to use for a given
    /// primary base. For example "h" is the code for 5hmC, so is appropriate
    /// for cytosine bases, but not adenine bases. However, if your
    /// bedMethyl file contains custom codes or codes that are not part of
    /// the specification, you can specify which primary base they
    /// belong to here with --assign-code x:C meaning associate modification
    /// code "x" with cytosine (C) primary sequence bases. If a code is
    /// encountered that is not part of the specification, the bedMethyl
    /// record will not be used, this will be logged.
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,

    /// File to write logs to, it's recommended to use this option.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Control the  batch size. The batch size is the number of regions to
    /// load at a time. Each region will be processed concurrently. Loading
    /// more regions at a time will decrease IO to load data, but will use
    /// more memory. Default will be 50% more than the number of
    /// threads assigned.
    #[arg(long, alias = "batch")]
    batch_size: Option<usize>,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// Path to tabix index associated with -a (--control-bed-methyl) bedMethyl
    /// file.
    #[arg(long)]
    index_a: Option<Vec<PathBuf>>,
    /// Path to tabix index associated with -b (--exp-bed-methyl) bedMethyl
    /// file.
    #[arg(long)]
    index_b: Option<Vec<PathBuf>>,
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::warn)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
    /// Prior distribution for estimating MAP-based p-value. Should be two
    /// arguments for alpha and beta (e.g. 1.0 1.0). See
    /// `dmr_scoring_details.md` for additional details on how the metric
    /// is calculated.
    #[arg(
        long,
        num_args = 2,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    prior: Option<Vec<f64>>,
    /// Consider only effect sizes greater than this when calculating the
    /// MAP-based p-value.
    #[arg(
        long,
        default_value_t = 0.05,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    delta: f64,
    /// Sample this many reads when estimating the max coverage thresholds.
    #[arg(
        long,
        short='N',
        default_value_t = 10_042,
        conflicts_with_all = ["max_coverages", "regions_bed"],
    )]
    n_sample_records: usize,
    /// Max coverages to enforce when calculating estimated MAP-based p-value.
    #[arg(long, num_args = 2, conflicts_with = "regions_bed")]
    max_coverages: Option<Vec<usize>>,
    /// When using replicates, cap coverage to be equal to the maximum coverage
    /// for a single sample. For example, if there are 3 replicates with
    /// max_coverage of 30, the total coverage would normally be 90. Using
    /// --cap-coverages will down sample the data to 30X.
    #[arg(
        long,
        conflicts_with = "regions_bed",
        default_value_t = false,
        hide_short_help = true
    )]
    cap_coverages: bool,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u64,
}

impl PairwiseDmr {
    fn check_modified_bases(
        &self,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        Self::validate_modified_bases(
            &self.modified_bases,
            self.mod_code_assignments.as_ref(),
        )
    }

    fn is_single_site(&self) -> bool {
        self.regions_bed.is_none()
    }

    fn parse_raw_assignments(
        raw_mod_code_assignments: Option<&Vec<String>>,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        if let Some(raw_assignments) = raw_mod_code_assignments {
            let user_assignments = raw_assignments.iter().try_fold(
                FxHashMap::default(),
                |mut acc, next| {
                    if next.contains(':') {
                        let parts = next.split(':').collect::<Vec<&str>>();
                        if parts.len() != 2 {
                            bail!(
                                "invalid assignment {next}, should be \
                                 <code>:<DNA>, such as m:C"
                            )
                        } else {
                            let dna_base = parts[1]
                                .parse::<char>()
                                .map_err(|e| {
                                    anyhow!(
                                        "invalid DNA base, should be single \
                                         letter, {e}"
                                    )
                                })
                                .and_then(|raw| DnaBase::parse(raw))?;
                            let mod_code = ModCodeRepr::parse(parts[0])?;
                            debug!(
                                "assigning modification code {mod_code:?} to \
                                 {dna_base:?}"
                            );
                            acc.insert(mod_code, dna_base);
                            Ok(acc)
                        }
                    } else {
                        bail!(
                            "invalid assignment {next}, should be \
                             <code>:<DNA>, such as m:C"
                        )
                    }
                },
            )?;
            Ok(MOD_CODE_TO_DNA_BASE
                .clone()
                .into_iter()
                .chain(user_assignments.into_iter())
                .collect())
        } else {
            Ok(MOD_CODE_TO_DNA_BASE.clone())
        }
    }

    fn validate_modified_bases(
        bases: &[char],
        raw_mod_code_assignments: Option<&Vec<String>>,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        if bases.is_empty() {
            bail!("need to specify at least 1 modified base")
        }
        for b in bases.iter() {
            match *b {
                'A' | 'C' | 'G' | 'T' => {
                    debug!("using primary sequence base {b}");
                }
                _ => bail!("modified base needs to be A, C, G, or T."),
            }
        }

        Self::parse_raw_assignments(raw_mod_code_assignments)
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        if self.control_bed_methyl.is_empty() || self.exp_bed_methyl.is_empty()
        {
            bail!("need to provide at least 1 'a' sample and 'b' sample")
        }
        let code_lookup = self.check_modified_bases()?;

        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let modified_bases = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<anyhow::Result<Vec<DnaBase>>>()?;

        if self.regions_bed.is_some()
            & (self.control_bed_methyl.len() > 1
                || self.exp_bed_methyl.len() > 1)
        {
            bail!(
                "in order to perform multiple comparisons over regions use \
                 modkit dmr multi"
            )
        }

        let a_handlers = self
            .control_bed_methyl
            .iter()
            .enumerate()
            .map(|(i, fp)| {
                let index_fp =
                    self.index_a.as_ref().and_then(|indexes| indexes.get(i));
                IndexHandler::new_infer_index_filepath(fp, index_fp)
            })
            .collect::<anyhow::Result<Vec<IndexHandler>>>()?;
        let b_handlers = self
            .exp_bed_methyl
            .iter()
            .enumerate()
            .map(|(i, fp)| {
                let index_fp =
                    self.index_b.as_ref().and_then(|indexes| indexes.get(i));
                IndexHandler::new_infer_index_filepath(fp, index_fp)
            })
            .collect::<anyhow::Result<Vec<IndexHandler>>>()?;
        let handlers = a_handlers
            .into_iter()
            .chain(b_handlers)
            .collect::<Vec<IndexHandler>>();

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
        );

        let writer: Box<dyn Write> = {
            match self.out_path.as_ref() {
                None => Box::new(BufWriter::new(std::io::stdout())),
                Some(fp) => {
                    let p = Path::new(fp);
                    create_out_directory(p)?;
                    if p.exists() && !self.force {
                        bail!("refusing to overwrite existing file {}", fp)
                    } else {
                        let fh = File::create(p)?;
                        Box::new(BufWriter::new(fh))
                    }
                }
            }
        };

        info!("reading reference FASTA at {:?}", self.reference_fasta);
        let genome_positions = GenomePositions::new_from_sequences(
            &modified_bases,
            &self.reference_fasta,
            self.mask,
            &sample_index.all_contigs(),
            &mpb,
        )?;

        let batch_size =
            self.batch_size.as_ref().map(|x| *x).unwrap_or_else(|| {
                (self.threads as f32 * 1.5f32).floor() as usize
            });

        if self.is_single_site() {
            info!("running single-site analysis");
            let linear_transitions = if self.fine_grained {
                false
            } else {
                !self.log_transition_decay
            };
            return SingleSiteDmrAnalysis::new(
                sample_index,
                genome_positions,
                self.cap_coverages,
                self.control_bed_methyl.len(),
                self.exp_bed_methyl.len(),
                batch_size,
                self.interval_size,
                self.prior.as_ref(),
                self.max_coverages.as_ref(),
                self.delta,
                self.n_sample_records,
                self.header,
                self.segmentation_fp.as_ref(),
                &mpb,
            )?
            .run(
                mpb,
                pool,
                self.max_gap_size,
                self.dmr_prior,
                self.diff_stay,
                self.significance_factor,
                self.decay_distance,
                linear_transitions,
                writer,
            );
        }

        let sample_index = Arc::new(sample_index);
        let genome_positions = Arc::new(genome_positions);

        let regions_of_interest =
            if let Some(roi_bed) = self.regions_bed.as_ref() {
                let rois = parse_roi_bed(roi_bed).with_context(|| {
                    format!("failed to parse supplied regions at {roi_bed:?}")
                })?;
                info!("loaded {} regions", rois.len());
                rois
            } else {
                unreachable!(
                    "regions should always be available unless we're doing \
                     single-site analysis"
                )
            };

        info!("loading {batch_size} regions at a time");

        let pb = mpb.add(get_master_progress_bar(regions_of_interest.len()));
        pb.set_message(format!("regions processed"));
        let failures = mpb.add(get_ticker());
        failures.set_message(format!("regions failed to process"));

        let dmr_interval_iter = RoiIter::new(
            0,
            1,
            "a",
            "b",
            sample_index.clone(),
            regions_of_interest,
            batch_size,
            failures.clone(),
            self.handle_missing,
            genome_positions.clone(),
        )?;

        let success_count = run_pairwise_dmr(
            dmr_interval_iter,
            sample_index.clone(),
            pool,
            writer,
            pb,
            self.header,
            "a",
            "b",
            failures.clone(),
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
#[command(arg_required_else_help = true)]
pub struct MultiSampleDmr {
    /// Two or more named samples to compare. Two arguments are required <path>
    /// <name>. This option should be repeated at least two times.
    #[arg(short = 's', long = "sample", num_args = 2)]
    samples: Vec<String>,
    /// Optional, paths to tabix indices associated with named samples. Two
    /// arguments are required <path> <name> where <name> corresponds to
    /// the name of the sample given to the -s/--sample argument.
    #[arg(short = 'i', long = "index", num_args = 2)]
    indices: Vec<String>,
    /// BED file of regions over which to compare methylation levels. Should be
    /// tab-separated (spaces allowed in the "name" column). Requires
    /// chrom, chromStart and chromEnd. The Name column is optional. Strand
    /// is currently ignored.
    #[arg(long, short = 'r', alias = "regions")]
    regions_bed: PathBuf,
    /// Include header in output
    #[arg(long, default_value_t = false)]
    header: bool,
    /// Directory to place output DMR results in BED format.
    #[arg(short = 'o', long)]
    out_dir: PathBuf,
    /// Prefix files in directory with this label
    #[arg(short = 'p', long)]
    prefix: Option<String>,
    /// Path to reference fasta for the pileup.
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Bases to use to calculate DMR, may be multiple. For example, to
    /// calculate differentially methylated regions using only cytosine
    /// modifications use --base C.
    #[arg(short, long="base", alias = "modified-bases", action=clap::ArgAction::Append)]
    modified_bases: Vec<char>,
    /// Extra assignments of modification codes to their respective primary
    /// bases. In general, modkit dmr will use the SAM specification to
    /// know which modification codes are appropriate to use for a given
    /// primary base. For example "h" is the code for 5hmC, so is appropriate
    /// for cytosine bases, but not adenine bases. However, if your
    /// bedMethyl file contains custom codes or codes that are not part of
    /// the specification, you can specify which primary base they
    /// belong to here with --assign-code x:C meaning associate modification
    /// code "x" with cytosine (C) primary sequence bases. If a code is
    /// encountered that is not part of the specification, the bedMethyl
    /// record will not be used, this will be logged.
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,
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
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::warn)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
}

impl MultiSampleDmr {
    fn get_writer(
        &self,
        a_name: &str,
        b_name: &str,
    ) -> anyhow::Result<Box<BufWriter<File>>> {
        let fp = if let Some(p) = self.prefix.as_ref() {
            self.out_dir.join(format!("{}_{}_{}.bed", p, a_name, b_name))
        } else {
            self.out_dir.join(format!("{}_{}.bed", a_name, b_name))
        };
        if fp.exists() && !self.force {
            bail!(
                "refusing to overwrite {:?}",
                fp.to_str().unwrap_or("failed decode")
            )
        } else {
            let fh = File::create(fp.clone()).with_context(|| {
                format!("failed to make output file at {fp:?}")
            })?;
            Ok(Box::new(BufWriter::new(fh)))
        }
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if !self.out_dir.exists() {
            info!("creating directory at {:?}", &self.out_dir);
            std::fs::create_dir_all(&self.out_dir)?;
        }
        let code_lookup = PairwiseDmr::validate_modified_bases(
            &self.modified_bases,
            self.mod_code_assignments.as_ref(),
        )?;
        let index_fps = self
            .indices
            .chunks(2)
            .filter_map(|raw| {
                if raw.len() != 2 {
                    error!(
                        "illegal index pair {:?}, should be length 2 of the \
                         form <path> <name>",
                        raw
                    );
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

        let handlers = self
            .samples
            .chunks(2)
            .filter_map(|raw| {
                if raw.len() != 2 {
                    error!(
                        "illegal sample pair {:?}, should be length 2 of the \
                         form <path> <name>",
                        raw
                    );
                    None
                } else {
                    let fp = Path::new(raw[0].as_str()).to_path_buf();
                    let name = raw[1].to_string();
                    if fp.exists() {
                        let specified_index_fp = index_fps.get(&name);
                        match IndexHandler::new_infer_index_filepath(
                            &fp,
                            specified_index_fp,
                        ) {
                            Ok(handler) => Some((name, handler)),
                            Err(e) => {
                                error!("failed to load {name}, {e}");
                                None
                            }
                        }
                    } else {
                        error!("bedMethyl for {name} at {} not found", &raw[0]);
                        None
                    }
                }
            })
            .collect::<Vec<(String, IndexHandler)>>();

        let mpb = MultiProgress::new();

        let motifs = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<anyhow::Result<Vec<DnaBase>>>()
            .context("failed to parse modified base")?;

        let (names, handlers) = handlers.into_iter().fold(
            (Vec::new(), Vec::new()),
            |(mut names, mut handlers), (name, handler)| {
                names.push(name);
                handlers.push(handler);
                (names, handlers)
            },
        );

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
        );

        let genome_positions = GenomePositions::new_from_sequences(
            &motifs,
            &self.reference_fasta,
            self.mask,
            &sample_index.all_contigs(),
            &mpb,
        )?;

        let regions_of_interest = parse_roi_bed(&self.regions_bed)?;

        let sample_index = Arc::new(sample_index);
        let genome_positions = Arc::new(genome_positions);

        info!("loaded {} regions", regions_of_interest.len());

        let chunk_size = (self.threads as f32 * 1.5f32).floor() as usize;
        info!("processing {chunk_size} regions concurrently");

        let sample_pb =
            mpb.add(get_master_progress_bar(sample_index.num_combinations()?));

        for pair in sample_index
            .samples()
            .iter()
            .combinations(2)
            .progress_with(sample_pb.clone())
        {
            let a_index = *pair[0];
            let b_index = *pair[1];
            let a_name = &names[a_index];
            let b_name = &names[b_index];
            sample_pb
                .set_message(format!("comparing {} and {}", a_name, b_name));
            let pb =
                mpb.add(get_subroutine_progress_bar(regions_of_interest.len()));
            pb.set_message("regions processed");
            let failures = mpb.add(get_ticker());
            failures.set_message("regions failed to process");

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()?;

            match RoiIter::new(
                a_index,
                b_index,
                a_name,
                b_name,
                sample_index.clone(),
                // todo remove this clone.. but at most there will be 2 copies
                // around
                regions_of_interest.clone(),
                chunk_size,
                failures.clone(),
                self.handle_missing,
                genome_positions.clone(),
            ) {
                Ok(dmr_interval_iter) => {
                    let writer = self.get_writer(a_name, b_name)?;
                    let success_count = run_pairwise_dmr(
                        dmr_interval_iter,
                        sample_index.clone(),
                        pool,
                        writer,
                        pb,
                        self.header,
                        a_name,
                        b_name,
                        failures.clone(),
                    )?;
                    debug!(
                        "{} regions processed successfully and {} regions \
                         failed for pair {} {}",
                        success_count,
                        failures.position(),
                        &a_index,
                        &b_index,
                    );
                }
                Err(e) => {
                    error!(
                        "pair {} {} failed to process, {}",
                        &a_index,
                        &b_index,
                        e.to_string()
                    );
                    if self.handle_missing == HandleMissing::fail {
                        return Err(e);
                    }
                }
            }
        }

        Ok(())
    }
}
