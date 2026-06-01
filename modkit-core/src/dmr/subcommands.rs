use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::Neg;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::{Args, Subcommand};
use indexmap::IndexMap;
use indicatif::{MultiProgress, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use prettytable::row;
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::dmr::isoform::{
    build_gene_common_coords, build_gene_models, group_bedmethyl_records_by_tx,
    parse_gtf, plot_isoform_metrics, run_isoform_dmr_on_gene, volcano_svg,
    DmrMetric, GeneCommonCoord, GeneDmrScore, GeneIsoformDmr,
    GeneIsoformDmrRecord, GeneIsoformDmrScore, GeneIsoformDmrWorker, GeneTxDmr,
    GeneTxDmrWorker, GtfGene, GtfId, GtfTranscript,
    PairTranscriptomeBedmethylHandler, PlottingArgs,
    TranscriptBedmethylHandler, TranscriptModel,
};
use crate::dmr::pairwise::run_pairwise_dmr;
use crate::dmr::single_site::SingleSiteDmrAnalysis;
use crate::dmr::tabix::MultiSampleIndex;
use crate::dmr::util::{parse_roi_bed, HandleMissing, RoiIter};
use crate::errs::MkResult;
use crate::genome_positions::GenomePositions;
use crate::mod_base_code::{DnaBase, ModCodeRepr, MOD_CODE_TO_DNA_BASE};
use crate::monoid::Moniod;
use crate::tabix::BedMethylTbxIndex;
use crate::util::{
    create_out_directory, format_errors_table, get_master_progress_bar,
    get_master_progress_bar_fancy, get_subroutine_progress_bar, get_ticker,
};
use modkit_logging::init_logging;

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
    /// Use a transcriptome-aligned bedMethyl table and GTF file to compare
    /// methylation across gene transcript isoforms. Run analysis on the entire
    /// transcriptome or on a single gene. Optionally also make plots for
    /// single genes.
    Isoform(EntryDmrIsoform),
    /// Uses two transcriptome-aligned bedmMethyl tables to compare methylation
    /// for two conditions at transcriptome sites aggregated at the gene
    /// level. This command will combine the methylation counts across
    /// isoforms into a single summary statistic per-site for the gene,
    /// then compare those counts between conditions at the gene level and
    /// report results on genomic coordinates. The transcript models are
    /// provided by a GTF file. Optionally make a volcano plot of the top-k
    /// differentially methylated sites per gene.
    #[command(
        name = "compare-tx-sites",
        alias = "compare-transcriptome-sites"
    )]
    GenePair(EntryGeneTx),
}

impl BedMethylDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Self::Pair(x) => x.run(),
            Self::Multi(x) => x.run(),
            Self::Isoform(x) => x.run(),
            Self::GenePair(x) => x.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct PairwiseDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There
    /// should be a tabix index with the same name and .tbi next to this
    /// file or the --index-a option must be provided.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 'a')]
    control_bed_methyl: Vec<PathBuf>,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample.
    /// There should be a tabix index with the same name and .tbi next to
    /// this file or the --index-b option must be provided.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 'b')]
    exp_bed_methyl: Vec<PathBuf>,
    /// Path to file to direct output, optional, no argument will direct output
    /// to stdout.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'o', long)]
    out_path: Option<String>,
    /// Include header in output
    #[clap(help_heading = "Output Options")]
    #[arg(long, alias = "with-header", default_value_t = false)]
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
    #[clap(help_heading = "Segmentation Options")]
    #[arg(long = "segment", conflicts_with = "regions_bed")]
    segmentation_fp: Option<PathBuf>,

    /// Maximum number of base pairs between modified bases for them to be
    /// segmented together.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(long, requires = "segmentation_fp", default_value_t = 5000)]
    max_gap_size: u64,
    /// Prior probability of a differentially methylated position
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.1,
        hide_short_help = true
    )]
    dmr_prior: f64,
    /// Maximum probability of continuing a differentially methylated block,
    /// decay will be dynamic based on proximity to the next position.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.9,
        hide_short_help = true
    )]
    diff_stay: f64,
    /// Significance factor, effective p-value necessary to favor the
    /// "Different" state.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.01,
        hide_short_help = true
    )]
    significance_factor: f64,
    /// Use logarithmic decay for "Different" stay probability
    #[clap(help_heading = "Segmentation Options")]
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
    #[clap(help_heading = "Segmentation Options")]
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
    #[clap(help_heading = "Segmentation Options")]
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
    #[clap(help_heading = "Sample Options")]
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
    #[clap(help_heading = "Sample Options")]
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,
    /// Calculate differential methylation for this modification code in
    /// isolation.
    #[clap(help_heading = "Sample Options")]
    #[arg(long = "single-code")]
    single_modification_code: Option<String>,

    /// Log out which sequences are in common between the samples and the
    /// reference FASTA, useful for debugging
    #[clap(help_heading = "Logging Options")]
    #[arg(
        long = "careful",
        requires = "log_filepath",
        hide_short_help = true,
        default_value_t = false
    )]
    careful: bool,
    /// File to write logs to, it's recommended to use this option.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of threads to use when for decompression.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: usize,
    /// Control the  batch size. The batch size is the number of regions to
    /// load at a time. Each region will be processed concurrently. Loading
    /// more regions at a time will decrease IO to load data, but will use
    /// more memory. Default will be 50% more than the number of
    /// threads assigned.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, alias = "batch")]
    batch_size: Option<usize>,
    /// Respect soft masking in the reference FASTA.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[clap(help_heading = "Logging Options")]
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::quiet)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
    /// Prior distribution for estimating MAP-based p-value. Should be two
    /// arguments for alpha and beta (e.g. 1.0 1.0). See
    /// `dmr_scoring_details.md` for additional details on how the metric
    /// is calculated.

    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        num_args = 2,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    prior: Option<Vec<f64>>,
    /// Consider only effect sizes greater than this when calculating the
    /// MAP-based p-value.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        default_value_t = 0.05,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    delta: f64,
    /// Sample this many reads when estimating the max coverage thresholds.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        short='N',
        default_value_t = 10_042,
        conflicts_with_all = ["max_coverages", "regions_bed"],
    )]
    n_sample_records: usize,
    /// Max coverages to enforce when calculating estimated MAP-based p-value.
    #[clap(help_heading = "Single-site Options")]
    #[arg(long, num_args = 2, conflicts_with = "regions_bed")]
    max_coverages: Option<Vec<usize>>,
    /// When using replicates, cap coverage to be equal to the maximum coverage
    /// for a single sample. For example, if there are 3 replicates with
    /// max_coverage of 30, the total coverage would normally be 90. Using
    /// --cap-coverages will down sample the data to 30X.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        conflicts_with = "regions_bed",
        default_value_t = false,
        hide_short_help = true
    )]
    cap_coverages: bool,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[clap(help_heading = "Compute Options")]
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
                                .and_then(|raw| {
                                    DnaBase::parse(raw).map_err(|e| e.into())
                                })?;
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
            .collect::<MkResult<Vec<DnaBase>>>()?;

        if self.regions_bed.is_some()
            & (self.control_bed_methyl.len() > 1
                || self.exp_bed_methyl.len() > 1)
        {
            info!(
                "multiple samples will be combined and DMR will be performed \
                 over regions"
            );
        }

        let a_handlers = self
            .control_bed_methyl
            .iter()
            .map(|fp| BedMethylTbxIndex::from_path(fp))
            .collect::<anyhow::Result<Vec<BedMethylTbxIndex>>>()?;
        let b_handlers = self
            .exp_bed_methyl
            .iter()
            .map(|fp| BedMethylTbxIndex::from_path(fp))
            .collect::<anyhow::Result<Vec<BedMethylTbxIndex>>>()?;
        let handlers = a_handlers
            .into_iter()
            .chain(b_handlers)
            .collect::<Vec<BedMethylTbxIndex>>();

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
            self.io_threads,
        );
        let total = self.control_bed_methyl.len() + self.exp_bed_methyl.len();
        let control_idxs =
            (0..self.control_bed_methyl.len()).collect::<Vec<usize>>();
        let exp_idxs =
            (self.control_bed_methyl.len()..total).collect::<Vec<usize>>();

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
        let mut tab = prettytable::Table::new();
        tab.set_format(
            *prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE,
        );
        tab.set_titles(row!["contig", "a_contains", "b_contains", "both"]);
        let mut common_contigs = 0usize;
        for (name, _) in genome_positions.contig_sizes() {
            let a_contains =
                control_idxs.iter().any(|i| sample_index.has_contig(*i, name));
            let b_contains =
                exp_idxs.iter().any(|i| sample_index.has_contig(*i, name));
            tab.add_row(row![
                name,
                a_contains,
                b_contains,
                a_contains && b_contains
            ]);
            if a_contains && b_contains {
                common_contigs += 1;
            }
        }
        if self.careful || common_contigs == 0 {
            debug!("contig breakdown:\n{tab}");
        }
        mpb.suspend(|| {
            info!(
                "{common_contigs} common sequence(s) between FASTA and both \
                 samples"
            );
        });

        let batch_size =
            self.batch_size.as_ref().map(|x| *x).unwrap_or_else(|| {
                (self.threads as f32 * 1.5f32).floor() as usize
            });

        let single_code_op = self
            .single_modification_code
            .as_ref()
            .map(|x| ModCodeRepr::parse(&x))
            .transpose()?;
        if let Some(code) = single_code_op.as_ref() {
            info!("calculating differntial methylation for {code} in isolation")
        }

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
                single_code_op,
                mpb.clone(),
                &pool,
            )?
            .run(
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
        pb.set_message("regions processed");
        let failures = mpb.add(get_ticker());
        failures.set_message("regions failed to process");
        let batch_failures = mpb.add(get_ticker());
        batch_failures.set_message("failed batches");

        let dmr_interval_iter = RoiIter::new(
            control_idxs.as_slice(),
            exp_idxs.as_slice(),
            "a",
            "b",
            sample_index.clone(),
            regions_of_interest,
            batch_size,
            self.handle_missing,
            genome_positions.clone(),
            &mpb,
        )?;

        let (success_count, region_errors) = run_pairwise_dmr(
            dmr_interval_iter,
            sample_index.clone(),
            pool,
            writer,
            pb,
            self.header,
            "a",
            "b",
            single_code_op,
            failures.clone(),
            batch_failures.clone(),
            mpb.clone(),
        )?;

        mpb.suspend(|| {
            info!(
                "{} regions processed successfully and {} regions failed",
                success_count,
                failures.position()
            );
            if !region_errors.is_empty() {
                let tab = format_errors_table(&region_errors);
                error!("region errors:\n{tab}");
            }
        });

        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct MultiSampleDmr {
    /// Two or more named samples to compare. Two arguments are required <path>
    /// <name>. This option should be repeated at least two times. When two
    /// samples have the same name, they will be combined.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 's', long = "sample", num_args = 2)]
    samples: Vec<String>,
    /// BED file of regions over which to compare methylation levels. Should be
    /// tab-separated (spaces allowed in the "name" column). Requires
    /// chrom, chromStart and chromEnd. The Name column is optional. Strand
    /// is currently ignored.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'r', alias = "regions")]
    regions_bed: PathBuf,
    /// Include header in output
    #[clap(help_heading = "Output Options")]
    #[arg(long, alias = "with-header", default_value_t = false)]
    header: bool,
    /// Directory to place output DMR results in BED format.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'o', long)]
    out_dir: PathBuf,
    /// Prefix files in directory with this label
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'p', long)]
    prefix: Option<String>,
    /// Path to reference fasta for the pileup.
    #[clap(help_heading = "Sample Options")]
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Bases to use to calculate DMR, may be multiple. For example, to
    /// calculate differentially methylated regions using only cytosine
    /// modifications use --base C.
    #[clap(help_heading = "Sample Options")]
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
    #[clap(help_heading = "Sample Options")]
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,
    /// Calculate differential methylation for this modification code in
    /// isolation.
    #[clap(help_heading = "Sample Options")]
    #[arg(long = "single-code")]
    single_modification_code: Option<String>,

    /// File to write logs to, it's recommended to use this option.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of threads to use when for decompression.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: usize,
    /// Respect soft masking in the reference FASTA.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[clap(help_heading = "Logging Options")]
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::quiet)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
    /// Only performs comparisons between this sample and all other samples.
    #[clap(help_heading = "Sample Options")]
    #[arg(long)]
    ref_sample: Option<String>,
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

        let handlers = self
            .samples
            .chunks(2)
            .enumerate()
            .filter_map(|(i, raw)| {
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
                        match BedMethylTbxIndex::from_path(&fp) {
                            Ok(handler) => Some((i, name, handler)),
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
            .collect::<Vec<(usize, String, BedMethylTbxIndex)>>();

        let mpb = MultiProgress::new();

        let motifs = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<MkResult<Vec<DnaBase>>>()
            .context("failed to parse modified base")?;

        let (names, handlers) = handlers.into_iter().fold(
            (HashMap::new(), Vec::new()),
            |(mut names, mut handlers), (sample_id, name, handler)| {
                names.entry(name).or_insert_with(Vec::new).push(sample_id);
                handlers.push(handler);
                (names, handlers)
            },
        );
        for (name, ids) in &names {
            if ids.len() > 1 {
                info!(
                    "sample {name} has {} replicates, they will be combined",
                    ids.len()
                );
            }
        }

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
            self.io_threads,
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

        let single_code_op = self
            .single_modification_code
            .as_ref()
            .map(|x| ModCodeRepr::parse(&x))
            .transpose()?;
        if let Some(code) = single_code_op.as_ref() {
            info!("calculating differntial methylation for {code} in isolation")
        }

        // Generate sample pairs based on whether a reference sample is specified
        let sample_pairs: Vec<(&String, &String)> = if let Some(ref_name) =
            &self.ref_sample
        {
            // Validate that the reference sample exists
            if !names.contains_key(ref_name) {
                bail!(
                    "reference sample '{}' not found in provided samples. Available samples: {}",
                    ref_name,
                    names.keys().sorted().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
                );
            }
            info!("using reference sample: {}", ref_name);

            // Compare all samples against the reference
            names
                .keys()
                .filter(|name| *name != ref_name)
                .map(|name| (ref_name, name))
                .collect()
        } else {
            // All pairwise comparisons
            let samples = names.keys().sorted().collect::<Vec<&String>>();
            samples
                .into_iter()
                .combinations(2)
                .map(|pair| (pair[0], pair[1]))
                .collect()
        };

        let sample_pb =
            mpb.add(get_master_progress_bar(sample_pairs.len() as u64));

        for (a_name, b_name) in
            sample_pairs.into_iter().progress_with(sample_pb.clone())
        {
            let a_idxs = names.get(a_name).unwrap();
            let b_idxs = names.get(b_name).unwrap();

            sample_pb
                .set_message(format!("comparing {} and {}", a_name, b_name));
            let pb =
                mpb.add(get_subroutine_progress_bar(regions_of_interest.len()));
            pb.set_message("regions processed");
            let failures = mpb.add(get_ticker());
            failures.set_message("regions failed to process");
            let batch_failures = mpb.add(get_ticker());
            batch_failures.set_message("failed batches");

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()?;

            debug!("running {a_name} as control and {b_name} as experiment");
            let mut all_region_errors = FxHashMap::default();
            match RoiIter::new(
                a_idxs,
                b_idxs,
                a_name,
                b_name,
                sample_index.clone(),
                regions_of_interest.clone(),
                chunk_size,
                self.handle_missing,
                genome_positions.clone(),
                &mpb,
            ) {
                Ok(dmr_interval_iter) => {
                    let writer = self.get_writer(a_name, b_name)?;
                    let (success_count, region_errors) = run_pairwise_dmr(
                        dmr_interval_iter,
                        sample_index.clone(),
                        pool,
                        writer,
                        pb,
                        self.header,
                        a_name,
                        b_name,
                        single_code_op,
                        failures.clone(),
                        batch_failures.clone(),
                        mpb.clone(),
                    )?;
                    mpb.suspend(|| {
                        info!(
                            "{} regions processed successfully and {} regions \
                             failed for pair {} {}",
                            success_count,
                            failures.position(),
                            &a_name,
                            &b_name,
                        );
                        if !region_errors.is_empty() {
                            let tab = format_errors_table(&region_errors);
                            error!("region errors:\n{tab}");
                            all_region_errors.op_mut(region_errors);
                        }
                    });
                }
                Err(e) => {
                    mpb.suspend(|| {
                        error!(
                            "pair {} {} failed to process, {}",
                            &a_name,
                            &b_name,
                            e.to_string()
                        );
                    });
                    if self.handle_missing == HandleMissing::fail {
                        return Err(e);
                    }
                }
            }
        }

        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryDmrIsoform {
    /// Input bedmethyl table. Expected to be bgzip-compressed and have a .tbi
    /// tabix index. This bedMethyl table must be generated from a
    /// transcriptome-aligned modBAM, NOT a genome-aligned modBAM. Expected
    /// that this table was generated by pileup with `--modified-bases`
    /// argument.
    input_bedmethyl: PathBuf,
    /// Path to output BED file or '-' for standard output.
    output_table: String,
    /// Path to GTF file to use for transcript models. May be compressed.
    #[arg(long)]
    gtf: PathBuf,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[arg(long, default_value_t = 1)]
    min_valid_coverage: u64,
    /// Calculate differential methylation for this modification code in
    /// isolation. Unmodified and other modification counts will be combined.
    #[arg(long = "single-code")]
    single_modification_code: Option<String>,
    /// Include per-modification combined proportions per position and
    /// per-isoform proportions in output. Adding this flag will make the
    /// output file substantially larger.
    #[arg(long = "full", default_value_t = false)]
    emit_full_results: bool,
    /// Path to optional debug log (recommended).
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Only process this Gene Id
    #[arg(long, short = 'g')]
    gene_id: Option<String>,
    /// Only process this Gene Name
    #[arg(long, short = 'n', conflicts_with = "gene_id")]
    gene_name: Option<String>,
    /// Don't show the progress bar.
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Number of concurrent processors to use. Only used when running on all
    /// genes.
    #[arg(long, short='t', default_value_t = 4, conflicts_with_all = ["gene_id", "gene_name"])]
    threads: usize,
    /// Generate an SVG plot displaying the methylation positions on each
    /// isoform (filtered by --min-score or --max-pval) and the common
    /// exons of the gene. Also displays negative log of the p-value of the
    /// DMR metric. See the online documentation for an example.
    #[arg(long)]
    plot: Option<PathBuf>,
    #[command(flatten)]
    svg_args: PlottingArgs,
}

impl EntryDmrIsoform {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());
        if let Some(p) = self.plot.as_ref() {
            if p.is_file() {
                bail!("path for --plot should be a directory, got a file {p:?}")
            }
        }
        if self.plot.is_some()
            && (self.gene_id.is_none() && self.gene_name.is_none())
        {
            bail!("--plot requires --gene-id or --gene-name")
        }
        if self.plot.is_some() && self.single_modification_code.is_none() {
            info!(
                "plotting without --single-code will mean that the DMR metric \
                 track represents changes in proportion of all modifications, \
                 not just the one plotted"
            );
        }
        let multi_progress = MultiProgress::new();
        if self.suppress_progress {
            multi_progress
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let single_mod_code = self
            .single_modification_code
            .as_ref()
            .map(|x| ModCodeRepr::parse(x))
            .transpose()?;

        let sorted_transcript_models = parse_gtf(&self.gtf, &multi_progress)
            .context("failed to read GTF")?;
        multi_progress.suspend(|| {
            info!(
                "parsed {} transcript models",
                sorted_transcript_models.len()
            );
        });
        let sorted_by_gene_common_coordinates =
            build_gene_common_coords(&sorted_transcript_models)?;
        let transcript_models =
            sorted_transcript_models.into_iter().collect::<FxHashMap<_, _>>();
        let bedmethyl_handler = TranscriptBedmethylHandler::from_path(
            &self.input_bedmethyl,
            &transcript_models,
            self.min_valid_coverage,
            multi_progress.clone(),
        )
        .context("failed to build bedmethyl handler")?;

        let genes_in_common = bedmethyl_handler
            .gene_id_to_transcript_ids
            .iter()
            .filter(|(_, x)| !x.is_empty())
            .map(|(x, _)| {
                if sorted_by_gene_common_coordinates.contains_key(x) {
                    1usize
                } else {
                    0usize
                }
            })
            .sum::<usize>();

        multi_progress.suspend(|| {
            info!(
                "parsed bedmethyl header and found transcripts corresponding \
                 to {} genes, {genes_in_common} of those are in common with \
                 GTF ({}%)",
                bedmethyl_handler.gene_id_to_transcript_ids.len(),
                genes_in_common as f32
                    / bedmethyl_handler.gene_id_to_transcript_ids.len() as f32
                    * 100f32
            );
        });

        let requested_gene_id =
            match (self.gene_id.as_ref(), self.gene_name.as_ref()) {
                (None, None) => None,
                (None, Some(name)) => {
                    multi_progress
                        .suspend(|| info!("searching for gene name {name}"));
                    sorted_by_gene_common_coordinates.values().find_map(|gc| {
                        if let Some(gn) = &gc.gene_name {
                            if gn == name {
                                Some(gc.gene_id.to_owned())
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    })
                }
                (Some(gene_id), None) => Some(GtfGene::from_str(gene_id)?),
                (Some(_), Some(_)) => unreachable!(),
            };

        if self.gene_name.is_some() && requested_gene_id.is_none() {
            bail!("failed to find requested gene name in GTF");
        }

        if let Some(requested_gene_id) = requested_gene_id.as_ref() {
            let Some(gene) =
                sorted_by_gene_common_coordinates.get(requested_gene_id)
            else {
                bail!("didn't find {requested_gene_id} in GTF")
            };
            self.run_single_gene(
                gene,
                transcript_models,
                bedmethyl_handler,
                single_mod_code,
                multi_progress,
            )?;
        } else {
            self.run_on_all_genes(
                transcript_models,
                sorted_by_gene_common_coordinates,
                bedmethyl_handler,
                single_mod_code,
                multi_progress,
            )?;
        }
        Ok(())
    }

    fn run_single_gene(
        &self,
        gene: &GeneCommonCoord,
        transcript_models: FxHashMap<GtfTranscript, TranscriptModel>,
        mut bedmethyl_handler: TranscriptBedmethylHandler,
        single_mod_code: Option<ModCodeRepr>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<()> {
        let mut writer = self.get_writer(&multi_progress)?;
        let Some(records) = run_isoform_dmr_on_gene::<6>(
            &transcript_models,
            gene,
            &mut bedmethyl_handler,
            single_mod_code,
            self.emit_full_results,
        )?
        else {
            bail!("didn't find {} in bedMethyl", &gene.gene_id)
        };

        for row in records.iter().map(|rec| {
            rec.to_row(
                &gene.chrom,
                &gene.gene_id,
                gene.gene_name.as_ref(),
                self.emit_full_results,
            )
        }) {
            writer.write(row.as_bytes())?;
        }
        if let Some(plot_dir) = self.plot.as_ref() {
            if !plot_dir.exists() {
                multi_progress.suspend(|| {
                    info!("creating directory for plotting output")
                });
                std::fs::create_dir_all(plot_dir)
                    .context("failed to create plotting output directory")?;
            }
            let gene_models = build_gene_models(&transcript_models)?;
            let Some(gene_model) = gene_models.get(&gene.gene_id) else {
                bail!("didn't find gene-id {}", gene.gene_id)
            };
            let Some(bedmethyl_records) = bedmethyl_handler
                .get_read_bedmethyl_for_gene(&gene.gene_id, None)?
            else {
                bail!("didn't find any bedmethyl records for {}", gene.gene_id)
            };
            if let Some(c) = single_mod_code {
                let dmr_metrics = records
                    .iter()
                    .map(|rec| DmrMetric {
                        start: rec.start,
                        value: rec.score(true).ln().neg(),
                    })
                    .collect::<Vec<DmrMetric>>();
                let meth_by_tx =
                    group_bedmethyl_records_by_tx(bedmethyl_records);
                plot_isoform_metrics(
                    meth_by_tx,
                    &dmr_metrics,
                    c,
                    plot_dir,
                    gene_model,
                    gene.gene_name.as_ref(),
                    &transcript_models,
                    &self.svg_args,
                    &multi_progress,
                )?;
            } else {
                let mod_codes = bedmethyl_records
                    .iter()
                    .map(|x| x.raw_mod_code)
                    .unique()
                    .collect::<Vec<ModCodeRepr>>();
                if mod_codes.len() > 1 {
                    multi_progress.suspend(|| {
                        info!(
                            "found {} modification codes, {}, plotting them \
                             all separately",
                            mod_codes.len(),
                            mod_codes.iter().map(|x| x.to_string()).join(",")
                        );
                    });
                }
                for mod_code in mod_codes {
                    let dmr_metrics = records
                        .iter()
                        .filter(|rec| rec.has_code(mod_code))
                        .map(|rec| DmrMetric {
                            start: rec.start,
                            value: rec.score(true).ln().neg(),
                        })
                        .collect::<Vec<DmrMetric>>();
                    let bedmethyl_records_for_code = bedmethyl_records
                        .iter()
                        .filter(|rec| rec.raw_mod_code == mod_code)
                        .cloned()
                        .collect::<Vec<BedMethylLine>>();
                    let meth_by_tx = group_bedmethyl_records_by_tx(
                        bedmethyl_records_for_code,
                    );
                    plot_isoform_metrics(
                        meth_by_tx,
                        &dmr_metrics,
                        mod_code,
                        plot_dir,
                        gene_model,
                        gene.gene_name.as_ref(),
                        &transcript_models,
                        &self.svg_args,
                        &multi_progress,
                    )?;
                }
            }
        }
        Ok(())
    }

    fn run_on_all_genes(
        &self,
        transcript_models: FxHashMap<GtfTranscript, TranscriptModel>,
        sorted_gene_common_coords: IndexMap<GtfGene, GeneCommonCoord>,
        bedmethyl_handler: TranscriptBedmethylHandler,
        single_mod_code: Option<ModCodeRepr>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<()> {
        multi_progress.suspend(|| {
            info!(
                "processing all genes, found {} in GTF",
                sorted_gene_common_coords.len()
            );
        });
        let mut writer = self.get_writer(&multi_progress)?;
        let transcript_models = Arc::new(transcript_models);

        let (empties_tx, empties_rx) = crossbeam_channel::unbounded();
        let (jobs_tx, jobs_rx) = crossbeam_channel::bounded(self.threads * 2);
        let (results_tx, results_rx) = crossbeam_channel::unbounded();
        let (records_tx, records_rx) =
            crossbeam_channel::bounded(self.threads * 2);

        let pb = multi_progress.add(get_master_progress_bar_fancy(
            sorted_gene_common_coords.len(),
        ));
        pb.set_message("genes processed");
        let records_counter = multi_progress.add(get_ticker());
        records_counter.set_message("records written");
        let mut workers = Vec::with_capacity(self.threads);
        for _ in 0..self.threads {
            workers.push(GeneIsoformDmrWorker::new(
                &self.input_bedmethyl,
                transcript_models.clone(),
                multi_progress.clone(),
                bedmethyl_handler.gene_id_to_transcript_ids.clone(),
                self.min_valid_coverage,
                single_mod_code,
                self.emit_full_results,
            )?);
            empties_tx
                .send(GeneIsoformDmr::empty())
                .context("should stage memory 1")?;
            empties_tx
                .send(GeneIsoformDmr::empty())
                .context("should stage memory 2")?;
        }
        multi_progress
            .suspend(|| info!("workers staged, starting processing.."));

        let source_thread = std::thread::spawn({
            let results_handle = results_tx.clone();
            move || {
                let get_mem = || -> Result<GeneIsoformDmr, ()> {
                    match empties_rx.recv() {
                        Ok(mem) => Ok(mem),
                        Err(_) => Err(()),
                    }
                };
                let mut seq = 0usize;
                for (_gene_id, gene) in sorted_gene_common_coords {
                    match get_mem() {
                        Ok(mem) => {
                            if jobs_tx.send((seq, gene, mem)).is_err() {
                                break;
                            }
                            seq = seq.wrapping_add(1)
                        }
                        Err(_) => break,
                    }
                }
                drop(jobs_tx);
                drop(results_handle);
            }
        });

        let mut handles = Vec::with_capacity(workers.len());
        for mut worker_state in workers {
            let jobs_rx = jobs_rx.clone();
            let results_tx = results_tx.clone();
            handles.push(std::thread::spawn(move || {
                while let Ok((seq, gene, mem)) = jobs_rx.recv() {
                    let out = worker_state.process_gene::<6>(mem, gene);
                    if results_tx.send((seq, out)).is_err() {
                        break;
                    }
                }
            }))
        }
        drop(results_tx);
        drop(jobs_rx);

        let aggregator = std::thread::spawn(move || {
            let mut next_seq = 0usize;
            let mut buffer: BTreeMap<usize, MkResult<GeneIsoformDmr>> =
                BTreeMap::new();
            while let Ok((seq, out)) = results_rx.recv() {
                buffer.insert(seq, out);

                while let Some(v) = buffer.remove(&next_seq) {
                    if records_tx.send(v).is_err() {
                        return;
                    }
                    next_seq = next_seq.wrapping_add(1);
                }
            }
            drop(records_tx);
        });
        let mut errs = HashMap::new();
        for result in records_rx {
            match result {
                Ok(mut gene_isoform_dmr) => {
                    let records_written = gene_isoform_dmr
                        .write(&mut writer, self.emit_full_results)?;
                    gene_isoform_dmr.clear();
                    let _ = empties_tx.send(gene_isoform_dmr);
                    records_counter.inc(records_written as u64);
                    pb.inc(1);
                }
                Err(e) => {
                    let c = errs.entry(e.to_string()).or_insert(0usize);
                    *c = c.saturating_add(1);
                }
            }
        }
        source_thread.join().expect("source thread paniced");
        for (i, worker_thread) in handles.into_iter().enumerate() {
            worker_thread.join().expect(&format!("worker {i} paniced"));
        }
        aggregator.join().expect("aggregator theread paniced");

        multi_progress.suspend(|| {
            info!("finished, processed {} genes", pb.position());
        });

        Ok(())
    }

    fn get_writer(
        &self,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Box<dyn Write>> {
        let mut writer: Box<dyn Write> = match self.output_table.as_ref() {
            "-" | "stdout" => {
                multi_progress.suspend(|| info!("writing to stdout"));
                Box::new(BufWriter::new(std::io::stdout()))
            }
            p @ _ => {
                multi_progress.suspend(|| info!("writing to file at {p}"));
                let fh = File::create(Path::new(p))
                    .context("failed to open output file")?;
                Box::new(BufWriter::new(fh))
            }
        };
        writer.write(
            GeneIsoformDmrRecord::<GeneIsoformDmrScore>::header(
                self.emit_full_results,
            )
            .as_bytes(),
        )?;
        Ok(writer)
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryGeneTx {
    /// Input bedmethyl table for first condition (condition 'A').
    /// Expected to be bgzip-compressed and have a .tbi
    /// tabix index. This bedMethyl table must be generated from a
    /// transcriptome-aligned modBAM, NOT a genome-aligned modBAM. Expected
    /// that this table was generated by pileup with `--modified-bases`
    /// argument.
    #[arg(short = 'a')]
    cond_a: PathBuf,
    #[arg(short = 'b')]
    /// Input bedmethyl table for first condition (condition 'B').
    cond_b: PathBuf,
    /// Path to output BED file or '-' for stdout.
    #[arg(long = "out", short = 'o')]
    output_table: String,
    /// Include combined proportions and counts for each condition. This will
    /// increase the size of the output file.
    #[arg(long = "full", default_value_t = false)]
    emit_full_results: bool,
    /// Path to GTF file to use for transcript models. Can be compressed.
    #[arg(long)]
    gtf: PathBuf,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[arg(long, default_value_t = 1)]
    min_valid_coverage: u64,
    /// Calculate differential methylation for this modification code in
    /// isolation. Unmodified and other modification counts will be combined.
    #[arg(long = "single-code")]
    single_modification_code: Option<String>,
    /// Path to optional debug log (recommended).
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Don't show the progress bar.
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Number of concurrent processors to use.
    #[arg(long, short = 't', default_value_t = 4)]
    threads: usize,

    /// Path to file to create an SVG volcano plot where the x-axis is log2
    /// fold change and the y-axis is negative log p-value. See the online
    /// documentation for and example.
    #[arg(long, requires = "single_modification_code")]
    plot: Option<PathBuf>,
    /// For each gene, sort the positions in decreasing negative log p-value
    /// and take this many for plotting. Setting a large number will
    /// increase the size of the SVG plot.
    #[arg(long, requires = "plot", default_value_t = 1)]
    top_k: usize,
    /// Discard positions with an effect size less than this before considering
    /// them for plotting.
    #[arg(long, requires = "plot", default_value_t = 0.1f64)]
    min_effect_size: f64,
    /// Add a title to the plot
    #[arg(long, requires = "plot")]
    plot_title: Option<String>,
    /// Path to plain text file of a list of gene names, one per line, to mark
    /// in the volcano plot.
    #[arg(long, requires = "plot")]
    gene_labels: Option<PathBuf>,
}

impl EntryGeneTx {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());
        if self.top_k == 0 {
            bail!("--top-k must be > 0")
        }
        let multi_progress = MultiProgress::new();
        if self.suppress_progress {
            multi_progress
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let single_mod_code = self
            .single_modification_code
            .as_ref()
            .map(|x| ModCodeRepr::parse(x))
            .transpose()?;

        let sorted_transcript_models = parse_gtf(&self.gtf, &multi_progress)
            .context("failed to read GTF")?;
        multi_progress.suspend(|| {
            info!(
                "parsed {} transcript models",
                sorted_transcript_models.len()
            );
        });
        let mut sorted_by_gene_common_coordinates =
            build_gene_common_coords(&sorted_transcript_models)?;
        let n_genes = sorted_by_gene_common_coordinates.len();
        let transcript_models =
            sorted_transcript_models.into_iter().collect::<FxHashMap<_, _>>();

        let cond_a_bm_handler = TranscriptBedmethylHandler::from_path(
            &self.cond_a,
            &transcript_models,
            self.min_valid_coverage,
            multi_progress.clone(),
        )?;
        let cond_b_bm_handler = TranscriptBedmethylHandler::from_path(
            &self.cond_b,
            &transcript_models,
            self.min_valid_coverage,
            multi_progress.clone(),
        )?;
        let pair_handler = PairTranscriptomeBedmethylHandler::from_handlers(
            cond_a_bm_handler,
            cond_b_bm_handler,
            &mut sorted_by_gene_common_coordinates,
        )?;

        let mut writer = self.get_writer(single_mod_code, &multi_progress)?;
        let transcript_models = Arc::new(transcript_models);

        let (empties_tx, empties_rx) = crossbeam_channel::unbounded();
        let (jobs_tx, jobs_rx) = crossbeam_channel::bounded(self.threads * 2);
        let (results_tx, results_rx) = crossbeam_channel::unbounded();
        let (records_tx, records_rx) =
            crossbeam_channel::bounded(self.threads * 2);

        let pb = multi_progress.add(get_master_progress_bar_fancy(
            sorted_by_gene_common_coordinates.len(),
        ));
        pb.set_message("genes processed");
        let records_counter = multi_progress.add(get_ticker());
        records_counter.set_message("records written");

        let mut workers = Vec::with_capacity(self.threads);
        for _ in 0..self.threads {
            workers.push(GeneTxDmrWorker::new(
                &self.cond_a,
                &self.cond_b,
                &pair_handler,
                self.min_valid_coverage,
                transcript_models.clone(),
                single_mod_code,
                multi_progress.clone(),
            )?);
            empties_tx
                .send(GeneTxDmr::empty())
                .context("should stage memory 1")?;
            empties_tx
                .send(GeneTxDmr::empty())
                .context("should stage memory 2")?;
        }

        multi_progress
            .suspend(|| info!("workers staged, starting processing.."));

        let source_thread = std::thread::spawn({
            let results_handle = results_tx.clone();
            move || {
                let get_mem = || -> Result<GeneTxDmr, ()> {
                    match empties_rx.recv() {
                        Ok(mem) => Ok(mem),
                        Err(_) => Err(()),
                    }
                };
                let mut seq = 0usize;
                for (_gene_id, gene) in sorted_by_gene_common_coordinates {
                    match get_mem() {
                        Ok(mem) => {
                            if jobs_tx.send((seq, gene, mem)).is_err() {
                                break;
                            }
                            seq = seq.wrapping_add(1)
                        }
                        Err(_) => break,
                    }
                }
                drop(jobs_tx);
                drop(results_handle);
            }
        });
        let mut handles = Vec::with_capacity(workers.len());
        for mut worker_state in workers {
            let jobs_rx = jobs_rx.clone();
            let results_tx = results_tx.clone();
            handles.push(std::thread::spawn(move || {
                while let Ok((seq, gene, mem)) = jobs_rx.recv() {
                    let out = worker_state.process_gene::<6>(mem, gene);
                    if results_tx.send((seq, out)).is_err() {
                        break;
                    }
                }
            }))
        }
        drop(results_tx);
        drop(jobs_rx);

        let aggregator = std::thread::spawn(move || {
            let mut next_seq = 0usize;
            let mut buffer: BTreeMap<usize, MkResult<GeneTxDmr>> =
                BTreeMap::new();
            while let Ok((seq, out)) = results_rx.recv() {
                buffer.insert(seq, out);

                while let Some(v) = buffer.remove(&next_seq) {
                    if records_tx.send(v).is_err() {
                        return;
                    }
                    next_seq = next_seq.wrapping_add(1);
                }
            }
            drop(records_tx);
        });

        let mut errs = FxHashMap::default();
        let mut plot_points = Vec::with_capacity(n_genes * self.top_k);
        let gene_labels = self.get_gene_labels(&multi_progress)?;
        for result in records_rx {
            match result {
                Ok(mut gene_tx_dmr) => {
                    let records_written = gene_tx_dmr.write(
                        &mut writer,
                        single_mod_code,
                        self.emit_full_results,
                    )?;
                    if self.plot.is_some() {
                        let points = gene_tx_dmr.topk_records(
                            self.top_k,
                            self.min_effect_size,
                            &gene_labels,
                        );
                        plot_points.extend(points);
                    }
                    gene_tx_dmr.clear();
                    let _ = empties_tx.send(gene_tx_dmr);
                    records_counter.inc(records_written as u64);
                    pb.inc(1);
                }
                Err(e) => {
                    let c = errs.entry(e.to_string()).or_insert(0usize);
                    *c = c.saturating_add(1);
                }
            }
        }

        if let Some(fp) = self.plot.as_ref() {
            multi_progress.suspend(|| {
                info!("plotting {} points to {fp:?}", plot_points.len())
            });
            let svg = volcano_svg(&plot_points, self.plot_title.as_ref());
            std::fs::write(fp, svg)?;
        }

        source_thread.join().expect("source thread paniced");
        for (i, worker_thread) in handles.into_iter().enumerate() {
            worker_thread.join().expect(&format!("worker {i} paniced"));
        }
        aggregator.join().expect("aggregator theread paniced");

        multi_progress.suspend(|| {
            info!("finished, {} errors", errs.len());
        });

        if !errs.is_empty() {
            let err_table = format_errors_table(&errs);
            multi_progress.suspend(|| info!("{err_table}"));
        }

        Ok(())
    }

    fn get_writer(
        &self,
        single_mod_code: Option<ModCodeRepr>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Box<dyn Write>> {
        let mut writer: Box<dyn Write> = match self.output_table.as_ref() {
            "-" | "stdout" => {
                multi_progress.suspend(|| info!("writing to stdout"));
                Box::new(BufWriter::new(std::io::stdout()))
            }
            p @ _ => {
                multi_progress.suspend(|| info!("writing to file at {p}"));
                let fh = File::create(Path::new(p))
                    .context("failed to open output file")?;
                Box::new(BufWriter::new(fh))
            }
        };
        writer.write(
            GeneIsoformDmrRecord::<GeneDmrScore>::header(
                single_mod_code,
                self.emit_full_results,
            )
            .as_bytes(),
        )?;
        Ok(writer)
    }

    fn get_gene_labels(
        &self,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<HashSet<String>> {
        if let Some(fp) = self.gene_labels.as_ref() {
            let reader = BufReader::new(File::open(fp)?);
            let names = reader
                .lines()
                .filter_ok(|x| x.as_str() != "")
                .filter_ok(|x| !x.starts_with("#"))
                .collect::<Result<HashSet<String>, _>>()?;
            if names.is_empty() {
                bail!("failed to parse any gene names from {fp:?}")
            } else {
                multi_progress
                    .suspend(|| info!("parsed {} genes to label", names.len()));
                debug!("labeling genes: {}", names.iter().sorted().join(","));
            };

            Ok(names)
        } else {
            Ok(HashSet::new())
        }
    }
}
