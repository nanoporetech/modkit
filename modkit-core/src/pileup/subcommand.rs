use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::Args;
use common_macros::hash_set;
use crossbeam_channel::bounded;
use indicatif::{MultiProgress, ParallelProgressIterator};
use itertools::Itertools;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, HeaderView, Read};

use modkit_logging::init_logging;

use crate::command_utils::{
    get_motif_lookup_from_parts, get_threshold_from_options,
    parse_edge_filter_input, parse_per_base_thresholds,
    parse_per_mod_thresholds, parse_raw_motifs,
    parse_raw_thresholds_string_with_default, parse_thresholds,
    parse_thresholds_values,
};
use crate::fasta::MotifLocationsLookup;
use crate::interval_chunks::{
    ChromCoordinatesFeeder, ReferenceIntervalBatchesFeeder, TotalLength,
};
use crate::mod_bam::{CollapseMethod, EdgeFilter, MmTagInfo};
use crate::mod_base_code::{
    DnaBase, ModCodeRepr, ModifiedBasesOptions, ANY_ADENINE, ANY_CYTOSINE,
    ANY_GUANINE, ANY_THYMINE, METHYL_CYTOSINE, SIX_METHYL_ADENINE,
};
use crate::motifs::motif_bed::{MotifInfo, RegexMotif};
use crate::pileup::bedrmod::BedRModArgs;
use crate::pileup::duplex::{process_region_duplex_batch, DuplexModBasePileup};
use crate::pileup::pileup_processor::{
    CountsMatrix, DnaAllContext, DnaCpGCombineStrands, DnaCytosineCombine,
    DnaModOption, DnaPileupWorker, Dynamic, GenericPileupWorker, PileupWorker,
    DNA_BASES_CYTOSINE_FIRST,
};
use crate::pileup::{ModBasePileup2, PileupNumericOptions};
use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::sampling_schedule::IdxStats;
use crate::sample_probs::{
    calculate_reads_per_contig, run_extract_probs_workers,
    AlignedBaseArgmaxProbs, BaseArgmaxProbs, ExtractProbsWorker,
    ProbsExtractor, RegionMleProbs,
};
use crate::util::{
    create_out_directory, filter_reference_records, get_master_progress_bar,
    get_master_progress_bar_fancy, get_subroutine_progress_bar, get_targets,
    get_ticker, reader_is_bam, reader_is_cram, record_is_not_primary,
    unindexed_reader_is_cram, Region, Strand,
};
use crate::writers::{
    BedMethylWriter, BedMethylWriter2, MultipleMotifBedmethylWriter,
    PhasedBedMethylWriter, PileupWriter,
};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct ModBamPileup {
    // running args
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,
    /// Output file to write results into. Specify "-" or "stdout" to direct
    /// output to stdout.
    out_bed: String,
    /// Output bgzf-compressed bedmethyl files suitable for Tabix-indexing
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false, hide_short_help = false)]
    bgzf: bool,
    /// Specify the number of threads to dedicate to BGZF compression.
    #[clap(help_heading = "Compute Options")]
    #[arg(
        long,
        default_value_t = 4usize,
        hide_short_help = true,
        requires = "bgzf"
    )]
    bgzf_threads: usize,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Process only the specified region of the BAM when performing pileup.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>. Commas are
    /// allowed.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,

    // processing args
    /// Number of threads to use while processing chunks concurrently.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Number of IO/decompression threads to use
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: u32,
    /// Use this many threads when performing threshold estimation, setting
    /// this to a higher value can speed up execution.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, hide_short_help = true)]
    sampling_threads: Option<usize>,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(
        short = 'i',
        long,
        default_value_t = 1_000_000,
        hide_short_help = true
    )]
    interval_size: u32,
    /// Size of queue for writing records, default will be the number of
    /// threads.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, hide_short_help = true)]
    queue_size: Option<usize>,

    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Maximum number of reads to use whilest tabulating counts at a given
    /// position.
    #[clap(help_heading = "Compute Options")]
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = 8_000u16,
        requires = "high_depth"
    )]
    max_depth: u16,
    /// Flag to indicate that modBAM has exceptionally high depth (>65,000X)
    /// and should be downsampled during pileup.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, hide_short_help = false, default_value_t = false)]
    high_depth: bool,

    // sampling args
    /// Sample this many reads when estimating the filtering threshold. Reads
    /// will be sampled evenly across aligned genome. If a region is
    /// specified, either with the --region option or the --sample-region
    /// option, then reads will be sampled evenly across the region given.
    /// This option is useful for large BAM files. In practice, 10-50
    /// thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Fail if we cannot sample at least the requested number of reads when
    /// estimating the pass threshold
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        long,
        alias = "strict",
        conflicts_with_all = ["no_filtering", "sampling_frac"],
        default_value_t = false,
        hide_short_help = true,
    )]
    strict_sampling: bool,
    /// Sample this fraction of the reads when estimating the pass-threshold.
    /// In practice, 10-100 thousand reads is sufficient to estimate the model
    /// output distribution and determine the filtering threshold. See
    /// filtering.md for details on filtering.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Set a random seed for deterministic running, the default is
    /// non-deterministic.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Do not perform any filtering, include all mod base calls in output. See
    /// filtering.md for details on filtering.
    #[clap(help_heading = "Filtering Options")]
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per-base. Global filter
    /// threshold can be specified with by a decimal number (e.g. 0.75).
    /// Per-base thresholds can be specified by colon-separated values, for
    /// example C:0.75 specifies a threshold value of 0.75 for cytosine
    /// modification calls. Additional per-base thresholds can be specified
    /// by repeating the option: for example --filter-threshold C:0.75
    /// --filter-threshold A:0.70 or specify a single base option and a
    /// default for all other bases with: --filter-threshold A:0.70
    /// --filter-threshold 0.9 will specify a threshold value of 0.70 for
    /// adenine and 0.9 for all other base modification calls.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long,
    group = "thresholds",
    action = clap::ArgAction::Append,
    alias = "pass_threshold"
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long = "mod-threshold",
    alias = "mod-thresholds",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Maximum threhsold value to use, if the estimated threshold value for
    /// any primary sequence base is greater than this value, use this value
    /// instead. To set a default value for all primary bases, use
    /// --max-threshold 0.8, for example. Or use per-base threshold values like
    /// --max-threshold C:0.8.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        long,
        group = "thresholds",
        alias = "max-threshold",
        action = clap::ArgAction::Append,
    )]
    max_filter_threshold: Option<Vec<String>>,
    /// Specify a region for sampling reads from when estimating the threshold
    /// probability. If this option is not provided, but --region is
    /// provided, the genomic interval passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[clap(help_heading = "Filtering Options")]
    #[arg(long)]
    sample_region: Option<String>,
    /// Interval chunk size in base pairs to process concurrently when
    /// estimating the threshold probability, can be larger than the pileup
    /// processing interval.
    #[clap(help_heading = "Filtering Options")]
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,
    /// BED file that will restrict threshold estimation and pileup results to
    /// positions overlapping intervals in the file. (alias: include-positions)
    #[clap(help_heading = "Selection Options")]
    #[arg(long, hide_short_help = true, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Include unmapped base modifications when estimating the pass threshold.
    #[clap(help_heading = "Selection Options")]
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = false,
        conflicts_with = "include_bed"
    )]
    include_unmapped: bool,

    #[arg(
        long,
        hide_long_help = true,
        hide_short_help = true,
        default_value_t = false
    )]
    use_dynamic: bool,

    /// Specify that modBAM contains duplex modification calls and these should
    /// be used.
    #[arg(long, hide_short_help = true, default_value_t = false)]
    duplex: bool,

    /// Output pileup counts for only sequence motifs provided. The first
    /// argument should be the sequence motif and the second argument is
    /// the 0-based offset to the base to pileup base modification counts
    /// for. For example: --motif CGCG 0 indicates to pileup counts for the
    /// first C on the top strand and the last C (complement to G) on
    /// the bottom strand. The --cpg argument is short hand for --motif CG 0.
    ///
    /// This argument can be passed multiple times. When more than one motif is
    /// used, the resulting output BED file will indicate the motif in the
    /// "name" field as <mod_code>,<motif>,<offset>. For example, given
    /// `--motif CGCG 2 --motif CG 0` there will be output lines with name
    /// fields such as "m,CG,0" and "m,CGCG,2". To use this option with
    /// `--combine-strands`, all motifs must be reverse-complement
    /// palindromic or an error will be raised.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, action = clap::ArgAction::Append, num_args = 2, requires = "reference_fasta")]
    motif: Option<Vec<String>>,
    /// Only output counts at CpG motifs.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, requires_all = ["reference_fasta", "modified_bases"], default_value_t = false)]
    cpg: bool,
    /// Reference sequence in FASTA format. Required for motif (e.g. CpG)
    /// filtering, requires FAI fasta index to be pre-generated. (alias: 'ref')
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long = "reference", alias = "ref", short = 'r')]
    reference_fasta: Option<PathBuf>,
    /// Preload the reference sequences, useful when working with many, short
    /// reference sequences such as a transcriptome.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, default_value_t = false)]
    preload_references: bool,
    #[clap(help_heading = "Modified Base Options")]
    /// Respect soft masking in the reference FASTA.
    #[arg(
        long,
        short = 'k',
        requires = "reference_fasta",
        default_value_t = false,
        hide_short_help = true
    )]
    mask: bool,
    /// Specify which modified bases to tablulate counts for. These can be
    /// the "long name" such as '5mC', '6mA' (or 'm6A' for RNA), or 'inosine'.
    /// You can also pass <primary_base>:<mod_code>, such as 'C:m'. Finally,
    /// when running with --combine-mods the reference base can be passed
    /// alone, such as 'C' or 'A'.
    #[arg(
        long,
        requires = "reference_fasta",
        conflicts_with_all = ["duplex"],
        value_parser = clap::value_parser!(ModifiedBasesOptions),
        num_args = 1..,
        value_delimiter = ' '
    )]
    modified_bases: Option<Vec<ModifiedBasesOptions>>,
    /// Combine base modification calls, all counts of modified bases are
    /// summed together. See collapse.md for details.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(
        long,
        default_value_t = false,
        group = "combine_args",
        hide_short_help = true
    )]
    combine_mods: bool,
    /// Use non-primary alignments in pileup (e.g. secondary and
    /// supplementaty). Keep in mind these mappings will be used _in
    /// addition_ to the primary mappings.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(
        long,
        default_value_t = false,
        hide_short_help = true,
        requires = "modified_bases"
    )]
    allow_non_primary: bool,
    /// When performing motif analysis (such as CpG), sum the counts from the
    /// positive and negative strands into the counts for the positive
    /// strand position. Requires that the BAM records have correct MN tags.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, default_value_t = false)]
    combine_strands: bool,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, hide_short_help = true, requires = "modified_bases")]
    edge_filter: Option<String>,

    /// Output bedMethyl where the delimiter of columns past column 10 are
    /// space-delimited instead of tab-delimited. This option can be useful
    /// for some browsers and parsers that don't expect the extra columns
    /// of the bedMethyl format.
    // #[clap(help_heading = "Output Options")]
    // #[arg(
    //     long = "mixed-delim",
    //     alias = "mixed-delimiters",
    //     default_value_t = false,
    //     hide_short_help = true
    // )]
    // mixed_delimiters: bool,
    /// Output a header with the bedMethyl
    #[clap(help_heading = "Output Options")]
    #[arg(
        long = "header",
        alias = "with-header",
        alias = "include_header",
        default_value_t = false
    )]
    with_header: bool,
    /// Prefix to prepend on phased output file names. Without this option
    /// the files will be <mod_code>_<strand>.bedmethyl
    #[clap(help_heading = "Output Options")]
    #[arg(long)]
    prefix: Option<String>,
    /// Produce separate bedMethyl tables for haplotype-labeled modBAM records.
    /// Currently only supports diploid genomes. Produces hp1.bedmethyl,
    /// hp2.bedmethyl, and combined.bedmethyl files (optionally
    /// compressed). hp1.bedmethyl and hp2.bedmethyl contain counts for
    /// records with HP=1 and HP=2 tags, respectively. combined.bedmethyl
    /// contains counts for all modBAM records.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    phased: bool,
    #[clap(flatten)]
    bedrmodargs: BedRModArgs,
}

impl ModBamPileup {
    fn parse_user_motifs(&self) -> Option<anyhow::Result<Vec<RegexMotif>>> {
        if let Some(raw_motif_parts) = &self.motif {
            Some(parse_raw_motifs(raw_motif_parts, self.cpg, false))
        } else if self.cpg {
            Some(Ok(vec![RegexMotif::parse_string("CG", 0).unwrap()]))
        } else {
            None
        }
    }

    fn use_cpg(motifs: &[RegexMotif]) -> bool {
        motifs == &[RegexMotif::parse_string("CG", 0).unwrap()]
    }

    fn is_5mc_5hmc_cpg(
        modified_bases: &[(DnaBase, ModCodeRepr)],
    ) -> anyhow::Result<bool> {
        if modified_bases.iter().any(|x| x.0 != DnaBase::C) {
            let non_cytosine = modified_bases
                .iter()
                .filter_map(|(b, _)| {
                    if *b != DnaBase::C {
                        Some(b.char())
                    } else {
                        None
                    }
                })
                .join(",");
            bail!(
                "cannot use CpG with primary bases other than Cytosine, got \
                 {non_cytosine}"
            )
        }
        let n_mod_codes = modified_bases.iter().map(|x| x.1).unique().count();
        Ok(n_mod_codes <= 2)
    }

    fn only_cytosine(modified_bases: &[(DnaBase, ModCodeRepr)]) -> bool {
        modified_bases.iter().map(|(b, _)| *b).unique().collect_vec()
            == vec![DnaBase::C]
    }

    // this Baroque method checks that we have A and/or C mods only, <= 2 C mods
    // and <= 1 A mod
    fn is_simple_dna(modified_bases: &[(DnaBase, ModCodeRepr)]) -> bool {
        assert!(!modified_bases.is_empty());
        if modified_bases.len() > 3 {
            return false;
        }

        let a_and_c = hash_set! { DnaBase::A, DnaBase::C };
        let primary_bases =
            modified_bases.iter().map(|x| x.0).collect::<HashSet<DnaBase>>();
        let is_subset_of_a_and_c = primary_bases.is_subset(&a_and_c);

        if !is_subset_of_a_and_c {
            return false;
        }
        if !primary_bases.contains(&DnaBase::C) {
            return false;
        }

        let simple_mod_choices = modified_bases
            .iter()
            .sorted_by_key(|(b, _)| *b)
            .group_by(|(b, _)| *b)
            .into_iter()
            .all(|(base, it)| match base {
                DnaBase::A => {
                    let a_mods =
                        it.map(|(_a, y)| *y).collect::<Vec<ModCodeRepr>>();
                    a_mods.len() == 1 && a_mods.contains(&SIX_METHYL_ADENINE)
                }
                DnaBase::C => {
                    let mut has_methyl_c = false;
                    let mut count = 0usize;
                    for (prim_base, code) in it {
                        debug_assert_eq!(*prim_base, DnaBase::C);
                        count = count.saturating_add(1);
                        if *code == METHYL_CYTOSINE {
                            has_methyl_c = true;
                        }
                    }
                    count <= 2usize && has_methyl_c
                }
                DnaBase::G => unreachable!(),
                DnaBase::T => unreachable!(),
            });

        simple_mod_choices
    }

    fn check_any_mod_codes(
        modified_bases: &[(DnaBase, ModCodeRepr)],
        combine_mods: bool,
    ) -> anyhow::Result<()> {
        for (base, codes) in &modified_bases.iter().group_by(|(x, _)| *x) {
            let has_anymod_code =
                codes.collect_vec().iter().any(|(_, x)| match base {
                    DnaBase::A => *x == ANY_ADENINE,
                    DnaBase::C => *x == ANY_CYTOSINE,
                    DnaBase::G => *x == ANY_GUANINE,
                    DnaBase::T => *x == ANY_THYMINE,
                });
            if has_anymod_code && !combine_mods {
                bail!(
                    "must explicitly pass --combine-mods when passing primary \
                     bases alone"
                )
            }
        }
        Ok(())
    }

    fn determine_preset(
        &self,
    ) -> anyhow::Result<(Option<Presets>, Option<Vec<RegexMotif>>)> {
        let mut regex_motifs =
            self.parse_user_motifs().transpose()?.unwrap_or_else(Vec::new);
        if regex_motifs.len() > 1 {
            if self.combine_strands {
                bail!(
                    "multiple motifs and combine-strands not currently \
                     supported"
                )
            }
            info!(
                "more than one motif requires use of general pileup processor"
            );
            let motif_primary_bases = regex_motifs
                .iter()
                .map(|mot| mot.motif_info.primary_base)
                .collect::<HashSet<DnaBase>>();
            if let Some(modified_bases) = self.modified_bases.as_ref() {
                for primary_base in
                    modified_bases.iter().map(|x| x.primary_base)
                {
                    if !motif_primary_bases.contains(&primary_base) {
                        info!(
                            "adding single-base motif: '{} 0'",
                            primary_base.char()
                        );
                        regex_motifs.push(
                            RegexMotif::parse_string(
                                &primary_base.char().to_string(),
                                0,
                            )
                            .unwrap(),
                        );
                    }
                }
            }
            return Ok((None, Some(regex_motifs)));
        }
        if let Some(modified_bases) = self.modified_bases.as_ref() {
            let modified_bases = modified_bases
                .iter()
                .map(|x| (x.primary_base, x.mod_code))
                .sorted()
                .collect::<Vec<_>>();
            if modified_bases.is_empty() {
                bail!("must provide at least one modified base option")
            }
            Self::check_any_mod_codes(&modified_bases, self.combine_mods)?;
            if self.combine_mods {
                let dna_bases_list = modified_bases
                    .iter()
                    .map(|(x, _)| x.char())
                    .unique()
                    .join(",");
                info!(
                    "parsed {} base modification(s), combining modification \
                     counts together for primary sequence bases \
                     '{dna_bases_list}'",
                    modified_bases.len()
                );
            } else {
                let modified_bases_list = modified_bases
                    .iter()
                    .map(|(b, x)| format!("{b}:{x}"))
                    .join(",");
                info!(
                    "parsed {} base modification(s). Base modifications other \
                     than '{modified_bases_list}' will be counted as \
                     'N_other'.",
                    modified_bases.len()
                );
            }

            if self.combine_strands
                && Self::use_cpg(&regex_motifs)
                && Self::is_5mc_5hmc_cpg(&modified_bases)?
                && !self.high_depth
            {
                debug!("using optimized CpG combine strands processor");
                let preset = Presets::DnaCpGCombineStrands {
                    cytosine_other_mod: modified_bases.iter().find_map(
                        |(_, x)| {
                            if *x != METHYL_CYTOSINE {
                                Some(*x)
                            } else {
                                None
                            }
                        },
                    ),
                };
                return Ok((Some(preset), Some(regex_motifs)));
            }

            if self.combine_strands {
                if regex_motifs.is_empty() {
                    bail!(
                        "need to provide  one reverse-complement palindromic \
                         motif to combine strands"
                    )
                }
                let motif_offset = regex_motifs[0].motif_info.offset();
                if motif_offset < 0 {
                    bail!("invalid palindromic motif");
                }
                let motif_offset = motif_offset as u32;
                let motif_bases = [
                    regex_motifs[0].motif_info.primary_base,
                    DnaBase::C,
                    DnaBase::C,
                    DnaBase::C,
                ];
                let preset = Presets::DynamicAllContext {
                    mod_codes: modified_bases,
                    motif_bases,
                    motif_offset: Some(motif_offset),
                };
                debug!(
                    "using optimized non-CpG single-motif combine strands \
                     processor"
                );
                return Ok((Some(preset), Some(regex_motifs)));
            }

            if self.combine_mods
                && Self::only_cytosine(&modified_bases)
                && !self.high_depth
            {
                debug!("using optimized Cytosine combine mods processor");
                let preset = Presets::DnaCytosineCombine;
                let c_motif = RegexMotif::parse_string("C", 0).unwrap();
                if !regex_motifs.contains(&c_motif) {
                    regex_motifs.push(c_motif);
                }
                return Ok((Some(preset), Some(regex_motifs)));
            }

            let mut motif_primary_bases = regex_motifs
                .iter()
                .map(|x| x.motif_info.primary_base)
                .collect::<BTreeSet<DnaBase>>();
            for (dna_base, _) in modified_bases.iter() {
                if !motif_primary_bases.contains(dna_base) {
                    info!("adding single-base motif: '{} 0'", dna_base.char());
                    regex_motifs.push(
                        RegexMotif::parse_string(
                            &dna_base.char().to_string(),
                            0,
                        )
                        .unwrap(),
                    );
                    motif_primary_bases.insert(*dna_base);
                }
            }
            let mut motif_bases =
                [motif_primary_bases.iter().nth(0).copied().unwrap(); 4];
            for (i, b) in motif_primary_bases.into_iter().enumerate() {
                motif_bases[i] = b;
            }

            if Self::is_simple_dna(&modified_bases)
                && !(self.use_dynamic || self.high_depth)
            {
                let other_cytosine_mod =
                    modified_bases.iter().find_map(|(b, c)| {
                        if *b == DnaBase::C && *c != METHYL_CYTOSINE {
                            Some(*c)
                        } else {
                            None
                        }
                    });

                let preset = Presets::DnaAllContext {
                    cytosine_mod_op: other_cytosine_mod,
                    motif_bases,
                };
                return Ok((Some(preset), Some(regex_motifs)));
            } else {
                let preset = Presets::DynamicAllContext {
                    mod_codes: modified_bases,
                    motif_bases,
                    motif_offset: None,
                };
                return Ok((Some(preset), Some(regex_motifs)));
            }
        } else {
            return Ok((
                None,
                if regex_motifs.is_empty() { None } else { Some(regex_motifs) },
            ));
        }
    }

    fn calc_per_base_thresholds(
        &self,
        regex_motifs: Option<Vec<RegexMotif>>,
        stranded_position_filter: Option<StrandedPositionFilter<()>>,
        bam_header: &HeaderView,
        thread_pool: &rust_htslib::tpool::ThreadPool,
        edge_filter: Option<&EdgeFilter>,
        preset: Option<&Presets>,
        sampling_frac: Option<f64>,
        sampling_region: Option<&Region>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<[f32; 4]> {
        let max_per_base_trehsolds = self
            .max_filter_threshold
            .as_ref()
            .map(|x| parse_raw_thresholds_string_with_default(x))
            .transpose()?;

        let motif_lookup = get_motif_lookup_from_parts(
            regex_motifs.clone(),
            self.reference_fasta.as_ref(),
            false,
            self.mask,
            self.preload_references,
        )?;
        let reference_records = get_targets(bam_header, sampling_region);
        let feeder = ChromCoordinatesFeeder::new(
            &reference_records,
            self.sampling_interval_size,
            motif_lookup,
            false,
            stranded_position_filter.clone(),
        )?;
        let mut workers: Vec<Box<dyn ExtractProbsWorker>> = Vec::new();
        let n_workers = self.sampling_threads.unwrap_or(self.threads);
        let (chrom_to_counts, rng_sample, sample_frac) =
            calculate_reads_per_contig(
                bam_header,
                sampling_frac,
                Some(self.num_reads),
                &multi_progress,
                self.sampling_interval_size,
                sampling_region,
            )?;
        let chrom_to_counts = chrom_to_counts.map(|x| Arc::new(x));

        if let Some(preset) = preset {
            let motif_bases = match preset {
                Presets::DnaCpGCombineStrands { cytosine_other_mod: _ }
                | Presets::DnaCytosineCombine => {
                    multi_progress.suspend(|| {
                        info!(
                            "calculating threshold value with probabilities \
                             aligned to reference bases cytosines",
                        )
                    });
                    [DnaBase::C; 4]
                }
                Presets::DnaAllContext { cytosine_mod_op: _, motif_bases }
                | Presets::DynamicAllContext {
                    mod_codes: _,
                    motif_bases,
                    motif_offset: _,
                } => {
                    multi_progress.suspend(|| {
                        info!(
                            "calculating threshold value with probabilities \
                             aligned to reference base(s) {}",
                            motif_bases.iter().unique().join(",")
                        )
                    });
                    *motif_bases
                }
            };
            for i in 0..n_workers {
                let worker = RegionMleProbs::<
                    AlignedBaseArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    &self.in_bam,
                    self.reference_fasta.as_ref(),
                    false,
                    motif_bases,
                    edge_filter,
                    thread_pool,
                    self.seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?;
                workers.push(Box::new(worker));
            }
        } else if let Some(motif_bases) = feeder.get_motif_bases() {
            multi_progress.suspend(|| {
                info!(
                    "calculating threshold value with probabilities aligned \
                     to reference base(s) {}",
                    motif_bases.iter().unique().join(",")
                )
            });
            for i in 0..n_workers {
                let worker = RegionMleProbs::<
                    AlignedBaseArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    &self.in_bam,
                    self.reference_fasta.as_ref(),
                    false,
                    motif_bases,
                    edge_filter,
                    thread_pool,
                    self.seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?;
                workers.push(Box::new(worker));
            }
        } else {
            multi_progress.suspend(|| {
                info!(
                    "calculating threshold value with probabilites from reads",
                )
            });
            for i in 0..n_workers {
                let worker =
                    RegionMleProbs::<BaseArgmaxProbs, ProbsExtractor>::new(
                        &self.in_bam,
                        self.reference_fasta.as_ref(),
                        false,
                        [DnaBase::A; 4],
                        edge_filter,
                        thread_pool,
                        self.seed.unwrap_or(i as u64),
                        rng_sample,
                        sample_frac,
                        chrom_to_counts.clone(),
                    )?;
                workers.push(Box::new(worker));
            }
        }

        let qual_hist = run_extract_probs_workers(
            workers,
            multi_progress.clone(),
            feeder.clone(),
        )?;
        if qual_hist.ok_records < self.num_reads
            && sampling_frac.is_none()
            && self.strict_sampling
        {
            let sample_frac = 0.2f64;
            multi_progress.suspend(|| {
                warn!(
                    "failed to sample enough records, got {}, looking for {}, \
                     trying again sampling {:.0}% of reads",
                    qual_hist.ok_records,
                    self.num_reads,
                    sample_frac * 100f64
                );
            });
            return self.calc_per_base_thresholds(
                regex_motifs.clone(),
                stranded_position_filter.clone(),
                bam_header,
                thread_pool,
                edge_filter,
                preset,
                Some(sample_frac),
                sampling_region,
                multi_progress,
            );
        }
        if qual_hist.ok_records == 0 {
            bail!(
                "Failed to sample _any_ records to estimate threshold. Use \
                 --sample-frac with a value greater than {sample_frac:.2}, \
                 lower --num-reads, or use --no-filtering"
            )
        }
        qual_hist.get_base_thresholds(
            self.filter_percentile,
            max_per_base_trehsolds,
            &multi_progress,
        )
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let master_progress = MultiProgress::new();
        master_progress
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());

        if !self.suppress_progress {
            master_progress
                .set_draw_target(indicatif::ProgressDrawTarget::stderr());
        }

        let io_thread_pool =
            rust_htslib::tpool::ThreadPool::new(self.io_threads)?;

        if self.filter_percentile > 1.0 {
            bail!("filter percentile must be <= 1.0")
        }

        // do this first so we fail when the file isn't readable
        let (header, check_idx) = bam::IndexedReader::from_path(&self.in_bam)
            .map_err(|e| e.into())
            .and_then(|reader| {
                let check = reader_is_bam(&reader);
                if self.reference_fasta.is_none() && reader_is_cram(&reader) {
                    bail!("--reference is required when using CRAM")
                } else {
                    Ok((reader.header().to_owned(), check))
                }
            })?;

        // options parsing below
        let region = self
            .region
            .as_ref()
            .map(|raw_region| {
                info!("parsing region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let sampling_region = self
            .sample_region
            .as_ref()
            .map(|raw_region| {
                info!("parsing sample region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|trims| parse_edge_filter_input(trims, false))
            .transpose()?;
        let per_mod_thresholds = self
            .mod_thresholds
            .as_ref()
            .map(|raw_per_mod_thresholds| {
                parse_per_mod_thresholds(raw_per_mod_thresholds)
            })
            .transpose()?;
        let reference_records = get_targets(&header, region.as_ref());
        if self.include_bed.is_some() {
            if !(self.modified_bases.is_some() || self.motif.is_some()) {
                bail!(
                    "currently, --include-bed requires a --motif or \
                     --modified-bases. This limitation may be removed in the \
                     future."
                )
            }
        }
        let position_filter = self
            .include_bed
            .as_ref()
            .map(|bed_fp| {
                let chrom_to_tid = reference_records
                    .iter()
                    .map(|reference_record| {
                        (reference_record.name.as_str(), reference_record.tid)
                    })
                    .collect::<HashMap<&str, u32>>();
                StrandedPositionFilter::from_bed_file(
                    bed_fp,
                    &chrom_to_tid,
                    self.suppress_progress,
                )
            })
            .transpose()?;

        let reference_records = if check_idx {
            filter_reference_records(
                reference_records,
                &self.in_bam,
                region.as_ref(),
                position_filter.as_ref(),
                &master_progress,
            )?
        } else {
            info!(
                "not checking for mapped reads or filtering contigs with CRAM \
                 input"
            );
            reference_records
        };

        let (preset, regex_motifs) = self.determine_preset()?;
        // The optimized workers do not support modification calls on the
        // opposite strand of the primary sequence base (e.g. PacBio 6mA
        // `T-a.` sub-tags), check the first records and fall back to the
        // general workers when such calls are found.
        let preset = if preset.is_some()
            && bam_has_negative_strand_mods(
                &self.in_bam,
                self.reference_fasta.as_ref(),
                1_000,
            )? {
            info!(
                "found base modification calls on the opposite strand of the \
                 primary sequence base (e.g. PacBio 6mA 'T-a.'), these are \
                 not supported by the optimized pileup workers, using general \
                 workers"
            );
            None
        } else {
            preset
        };

        let (pileup_options, combine_strands) = match &preset {
            Some(preset) => match preset {
                Presets::DnaCpGCombineStrands { .. } => {
                    (PileupNumericOptions::Passthrough, true)
                }
                _ => (PileupNumericOptions::Passthrough, false),
            },
            None => {
                let options = if self.combine_mods {
                    PileupNumericOptions::Combine
                } else {
                    PileupNumericOptions::Passthrough
                };
                (options, self.combine_strands)
            }
        };

        let (empties_tx, empties_rx) = crossbeam_channel::unbounded();
        let out_fp_str = self.out_bed.clone();
        let motif_labels = regex_motifs
            .as_ref()
            .map(|regex_motifs| {
                regex_motifs
                    .iter()
                    .map(|mot| format!("{}", mot))
                    .collect::<Vec<String>>()
            })
            .unwrap_or(Vec::new());

        // setup the writer here so we fail before doing any work (if there are
        // problems).
        let mut writer: Box<dyn PileupWriter<ModBasePileup2>> = if self.phased {
            let out_dir = Path::new(&self.out_bed);
            if out_dir.exists() && !out_dir.is_dir() {
                bail!("phased output needs to point to a directory");
            }
            if !out_dir.exists() {
                info!("creating directory at {out_dir:?}");
                std::fs::create_dir_all(out_dir)?;
            }
            info!("producing phased output");
            if self.bgzf {
                info!(
                    "using bgzf compression with {} compression threads",
                    self.bgzf_threads
                );
                Box::new(PhasedBedMethylWriter::new_bgzf(
                    &out_dir.to_path_buf(),
                    self.prefix.as_ref(),
                    true, // TODO
                    master_progress.clone(),
                    empties_tx.clone(),
                    self.bgzf_threads,
                )?)
            } else {
                Box::new(PhasedBedMethylWriter::new_file(
                    &out_dir.to_path_buf(),
                    self.prefix.as_ref(),
                    true, // TODO
                    master_progress.clone(),
                    empties_tx.clone(),
                )?)
            }
        } else {
            let have_multiple_motifs =
                regex_motifs.as_ref().map(|x| x.len() > 1).unwrap_or(false);
            let use_special_writer = have_multiple_motifs && preset.is_none();
            match out_fp_str.as_str() {
                "stdout" | "-" => {
                    if self.bgzf {
                        bail!("bgzf compression requires file output")
                    }
                    if have_multiple_motifs && use_special_writer {
                        debug!("using multiple-motif stdout writer");
                        Box::new(MultipleMotifBedmethylWriter::new_stdout(
                            self.with_header,
                            &self.bedrmodargs,
                            &header,
                            self.modified_bases.as_ref(),
                            empties_tx.clone(),
                            master_progress.clone(),
                        )?)
                    } else {
                        debug!("using standard stdout writer");
                        Box::new(BedMethylWriter2::new_stdout(
                            self.with_header,
                            &self.bedrmodargs,
                            &header,
                            self.modified_bases.as_ref(),
                            master_progress.clone(),
                            empties_tx.clone(),
                        )?)
                    }
                }
                _ => {
                    create_out_directory(&out_fp_str)?;
                    if use_special_writer {
                        debug!("using multiple-motif writer");
                        if self.bgzf {
                            info!(
                                "using bgzf compression with {} compression \
                                 threads",
                                self.bgzf_threads
                            );
                            Box::new(MultipleMotifBedmethylWriter::new_bgzf(
                                &Path::new(&out_fp_str).to_path_buf(),
                                self.with_header,
                                &self.bedrmodargs,
                                &header,
                                self.modified_bases.as_ref(),
                                master_progress.clone(),
                                empties_tx.clone(),
                                self.bgzf_threads,
                            )?)
                        } else {
                            Box::new(MultipleMotifBedmethylWriter::new_file(
                                &Path::new(&out_fp_str).to_path_buf(),
                                self.with_header,
                                &self.bedrmodargs,
                                &header,
                                self.modified_bases.as_ref(),
                                master_progress.clone(),
                                empties_tx.clone(),
                            )?)
                        }
                    } else {
                        if self.bgzf {
                            info!(
                                "using bgzf compression with {} compression \
                                 threads",
                                self.bgzf_threads
                            );
                            Box::new(BedMethylWriter2::new_bgzf(
                                &Path::new(&out_fp_str).to_path_buf(),
                                self.with_header,
                                &self.bedrmodargs,
                                &header,
                                self.modified_bases.as_ref(),
                                master_progress.clone(),
                                empties_tx.clone(),
                                self.bgzf_threads,
                            )?)
                        } else {
                            debug!("using standard file writer");
                            Box::new(BedMethylWriter2::new(
                                &Path::new(&out_fp_str).to_path_buf(),
                                self.with_header,
                                &self.bedrmodargs,
                                &header,
                                self.modified_bases.as_ref(),
                                master_progress.clone(),
                                empties_tx.clone(),
                            )?)
                        }
                    }
                }
            }
        };
        let motif_lookup = get_motif_lookup_from_parts(
            regex_motifs.clone(),
            self.reference_fasta.as_ref(),
            combine_strands,
            self.mask,
            self.preload_references,
        )?;

        // TODO: re-enable this after a think about the best way to "optimize"
        // let reference_records = if let Some(pf) = position_filter.as_ref() {
        //     info!(
        //         "optimizing reference records, started with {}",
        //         reference_records.len()
        //     );
        //     let rr = pf.optimize_reference_records(
        //         reference_records,
        //         self.interval_size,
        //     );
        //     info!("..finished with {}", rr.len());
        //     rr
        // } else {
        //     reference_records
        // };
        // start the actual work here
        let (default_threshold, base_thresholds) =
            if let Some(raw_threshold) = &self.filter_threshold {
                if preset.is_some() {
                    let (default_threshold, per_base_thresholds) =
                        parse_per_base_thresholds(&raw_threshold)?;
                    match default_threshold {
                        Some(t) if t <= 0.3f32 => {
                            warn!(
                                "--modified-bases uses simple thresholding \
                                 algorithm, using low threshold {t} may not \
                                 produce desired results"
                            )
                        }
                        _ => {}
                    }
                    for (b, t) in per_base_thresholds {
                        if t <= 0.3 {
                            warn!(
                                "--modified-bases uses simple thresholding \
                                 algorithm, using low threshold {t} for base \
                                 {b} may not produce desired results"
                            )
                        }
                    }
                }
                parse_thresholds_values(raw_threshold)?
            } else if self.no_filtering {
                info!("not performing filtering");
                (0f32, [0f32; 4])
            } else {
                (
                    0f32,
                    self.calc_per_base_thresholds(
                        regex_motifs.clone(),
                        position_filter.clone(),
                        &header,
                        &io_thread_pool,
                        edge_filter.as_ref(),
                        preset.as_ref(),
                        self.sampling_frac,
                        region.as_ref().or(sampling_region.as_ref()),
                        master_progress.clone(),
                    )?,
                )
            };
        let in_bam_fp = self.in_bam.clone();
        let per_mod_thresholds = per_mod_thresholds
            .map(|x| {
                x.into_iter()
                    .sorted_by_key(|(x, _)| {
                        if x.is_any() {
                            std::cmp::Ordering::Greater
                        } else {
                            std::cmp::Ordering::Less
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_else(Vec::new);

        debug!("per_mod_thresholds={per_mod_thresholds:?}");
        let mut workers: Vec<Box<dyn PileupWorker>> = Vec::new();
        match preset {
            Some(Presets::DnaCpGCombineStrands { cytosine_other_mod }) => {
                master_progress.suspend(|| {
                    info!(
                        "using optimized genome workers, combining strands at \
                         CpG motifs"
                    );
                });

                let ref_fasta = self.reference_fasta.as_ref().expect(
                    "modified bases option should require reference is \
                     provided",
                );

                let xs = (0..self.threads)
                    .map(|_| {
                        DnaPileupWorker::<
                            DnaCpGCombineStrands,
                            CountsMatrix,
                            1,
                            false,
                        >::new(
                            &in_bam_fp,
                            ref_fasta,
                            self.phased,
                            base_thresholds,
                            DNA_BASES_CYTOSINE_FIRST,
                            if self.combine_mods {
                                DnaModOption::Combine
                            } else {
                                if let Some(other) = cytosine_other_mod {
                                    DnaModOption::Other(other)
                                } else {
                                    DnaModOption::Single
                                }
                            },
                            self.allow_non_primary,
                            Vec::new(),
                            per_mod_thresholds.clone(),
                            edge_filter.as_ref(),
                            1,
                            u16::MAX,
                        )
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                for x in xs {
                    workers.push(Box::new(x));
                }
            }
            Some(Presets::DnaCytosineCombine) => {
                master_progress.suspend(|| {
                    info!(
                        "using optimized genome workers, summing Cytosine \
                         modification counts"
                    );
                });

                let ref_fasta = self.reference_fasta.as_ref().expect(
                    "modified bases option should require reference is \
                     provided",
                );
                let xs = (0..self.threads)
                    .map(|_| {
                        DnaPileupWorker::<
                            DnaCytosineCombine,
                            CountsMatrix,
                            2,
                            false,
                        >::new(
                            &in_bam_fp,
                            ref_fasta,
                            self.phased,
                            base_thresholds,
                            DNA_BASES_CYTOSINE_FIRST,
                            DnaModOption::Combine,
                            self.allow_non_primary,
                            Vec::new(),
                            per_mod_thresholds.clone(),
                            edge_filter.as_ref(),
                            0,
                            u16::MAX,
                        )
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                for x in xs {
                    workers.push(Box::new(x));
                }
            }
            Some(Presets::DnaAllContext { cytosine_mod_op, motif_bases }) => {
                master_progress.suspend(|| {
                    info!("using optimized workers for A,C all-context");
                });
                let ref_fasta = self.reference_fasta.as_ref().expect(
                    "modified bases option should require reference is \
                     provided",
                );
                let dna_mod_option = if self.combine_mods {
                    DnaModOption::Combine
                } else {
                    if let Some(other) = cytosine_mod_op {
                        DnaModOption::Other(other)
                    } else {
                        DnaModOption::Single
                    }
                };
                let xs = (0..self.threads)
                    .map(|_| {
                        DnaPileupWorker::<
                                DnaAllContext,
                                CountsMatrix,
                                2,
                                false,
                            >::new(
                                &in_bam_fp,
                                ref_fasta,
                                self.phased,
                                base_thresholds,
                                motif_bases,
                                dna_mod_option,
                                self.allow_non_primary,
                                Vec::new(),
                                per_mod_thresholds.clone(),
                                edge_filter.as_ref(),
                                0,
                                u16::MAX,
                            )
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                for x in xs {
                    workers.push(Box::new(x));
                }
            }
            Some(Presets::DynamicAllContext {
                mod_codes,
                motif_bases,
                motif_offset,
            }) => {
                master_progress.suspend(|| {
                    if motif_offset.is_some() {
                        info!(
                            "using optimized workers for all-context, \
                             combining modification counts across strands"
                        );
                    } else {
                        info!("using optimized workers for all-context");
                    }
                    if self.high_depth {
                        info!(
                            "using maximum of {} reads per position",
                            self.max_depth
                        );
                    }
                });
                let ref_fasta = self.reference_fasta.as_ref().expect(
                    "modified bases option should require reference is \
                     provided",
                );
                let mod_codes = if self.combine_mods {
                    mod_codes
                        .into_iter()
                        .map(|(base, _)| (base, ModCodeRepr::Code(base.char())))
                        .unique()
                        .sorted_by_key(|x| x.0)
                        .collect()
                } else {
                    mod_codes
                };
                let dna_mod_option = if self.combine_mods {
                    DnaModOption::Combine
                } else {
                    // doesn't matter
                    DnaModOption::Other(ANY_CYTOSINE)
                };
                for _ in 0..self.threads {
                    if let Some(motif_offset) = motif_offset {
                        let worker = DnaPileupWorker::<
                            Dynamic,
                            CountsMatrix,
                            1,
                            true,
                        >::new(
                            &in_bam_fp,
                            ref_fasta,
                            self.phased,
                            base_thresholds,
                            motif_bases,
                            dna_mod_option,
                            self.allow_non_primary,
                            mod_codes.clone(),
                            per_mod_thresholds.clone(),
                            edge_filter.as_ref(),
                            motif_offset as u32,
                            self.max_depth,
                        )?;
                        workers.push(Box::new(worker))
                    } else {
                        let worker = DnaPileupWorker::<
                            Dynamic,
                            CountsMatrix,
                            2,
                            true,
                        >::new(
                            &in_bam_fp,
                            ref_fasta,
                            self.phased,
                            base_thresholds,
                            motif_bases,
                            dna_mod_option,
                            self.allow_non_primary,
                            mod_codes.clone(),
                            per_mod_thresholds.clone(),
                            edge_filter.as_ref(),
                            0,
                            self.max_depth,
                        )?;
                        workers.push(Box::new(worker))
                    }
                }
            }
            _ => {
                master_progress.suspend(|| {
                    if self.duplex {
                        info!("using general workers for duplex");
                    } else {
                        info!("using general workers");
                    }
                });
                let motif_infos = motif_lookup
                    .as_ref()
                    .map(|x| {
                        x.motif_infos().copied().collect::<Vec<MotifInfo>>()
                    })
                    .unwrap_or_else(Vec::new);
                let xs = (0..self.threads)
                    .map(|_| {
                        GenericPileupWorker::new(
                            &self.in_bam,
                            self.reference_fasta.as_ref(),
                            motif_infos.clone(),
                            default_threshold,
                            base_thresholds,
                            &per_mod_thresholds,
                            pileup_options.clone(),
                            self.combine_strands,
                            self.max_depth,
                        )
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                for x in xs {
                    workers.push(Box::new(x));
                }
            }
        };

        let feeder = ChromCoordinatesFeeder::new(
            &reference_records,
            self.interval_size,
            motif_lookup,
            combine_strands,
            position_filter,
        )?;

        let tid_progress = master_progress
            .add(get_master_progress_bar_fancy(feeder.total_length()));
        tid_progress.set_message("genome positions");
        let write_progress = master_progress.add(get_ticker());
        write_progress.set_message("rows written");
        let erred_reads = master_progress.add(get_ticker());
        erred_reads.set_message("~records errored");

        let (jobs_tx, jobs_rx) = crossbeam_channel::bounded(workers.len() * 2);
        let (results_tx, results_rx) = crossbeam_channel::unbounded();
        let (records_tx, records_rx) = crossbeam_channel::bounded(
            self.queue_size.unwrap_or(workers.len() * 2),
        );

        for _ in 0..(workers.len() * 2) {
            empties_tx.send(ModBasePileup2::new_empty()).unwrap();
        }
        let mpb_handle = master_progress.clone();

        let source = std::thread::spawn({
            let results_handle = results_tx.clone();
            let feeder = feeder
                .into_iter()
                .inspect(move |r| match r {
                    Ok(_) => {}
                    Err(e) => {
                        mpb_handle.suspend(|| {
                            error!("failed to fetch sequence, {e}")
                        });
                    }
                })
                .filter_map(|r| r.ok());
            move || {
                let get_mem = || -> Result<ModBasePileup2, ()> {
                    match empties_rx.recv() {
                        Ok(mem) => Ok(mem),
                        Err(_) => Err(()),
                    }
                };
                let mut seq = 0usize;
                for chrom_coords in feeder {
                    match get_mem() {
                        Ok(mem) => {
                            if jobs_tx.send((seq, chrom_coords, mem)).is_err() {
                                break;
                            }
                            seq = seq.wrapping_add(1);
                        }
                        Err(_) => {
                            break;
                        }
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
                while let Ok((seq, item, mem)) = jobs_rx.recv() {
                    let out = worker_state.process(item, mem);
                    // If collector disappeared, we can exit.
                    if results_tx.send((seq, out)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(results_tx);
        drop(jobs_rx);

        let aggregator = std::thread::spawn(move || {
            let mut next_seq = 0usize;
            let mut buffer: BTreeMap<usize, anyhow::Result<ModBasePileup2>> =
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

        let mut n_failed_intervals = 0usize;
        for result in records_rx.into_iter() {
            match result {
                Ok(mod_base_pileup) => {
                    tid_progress.inc(mod_base_pileup.interval_width as u64);
                    erred_reads.inc(mod_base_pileup.failed_records as u64);
                    let rows_written =
                        writer.write(mod_base_pileup, &motif_labels)?;
                    write_progress.inc(rows_written);
                }
                Err(message) => {
                    n_failed_intervals += 1;
                    master_progress.suspend(|| {
                        error!("failed to process interval, {message}");
                    });
                }
            }
        }

        let rows_processed = write_progress.position();
        let n_failed_reads = erred_reads.position();

        if n_failed_reads > 0 {
            master_progress.suspend(|| {
                error!("~{n_failed_reads} records failed processing");
            });
        }
        if n_failed_intervals > 0 {
            master_progress.suspend(|| {
                error!(
                    "{n_failed_intervals} interval(s) failed processing, \
                     output is incomplete"
                );
            });
        }

        write_progress.finish_and_clear();
        erred_reads.finish_and_clear();
        master_progress.suspend(|| {
            info!("Done, processed {rows_processed} rows.");
            crate::mod_bam::report_conflict_summary();
        });
        tid_progress.finish_and_clear();

        source.join().expect("source thread paniced");
        for (i, worker_thread) in handles.into_iter().enumerate() {
            worker_thread.join().expect(&format!("worker thread {i} paniced"));
        }
        aggregator.join().expect("aggregator theread paniced");

        Ok(())
    }
}

/// Check the first `n_records` primary records of a modBAM for base
/// modification calls on the opposite strand of the primary sequence base
/// (`-` in the MM tag, e.g. PacBio `T-a.`).
fn bam_has_negative_strand_mods(
    in_bam: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    n_records: usize,
) -> anyhow::Result<bool> {
    let mut reader = bam::Reader::from_path(in_bam)?;
    if unindexed_reader_is_cram(&reader) {
        if let Some(reference_fp) = reference_fasta {
            reader.set_reference(reference_fp)?;
        } else {
            bail!("CRAM input requires reference")
        }
    }
    let mut record = bam::Record::new();
    let mut n_checked = 0usize;
    while let Some(result) = reader.read(&mut record) {
        result?;
        if record_is_not_primary(&record) {
            continue;
        }
        if let Ok(mm_tag_infos) = MmTagInfo::from_record(&record) {
            if mm_tag_infos.iter().any(|x| x.strand() == Strand::Negative) {
                return Ok(true);
            }
        }
        n_checked += 1;
        if n_checked >= n_records {
            break;
        }
    }
    Ok(false)
}

#[derive(Clone, Debug)]
enum Presets {
    /// CpG-special, combine strands, maybe combine mods
    DnaCpGCombineStrands { cytosine_other_mod: Option<ModCodeRepr> },
    /// Cytosine-special, not combine strands, but combine mods
    DnaCytosineCombine,
    /// Everything else
    DnaAllContext {
        cytosine_mod_op: Option<ModCodeRepr>,
        motif_bases: [DnaBase; 4],
    },
    DynamicAllContext {
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_bases: [DnaBase; 4],
        motif_offset: Option<u32>,
    },
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct DuplexModBamPileup {
    // running args
    /// Input BAM, should be sorted and have associated index available.
    in_bam: PathBuf,
    /// Output file to write results into. Will write to stdout if not
    /// provided.
    #[arg(short = 'o', long)]
    out_bed: Option<PathBuf>,
    /// Aggregate double-stranded base modifications for CpG dinucleotides.
    /// This flag is short-hand for --motif CG 0.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, group = "motif_options", default_value_t = false)]
    cpg: bool,

    /// Specify the sequence motif to pileup double-stranded base modification
    /// pattern counts for. The first argument should be the sequence motif
    /// and the second argument is the 0-based offset to the base to pileup
    /// base modification counts for. For example: --motif CG 0 indicates
    /// to generate pattern counts for the C on the top strand
    /// and the following C (opposite to G) on the negative strand. The motif
    /// must be reverse-complement palindromic or an error will be raised.
    /// See the documentation for more examples and details.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, group = "motif_options", num_args = 2)]
    motif: Option<Vec<String>>,
    /// Reference sequence in FASTA format.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long = "ref", alias = "reference", short = 'r')]
    reference_fasta: PathBuf,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Process only the specified region of the BAM when performing pileup.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>. Commas are
    /// allowed.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,
    /// Maximum number of records to use when calculating pileup. This argument
    /// is passed to the pileup engine. If you have high depth data,
    /// consider increasing this value substantially. Must be less than
    /// 2147483647 or an error will be raised.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, default_value_t = 8000, hide_short_help = true)]
    max_depth: u32,

    // processing args
    /// Number of threads to use while processing chunks concurrently.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
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
    /// Size of queue for writing records
    #[clap(help_heading = "Compute Options")]
    #[arg(long, hide_short_help = true, default_value_t = 1000)]
    queue_size: usize,

    /// Break contigs into chunks containing this many intervals (see
    /// `interval_size`). This option can be used to help prevent excessive
    /// memory usage, usually with no performance penalty. By default,
    /// modkit will set this value to 1.5x the number of threads specified,
    /// so if 4 threads are specified the chunk_size will be 6.
    /// A warning will be shown if this option is less than the number of
    /// threads specified.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, hide_short_help = true)]
    chunk_size: Option<usize>,
    #[clap(help_heading = "Logging Options")]
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

    // sampling args
    /// Sample this many reads when estimating the filtering threshold. Reads
    /// will be sampled evenly across aligned genome. If a region is
    /// specified, either with the --region option or the --sample-region
    /// option, then reads will be sampled evenly across the region given.
    /// This option is useful for large BAM files. In practice, 10-50
    /// thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'n',
        long,
        default_value_t = 10_042
    )]
    num_reads: usize,
    /// Sample this fraction of the reads when estimating the
    /// filter-percentile. In practice, 50-100 thousand reads is sufficient
    /// to estimate the model output distribution and determine the
    /// filtering threshold. See filtering.md for details on filtering.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        group = "sampling_options",
        short = 'f',
        long,
        hide_short_help = true
    )]
    sampling_frac: Option<f64>,
    /// Set a random seed for deterministic running, the default is
    /// non-deterministic.
    #[clap(help_heading = "Sampling Options")]
    #[arg(
        long,
        conflicts_with = "num_reads",
        requires = "sampling_frac",
        hide_short_help = true
    )]
    seed: Option<u64>,
    /// Do not perform any filtering, include all mod base calls in output. See
    /// filtering.md for details on filtering.
    #[clap(help_heading = "Filtering Options")]
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or per-base. Global filter
    /// threshold can be specified with by a decimal number (e.g. 0.75).
    /// Per-base thresholds can be specified by colon-separated values, for
    /// example C:0.75 specifies a threshold value of 0.75 for cytosine
    /// modification calls. Additional per-base thresholds can be specified
    /// by repeating the option: for example --filter-threshold C:0.75
    /// --filter-threshold A:0.70 or specify a single base option and a
    /// default for all other bases with: --filter-threshold A:0.70
    /// --filter-threshold 0.9 will specify a threshold value of 0.70 for
    /// adenine and 0.9 for all other base modification calls.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long,
    group = "thresholds",
    action = clap::ArgAction::Append,
    alias = "pass_threshold"
    )]
    filter_threshold: Option<Vec<String>>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[clap(help_heading = "Filtering Options")]
    #[arg(
    long = "mod-threshold",
    alias = "mod-thresholds",
    action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Specify a region for sampling reads from when estimating the threshold
    /// probability. If this option is not provided, but --region is
    /// provided, the genomic interval passed to --region will be used.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[clap(help_heading = "Filtering Options")]
    #[arg(long)]
    sample_region: Option<String>,
    /// Interval chunk size in base pairs to process concurrently when
    /// estimating the threshold probability, can be larger than the pileup
    /// processing interval.
    #[clap(help_heading = "Filtering Options")]
    #[arg(long, default_value_t = 1_000_000, hide_short_help = true)]
    sampling_interval_size: u32,
    /// BED file that will restrict threshold estimation and pileup results to
    /// positions overlapping intervals in the file. (alias: include-positions)
    #[clap(help_heading = "Selection Options")]
    #[arg(long, hide_short_help = true, alias = "include-positions")]
    include_bed: Option<PathBuf>,
    /// Include unmapped base modifications when estimating the pass threshold.
    #[clap(help_heading = "Selection Options")]
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = false,
        conflicts_with = "include_bed"
    )]
    include_unmapped: bool,

    // collapsing and combining args
    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(long, group = "combine_args", hide_short_help = true)]
    ignore: Option<String>,
    /// Force allow implicit-canonical mode. By default modkit does not allow
    /// pileup with the implicit mode (e.g. C+m, no '.' or '?'). The
    /// `update-tags` subcommand is provided to update tags to the new
    /// mode. This option allows the interpretation of implicit mode tags:
    /// residues without modified base probability will be interpreted as
    /// being the non-modified base.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = false,
        hide_short_help = true
    )]
    force_allow_implicit: bool,

    /// Respect soft masking in the reference FASTA.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(
        long,
        short = 'k',
        requires = "reference_fasta",
        default_value_t = false,
        hide_short_help = true
    )]
    mask: bool,
    /// Combine base modification calls, all counts of modified bases are
    /// summed together. See collapse.md for details.
    #[clap(help_heading = "Modified Base Options")]
    #[arg(
        long,
        default_value_t = false,
        group = "combine_args",
        hide_short_help = true
    )]
    combine_mods: bool,
    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, hide_short_help = true)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[clap(help_heading = "Selection Options")]
    #[arg(
        long,
        requires = "edge_filter",
        default_value_t = false,
        hide_short_help = true
    )]
    invert_edge_filter: bool,

    // output args
    /// **Deprecated** The default output has all tab-delimiters.
    /// For bedMethyl output, separate columns with only tabs. The default is
    /// to use tabs for the first 10 fields and spaces thereafter. The
    /// default behavior is more likely to be compatible with genome viewers.
    /// Enabling this option may make it easier to parse the output with
    /// tabular data handlers that expect a single kind of separator.
    #[clap(help_heading = "Output Options")]
    #[arg(
        long,
        hide_short_help = true,
        conflicts_with = "mixed_delimiters",
        default_value_t = false
    )]
    only_tabs: bool,
    /// Output bedMethyl where the delimiter of columns past column 10 are
    /// space-delimited instead of tab-delimited. This option can be useful
    /// for some browsers and parsers that don't expect the extra columns
    /// of the bedMethyl format.
    #[clap(help_heading = "Output Options")]
    #[arg(
        long = "mixed-delim",
        conflicts_with = "only_tabs",
        alias = "mixed-delimiters",
        default_value_t = false,
        hide_short_help = true
    )]
    mixed_delimiters: bool,
}

impl DuplexModBamPileup {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if self.only_tabs {
            warn!(
                "--only-tabs is deprecated. The default output format will \
                 have only tabs. In the next version, passing this flag will \
                 cause an error. To get the old version of the format use \
                 --mixed-delim"
            );
        }
        // do this first so we fail when the file isn't readable
        let header =
            bam::IndexedReader::from_path(&self.in_bam).map(|reader| {
                if !reader_is_bam(&reader) {
                    info!(
                        "\
                    detected non-BAM input format, please consider using BAM, \
                         CRAM may be unstable"
                    );
                }
                reader.header().to_owned()
            })?;

        // options parsing below
        let region = self
            .region
            .as_ref()
            .map(|raw_region| {
                info!("parsing region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let sampling_region = self
            .sample_region
            .as_ref()
            .map(|raw_region| {
                info!("parsing sample region {raw_region}");
                Region::parse_str(raw_region, &header)
            })
            .transpose()?;
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;
        let per_mod_thresholds = self
            .mod_thresholds
            .as_ref()
            .map(|raw_per_mod_thresholds| {
                parse_per_mod_thresholds(raw_per_mod_thresholds)
            })
            .transpose()?;
        let reference_records = get_targets(&header, region.as_ref());
        let position_filter = self
            .include_bed
            .as_ref()
            .map(|bed_fp| {
                let chrom_to_tid = reference_records
                    .iter()
                    .map(|reference_record| {
                        (reference_record.name.as_str(), reference_record.tid)
                    })
                    .collect::<HashMap<&str, u32>>();
                StrandedPositionFilter::from_bed_file(
                    bed_fp,
                    &chrom_to_tid,
                    self.suppress_progress,
                )
            })
            .transpose()?;
        // use the path here instead of passing the reader directly to avoid
        // potentially changing mutable internal state of the reader.
        IdxStats::check_any_mapped_reads(
            &self.in_bam,
            region.as_ref(),
            position_filter.as_ref(),
        )
        .context(
            "\
            did not find any mapped reads, perform alignment first or use \
             modkit extract and/or modkit summary to inspect unaligned modBAMs",
        )?;
        let chunk_size = if let Some(chunk_size) = self.chunk_size {
            if chunk_size < self.threads {
                warn!(
                    "chunk size {chunk_size} is less than number of threads \
                     ({}), this will limit parallelism",
                    self.threads
                );
            }
            chunk_size
        } else {
            let cs = (self.threads as f32 * 1.5).floor() as usize;
            info!(
                "calculated chunk size: {cs}, interval size {}, processing {} \
                 positions concurrently",
                self.interval_size,
                cs * self.interval_size as usize
            );
            cs
        };

        if self.filter_percentile > 1.0 {
            bail!("filter percentile must be <= 1.0")
        }
        let (pileup_options, collapse_method) =
            match (self.combine_mods, self.ignore.as_ref()) {
                (false, None) => (PileupNumericOptions::Passthrough, None),
                (true, _) => (PileupNumericOptions::Combine, None),
                (_, Some(raw_mod_code)) => {
                    let mod_code = ModCodeRepr::parse(&raw_mod_code)?;
                    info!("ignoring mod code {}", raw_mod_code);
                    let method = CollapseMethod::ReDistribute(mod_code);
                    (
                        PileupNumericOptions::Collapse(method.clone()),
                        Some(method),
                    )
                }
            };

        // motif handling
        let regex_motif = {
            if self.cpg {
                RegexMotif::parse_string("CG", 0)?
            } else {
                if self.motif.is_none() {
                    bail!(
                        "either --cpg or a --motif must be provided for \
                         pileup-hemi"
                    )
                }
                let raw_motif = self.motif.as_ref().unwrap();
                if raw_motif.len() != 2 {
                    bail!("motif arg should be length 2, eg. CG 0")
                }
                let motif_seq = &raw_motif[0];
                let focus_base =
                    raw_motif[1].parse::<usize>().map_err(|e| {
                        anyhow!(
                            "couldn't parse focus position, {}",
                            e.to_string()
                        )
                    })?;
                let regex_motif =
                    RegexMotif::parse_string(motif_seq, focus_base)?;
                regex_motif
            }
        };
        if !regex_motif.is_palendrome() {
            bail!("motif must be palindromic for pileup-hemi")
        }

        let mut writer: Box<dyn PileupWriter<DuplexModBasePileup>> =
            if let Some(out_fp) = self.out_bed.as_ref() {
                create_out_directory(out_fp)?;
                let fh = std::fs::File::create(out_fp)
                    .context("failed to make output file")?;
                let writer = BufWriter::new(fh);
                Box::new(BedMethylWriter::new(
                    writer,
                    self.mixed_delimiters,
                    false,
                )?)
            } else {
                let writer = BufWriter::new(std::io::stdout());
                Box::new(BedMethylWriter::new(
                    writer,
                    self.mixed_delimiters,
                    false,
                )?)
            };

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .with_context(|| "failed to make threadpool")?;

        // put this into it's own function
        let motif_info = regex_motif.motif_info;
        let motif_lookup = MotifLocationsLookup::from_paths(
            &self.reference_fasta,
            self.mask,
            None,
            vec![regex_motif],
            false,
        )?;

        // start the actual work here
        let threshold_caller =
            if let Some(raw_threshold) = &self.filter_threshold {
                parse_thresholds(raw_threshold, per_mod_thresholds)?
            } else {
                pool.install(|| {
                    get_threshold_from_options(
                        &self.in_bam,
                        self.threads,
                        self.sampling_interval_size,
                        self.sampling_frac,
                        self.num_reads,
                        self.no_filtering,
                        self.filter_percentile,
                        self.seed,
                        sampling_region.as_ref().or(region.as_ref()),
                        per_mod_thresholds,
                        edge_filter.as_ref(),
                        collapse_method.as_ref(),
                        position_filter.as_ref(),
                        !self.include_unmapped,
                        self.suppress_progress,
                    )
                })?
            };

        if !self.no_filtering {
            for (base, threshold) in threshold_caller.iter_thresholds() {
                let base = base.char();
                match (threshold * 100f32).ceil() as usize {
                    0..=60 => error!(
                        "Threshold of {threshold} for base {base} is very \
                         low. Consider increasing the filter-percentile or \
                         specifying a higher threshold."
                    ),
                    61..=70 => warn!(
                        "Threshold of {threshold} for base {base} is low. \
                         Consider increasing the filter-percentile or \
                         specifying a higher threshold."
                    ),
                    _ => info!(
                        "Using filter threshold {} for {base}.",
                        threshold
                    ),
                }
            }
            for (mod_code_repr, threshold) in
                threshold_caller.iter_mod_thresholds()
            {
                match (threshold * 100f32).ceil() as usize {
                    0..=60 => error!(
                        "Threshold of {threshold} for mod code \
                         {mod_code_repr} is very low. Consider increasing the \
                         filter-percentile or specifying a higher threshold."
                    ),
                    61..=70 => warn!(
                        "Threshold of {threshold} for mod code \
                         {mod_code_repr} is low. Consider increasing the \
                         filter-percentile or specifying a higher threshold."
                    ),
                    _ => info!(
                        "Using filter threshold {} for mod code \
                         {mod_code_repr}.",
                        threshold
                    ),
                }
            }
        }

        let (snd, rx) = bounded(self.queue_size);
        let reference_records = if let Some(pf) = position_filter.as_ref() {
            pf.optimize_reference_records(reference_records, self.interval_size)
        } else {
            reference_records
        };
        let feeder = ReferenceIntervalBatchesFeeder::new(
            reference_records,
            chunk_size,
            self.interval_size,
            true, // must be true for duplex
            Some(motif_lookup),
            position_filter,
        )?;

        let in_bam_fp = self.in_bam.clone();

        let master_progress = MultiProgress::new();
        if self.suppress_progress {
            master_progress
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        let tid_progress = master_progress
            .add(get_master_progress_bar(feeder.total_length() as usize));
        tid_progress.set_message("genome positions");
        let write_progress = master_progress.add(get_ticker());
        write_progress.set_message("rows written");
        let skipped_reads = master_progress.add(get_ticker());
        skipped_reads.set_message("~records skipped");
        let processed_reads = master_progress.add(get_ticker());
        processed_reads.set_message("~records processed");

        let force_allow = self.force_allow_implicit;
        let max_depth = self.max_depth;

        pool.spawn(move || {
            for multi_chrom_coords in feeder
                .inspect(|x| match x {
                    Ok(_) => {},
                    Err(e) => {
                        error!("fetching sequence failed, {e}");
                    }
                })
                .filter_map(|r| r.ok()) {
                let genome_length_in_batch = multi_chrom_coords.iter()
                    .map(|x| x.total_length())
                    .sum::<u64>();
                let n_intervals = multi_chrom_coords.len();
                let interval_progress = master_progress
                    .add(get_subroutine_progress_bar(n_intervals));

                interval_progress
                    .set_message(format!("processing {n_intervals} intervals"));
                for work_chunk in multi_chrom_coords.chunks(chunk_size) {
                    let mut result: Vec<anyhow::Result<DuplexModBasePileup>> = vec![];
                    let chunk_progress = master_progress.add(get_subroutine_progress_bar(work_chunk.len()));
                    chunk_progress.set_message("chunk progress");
                    let (res, _) = rayon::join(
                        || {
                            work_chunk
                                .into_par_iter()
                                .progress_with(chunk_progress)
                                .map(|multi_chrom_coords| {
                                    process_region_duplex_batch(
                                        multi_chrom_coords,
                                        &in_bam_fp,
                                        &threshold_caller,
                                        &pileup_options,
                                        force_allow,
                                        max_depth,
                                        motif_info,
                                        edge_filter.as_ref(),
                                    )
                                })
                                .flatten()
                                .collect::<Vec<anyhow::Result<DuplexModBasePileup>>>()
                        },
                        || {
                            result.into_iter().for_each(|mod_base_pileup| {
                                match snd.send(mod_base_pileup) {
                                    Ok(_) => {
                                        interval_progress.inc(1)
                                    }
                                    Err(e) => {
                                        error!("failed to send results, {}", e.to_string())
                                    },
                                }
                            });
                        },
                    );
                    result = res;
                    result.into_iter().for_each(|pileup| {
                        match snd.send(pileup) {
                            Ok(_) => {
                                interval_progress.inc(1)
                            }
                            Err(e) => {
                                error!("failed to send results, {}", e.to_string())
                            },
                        }
                    });
                }
                tid_progress.inc(genome_length_in_batch);
            }
            tid_progress.finish_and_clear();
        });

        for result in rx.into_iter() {
            match result {
                Ok(mod_base_pileup) => {
                    processed_reads
                        .inc(mod_base_pileup.processed_records as u64);
                    skipped_reads.inc(mod_base_pileup.skipped_records as u64);
                    let rows_written = writer.write(mod_base_pileup, &[])?;
                    write_progress.inc(rows_written);
                }
                Err(message) => {
                    debug!("> unexpected error {message}");
                }
            }
        }
        let rows_processed = write_progress.position();
        let n_skipped_reads = skipped_reads.position();
        let n_skipped_message = if n_skipped_reads == 0 {
            format!("zero reads")
        } else {
            format!("~{n_skipped_reads} reads")
        };
        let n_processed_reads = processed_reads.position();
        write_progress.finish_and_clear();
        processed_reads.finish_and_clear();
        skipped_reads.finish_and_clear();
        info!(
            "Done, processed {rows_processed} rows. Processed \
             ~{n_processed_reads} reads and skipped {n_skipped_message}."
        );
        crate::mod_bam::report_conflict_summary();
        Ok(())
    }
}
