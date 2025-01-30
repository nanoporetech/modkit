use clap::{arg, Args};
use std::path::PathBuf;

#[derive(Args, Clone, Debug)]
pub(super) struct InputArgs {
    /// Input bedmethyl table, can be used directly from modkit pileup.
    #[arg(short = 'i', long)]
    pub in_bedmethyl: PathBuf,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
    /// Number of tabix/bgzf IO threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 2, hide_short_help = true)]
    pub io_threads: usize,
    /// Reference sequence in FASTA format used for the pileup.
    #[arg(short = 'r', long = "ref")]
    pub reference_fasta: PathBuf,
    /// Use only bedMethyl records from this contig, requires that the
    /// bedMethyl be BGZIP-compressed and tabix-indexed.
    #[arg(long)]
    pub contig: Option<String>,
    /// Output log to this file.
    #[arg(long, alias = "log")]
    #[clap(help_heading = "Logging Options")]
    pub log_filepath: Option<PathBuf>,
    /// Disable the progress bars.
    #[arg(long, default_value_t = false)]
    #[clap(help_heading = "Logging Options")]
    pub suppress_progress: bool,
}

#[derive(Args, Clone, Debug)]
pub(super) struct RefineArgs {
    /// Fraction modified threshold below which consider a genome location to
    /// be "low modification".
    #[clap(help_heading = "Search Options")]
    #[arg(long = "low-thresh", default_value_t = 0.2)]
    pub low_threshold: f32,
    /// Fraction modified threshold above which consider a genome location to
    /// be "high modification" or enriched for modification.
    #[clap(help_heading = "Search Options")]
    #[arg(long = "high-thresh", default_value_t = 0.6)]
    pub high_threshold: f32,
    /// Minimum log-odds to consider a motif sequence to be enriched.
    #[clap(help_heading = "Search Options")]
    #[arg(long, default_value_t = 1.5)]
    pub min_log_odds: f32,
    /// Minimum log-odds to consider a motif sequence to be enriched when
    /// performing exhaustive search.
    #[clap(help_heading = "Search Options")]
    #[arg(long, conflicts_with = "skip_search", default_value_t = 2.5)]
    pub exhaustive_seed_min_log_odds: f32,
    /// Exhaustive search seed length, increasing this value increases
    /// computational time.
    #[clap(help_heading = "Search Options")]
    #[arg(long, conflicts_with = "skip_search", default_value_t = 3)]
    pub exhaustive_seed_len: usize,
    /// Skip the exhaustive search phase, saves time but the results may be
    /// less sensitive
    #[clap(help_heading = "Search Options")]
    #[arg(long, default_value_t = false)]
    pub skip_search: bool,
    /// Minimum coverage in the bedMethyl to consider a record valid.
    #[clap(help_heading = "Search Options")]
    #[arg(long, default_value_t = 5)]
    pub min_coverage: u64,
    /// Upstream and downstream number of bases to search for a motif sequence
    /// around a modified base. Example: --context-size 12 12.
    #[clap(help_heading = "Search Options")]
    #[arg(long, num_args=2, default_values_t=vec![12, 12])]
    pub context_size: Vec<u64>,
    /// Minimum number of total sites in the genome required for a motif to be
    /// considered.
    #[clap(help_heading = "Search Options")]
    #[arg(long, default_value_t = 300)]
    pub min_sites: u64,
    /// Minimum fraction of sites in the genome to be "high-modification" for a
    /// motif to be considered.
    #[clap(help_heading = "Search Options")]
    #[arg(long = "min-frac-mod", default_value_t = 0.85)]
    pub frac_sites_thresh: f32,
}

#[derive(Args, Clone, Debug)]
pub(super) struct KnownMotifsArgs {
    /// Format should be <sequence> <offset> <mod_code>.
    #[clap(help_heading = "Output Options")]
    #[arg(long="known-motif", num_args = 3, action = clap::ArgAction::Append)]
    pub known_motifs: Option<Vec<String>>,
    /// Path to known motifs in tabular format. Tab-separated values:
    /// <mod_code>\t<motif_seq>\t<offset>. May have the same header as the
    /// output table from this command.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "known-motifs-table")]
    pub known_motifs_table: Option<PathBuf>,
}
