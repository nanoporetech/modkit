use clap::Args;
use std::path::PathBuf;

#[derive(Args)]
pub(super) struct InputArgs {
    /// Path to modBAM file to extract read-level information from, or one of
    /// `-` or `stdin` to specify a stream from standard input. If a file
    /// is used it may be sorted and have associated index.
    pub in_bam: String,
    /// Path to output file, "stdout" or "-" will direct output to standard
    /// out.
    pub out_path: String,
    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,
    /// Write output as BGZF compressed file.
    #[arg(long, default_value_t = false)]
    pub bgzf: bool,
    /// Number of threads to use for parallel bgzf writing.
    #[arg(long, requires = "bgzf", default_value_t = 4)]
    pub out_threads: usize,

    /// Number of reads that can be in memory at a time. Increasing this value
    /// will increase thread usage, at the cost of memory usage.
    #[arg(short = 'q', long, default_value_t = 10_000)]
    pub queue_size: usize,
    /// Path to file to write run log.
    #[arg(long, alias = "log")]
    pub log_filepath: Option<PathBuf>,
    /// Include only mapped bases in output (alias: mapped).
    #[arg(long, alias = "mapped", default_value_t = false)]
    pub mapped_only: bool,
    /// Output aligned secondary and supplementary base modification
    /// probabilities as additional rows. The primary alignment will have
    /// all of the base modification probabilities (including soft-clipped
    /// ones, unless --mapped-only is used). The non-primary alignments
    /// will only have mapped bases in the output.
    #[arg(long, alias = "non-primary", default_value_t = false)]
    pub allow_non_primary: bool,
    /// Number of reads to use. Note that when using a sorted, indexed modBAM
    /// that the sampling algorithm will attempt to sample records evenly
    /// over the length of the reference sequence. The result is the final
    /// number of records used may be slightly more or less than the
    /// requested number. When piping from stdin or using a modBAM without
    /// an index, the requested number of reads will be the first `num_reads`
    /// records.
    #[arg(long)]
    pub num_reads: Option<usize>,
    /// Process only reads that are aligned to a specified region of the BAM.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    pub region: Option<String>,
    /// Force overwrite of output file
    #[arg(long, default_value_t = false)]
    pub force: bool,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    pub suppress_progress: bool,
    /// Set the query and reference k-mer size (if a reference is provided).
    /// Maximum number for this value is 50.
    #[arg(long, default_value_t = 5)]
    pub kmer_size: usize,
    /// Ignore the BAM index (if it exists) and default to a serial scan of the
    /// BAM.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    pub ignore_index: bool,

    /// Don't print the header lines in the output tables.
    #[arg(long, default_value_t = false)]
    pub no_headers: bool,

    /// BED file with regions to include (alias: include-positions). Implicitly
    /// only includes mapped sites.
    #[arg(long, alias = "include-positions")]
    pub include_bed: Option<PathBuf>,
    /// BED file with regions to _exclude_ (alias: exclude).
    #[arg(long, alias = "exclude", short = 'v')]
    pub exclude_bed: Option<PathBuf>,
    /// Output read-level base modification probabilities restricted to the
    /// reference sequence motifs provided. The first argument should be
    /// the sequence motif and the second argument is the 0-based offset to
    /// the base to pileup base modification counts for. For example:
    /// --motif CGCG 0 indicates include base modifications for which the read
    /// is aligned to the first C on the top strand and the last C
    /// (complement to G) on the bottom strand. The --cpg argument is short
    /// hand for --motif CG 0. This argument can be passed multiple times.
    #[arg(long, action = clap::ArgAction::Append, num_args = 2, requires = "reference")]
    pub motif: Option<Vec<String>>,
    /// Only output counts at CpG motifs. Requires a reference sequence to be
    /// provided.
    #[arg(long, requires = "reference", default_value_t = false)]
    pub cpg: bool,
    /// When using motifs, respect soft masking in the reference sequence.
    #[arg(
        long,
        short = 'k',
        requires = "motif",
        default_value_t = false,
        hide_short_help = true
    )]
    pub mask: bool,

    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    pub edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    pub invert_edge_filter: bool,

    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[arg(long, hide_short_help = true)]
    pub ignore: Option<String>,

    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    /// Only used when an indexed modBAM is provided.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    pub interval_size: u32,

    /// Ignore implicitly canonical base modification calls. When the `.`
    /// flag is used in the MM tag, this implies that bases missing a base
    /// modification probability are to be assumed canonical. Set this flag
    /// to omit those base modifications from the output. For additional
    /// details see the SAM spec: https://samtools.github.io/hts-specs/SAMtags.pdf.
    #[arg(long, hide_short_help = true)]
    pub ignore_implicit: bool,
}
