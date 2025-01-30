use std::string::FromUtf8Error;

pub type MkResult<T, E = MkError> = Result<T, E>;

#[derive(thiserror::Error, Debug)]
pub enum MkError {
    // Tag parsing
    #[error("invalid-MM-tag")]
    InvalidMm(String),
    #[error("invalid-ML-tag")]
    InvalidMl(String),
    #[error("MM-tag-missing")]
    MmMissing,
    #[error("ML-tag-missing")]
    MlMissing,
    #[error("invalid-MN-tag")]
    InvalidMn(String),
    #[error("invalid-MM-mode")]
    InvalidSkipMode,
    #[error("non-primary-no-MN")]
    NonPrimaryMissingMn,
    #[error("aux-data-missing")]
    AuxMissing,
    #[error("multiple-tag-instances")]
    MultipleTagInstances,
    #[error("conflict-{}", .0)]
    Conflict(#[from] ConflictError),
    #[error("HtsLib-error-{}", .0)]
    HtsLibError(#[from] rust_htslib::errors::Error),
    #[error("no-modbase-info")]
    NoModifiedBaseInformation,

    // DNA/RNA "type" symbol parsing
    #[error("invalid-DNA-RNA-base")]
    InvalidDnaBase,
    #[error("invalid-strand")]
    InvalidStrand,

    // Pileup-specific
    #[error("invalid-implicit-mode")]
    InvalidImplicitMode,

    // Adjust
    #[error("invalid-collapse-method")]
    InvalidCollapseMethod,

    // Misc
    #[error("invalid-record-name")]
    InvalidRecordName,
    #[error("invalid-cigar")]
    InvalidCigar,
    #[error("invalid-read-sequence")]
    InvalidReadSequence(#[from] FromUtf8Error),
    #[error("empty-read-sequence")]
    EmptyReadSequence,
    #[error("invalid region, {}, should be 'chrom' or 'chrom:start-stop'", .0)]
    InvalidRegion(String),
    #[error("contig-missing")]
    ContigMissing,

    // Entropy
    #[error("zero-reads")]
    EntropyZeroCoverage { chrom_id: u32, start: u64, end: u64 },
    #[error("insufficient-coverage")]
    EntropyInsufficientCoverage { chrom_id: u32, start: u64, end: u64 },

    // Maths
    #[error("not enough datapoints, got {}", .0)]
    PercentileNotEnoughDatapoints(usize),
    #[error("invalid quantile, got {}", .0)]
    PercentileInvalidQuantile(f32),
}

#[derive(thiserror::Error, Debug)]
pub enum ConflictError {
    #[error("inferred-prob-greater-than-one")]
    InferredSumGreaterThanOne,
    #[error("explicit-prob-greater-than-one")]
    ProbaGreaterThanOne,
    #[error("explicit-and-inferred")]
    ExplicitConflictInferred,
}
