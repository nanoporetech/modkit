use crate::mod_bam::{CollapseMethod, EdgeFilter};
use crate::monoid::Moniod;
use crate::reads_sampler::record_sampler::RecordSampler;
use rust_htslib::bam;

pub(crate) trait RecordProcessor {
    type Output: Moniod + WithRecords + Send;

    fn process_records<T: bam::Read>(
        records: bam::Records<T>,
        with_progress: bool,
        record_sampler: RecordSampler,
        collapse_method: Option<&CollapseMethod>,
        edge_filter: Option<&EdgeFilter>,
    ) -> anyhow::Result<Self::Output>;
}

pub(crate) trait WithRecords {
    /// Number of rows/records in this datastructure.
    /// For example, if this contains mod probs for
    /// a set of reads, total number of mod probs
    /// for all of the reads combined. May be the same
    /// number as `num_reads`.
    fn size(&self) -> u64;
    /// Number of reads in this data structure.
    fn num_reads(&self) -> usize;
    // todo(arand) consider adding num_failed?
}
