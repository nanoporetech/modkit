use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use derive_new::new;
use log::debug;
use rust_htslib::bam;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::mod_bam::CollapseMethod;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::util::StrandRule;

pub(crate) mod base_mods_adapter;
pub(crate) mod bedrmod;
pub(crate) mod duplex;
pub(super) mod pileup_processor;
pub mod subcommand;

#[derive(Debug, Copy, Clone, new)]
pub struct PileupFeatureCounts2 {
    pub position: u32,
    pub raw_strand: char,
    pub filtered_coverage: u16,
    pub mod_code: ModCodeRepr,
    pub n_canonical: u16,
    pub n_modified: u16,
    pub n_other_modified: u16,
    pub n_delete: u16,
    pub n_filtered: u16,
    pub n_diff: u16,
    pub n_nocall: u16,
    pub motif_idxs: u8,
}

impl PileupFeatureCounts2 {
    pub(super) fn is_valid(&self) -> bool {
        self.filtered_coverage > 0
    }

    pub(super) fn add(&mut self, other: &Self) {
        assert_eq!(self.position, other.position, "{:?} =/= {:?}", self, other);
        assert_eq!(self.raw_strand, other.raw_strand);
        assert_eq!(self.mod_code, other.mod_code);
        self.filtered_coverage =
            self.filtered_coverage.saturating_add(other.filtered_coverage);
        self.n_canonical = self.n_canonical.saturating_add(other.n_canonical);
        self.n_modified = self.n_modified.saturating_add(other.n_modified);
        self.n_other_modified =
            self.n_other_modified.saturating_add(other.n_other_modified);
        self.n_delete = self.n_delete.saturating_add(other.n_delete);
        self.n_filtered = self.n_filtered.saturating_add(other.n_filtered);
        self.n_diff = self.n_diff.saturating_add(other.n_diff);
        self.n_nocall = self.n_nocall.saturating_add(other.n_nocall);
    }
}

#[derive(Debug, Copy, Clone, new)]
pub struct PileupFeatureCounts {
    pub raw_strand: char,
    pub filtered_coverage: u32,
    pub raw_mod_code: ModCodeRepr,
    pub fraction_modified: f32,
    pub n_canonical: u32,
    pub n_modified: u32,
    pub n_other_modified: u32,
    pub n_delete: u32,
    pub n_filtered: u32,
    pub n_diff: u32,
    pub n_nocall: u32,
    pub motif_idx: Option<usize>,
}

impl From<BedMethylLine> for PileupFeatureCounts {
    fn from(item: BedMethylLine) -> PileupFeatureCounts {
        PileupFeatureCounts::new(
            item.strand.into(),
            item.valid_coverage.try_into().unwrap_or(0),
            item.raw_mod_code,
            item.frac_modified(),
            item.count_canonical.try_into().unwrap_or(0),
            item.count_methylated.try_into().unwrap_or(0),
            item.count_other.try_into().unwrap_or(0),
            item.count_delete.try_into().unwrap_or(0),
            item.count_fail.try_into().unwrap_or(0),
            item.count_diff.try_into().unwrap_or(0),
            item.count_nocall.try_into().unwrap_or(0),
            None,
        )
    }
}

#[derive(new)]
struct StrandPileup {
    pub(crate) bam_pileup: bam::pileup::Pileup,
    _strand_rule: StrandRule,
}

#[derive(new)]
struct PileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    start_pos: u32,
    end_pos: u32,
    focus_positions: &'a BitVec<usize, Lsb0>,
}

impl<'a> Iterator for PileupIter<'a> {
    type Item = StrandPileup;

    fn next(&mut self) -> Option<Self::Item> {
        let mut pileup: Option<Self::Item> = None;
        while let Some(Ok(plp)) = self.pileups.next() {
            let off_end = plp.pos() >= self.end_pos;
            if off_end {
                // we're done
                return None;
            } else if plp.pos() < self.start_pos {
                // advance into region we're looking at
                continue;
            } else {
                let pos = plp.pos();
                let adj_pos = pos
                    .checked_sub(self.start_pos)
                    .expect("genomic position should be greater than start");
                let st = (adj_pos * 2) as usize;
                let end = st.saturating_add(2) as usize;
                let bs = &self.focus_positions[st..end];
                match (bs[0], bs[1]) {
                    (true, _) => {
                        pileup =
                            Some(StrandPileup::new(plp, StrandRule::Positive));
                        break;
                    }
                    (_, true) => {
                        pileup =
                            Some(StrandPileup::new(plp, StrandRule::Negative));
                        break;
                    }
                    _ => {
                        continue;
                    }
                }
            }
        }
        pileup
    }
}

#[derive(Debug, Hash, Eq, PartialEq, Copy, Clone, Ord, PartialOrd)]
pub enum PartitionKey {
    NoKey,
    Key(usize),
}

fn get_forward_read_base(
    alignment: &bam::pileup::Alignment,
    record: &bam::Record,
) -> Option<DnaBase> {
    alignment.qpos().and_then(|pos| {
        if pos >= record.seq_len() {
            debug!("Record position is not included in sequence?");
            None
        } else {
            DnaBase::parse(record.seq()[pos] as char).ok()
        }
    })
}

#[derive(Debug)]
pub struct ModBasePileup2 {
    pub chrom_name: String,
    pub(crate) position_feature_counts: Vec<PileupFeatureCounts2>,
    pub(crate) interval_width: usize,
    pub(crate) stride: usize,
    pub(crate) failed_records: usize,
    pub(crate) reference_equal_records: usize,
    pub(crate) phased_feature_counts: [Vec<PileupFeatureCounts2>; 2],
}

impl ModBasePileup2 {
    fn new_empty() -> Self {
        Self {
            chrom_name: "EMPTY".to_string(),
            position_feature_counts: Vec::new(),
            interval_width: 0,
            stride: 0,
            failed_records: 0,
            reference_equal_records: 0,
            phased_feature_counts: [Vec::new(), Vec::new()],
        }
    }
}

#[derive(Clone, Debug)]
pub enum PileupNumericOptions {
    Passthrough,
    Combine,
    Collapse(CollapseMethod),
}

impl PileupNumericOptions {
    fn get_collapse_method(&self) -> Option<&CollapseMethod> {
        match self {
            Self::Collapse(method) => Some(method),
            _ => None,
        }
    }
}
