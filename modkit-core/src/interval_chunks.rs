use std::collections::{BTreeMap, BTreeSet, VecDeque};

use anyhow::{anyhow, bail};
use bitvec::bitvec;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use derive_new::new;
use itertools::Itertools;
use log::debug;
use rustc_hash::FxHashMap;

use crate::errs::MkResult;
use crate::fasta::MotifLocationsLookup;
use crate::mod_base_code::DnaBase;
use crate::motifs::motif_bed::MotifInfo;
use crate::position_filter::{GenomeIntervals, StrandedPositionFilter};
use crate::util::{ReferenceRecord, StrandRule};

pub fn slice_dna_sequence(str_seq: &str, start: usize, end: usize) -> String {
    str_seq
        .char_indices()
        .filter_map(
            |(pos, nt)| {
                if pos >= start && pos <= end {
                    Some(nt)
                } else {
                    None
                }
            },
        )
        .collect::<String>()
}

pub(crate) enum FocusPositions2 {
    MotifMask { mask: BitVec<usize, Lsb0>, num_motifs: u8 },
    MaskedPositions { mask: BitVec<usize, Lsb0> },
    AllPositions,
}

/// A "kitchen-sink" enum for different situations (mostly in pileup).
pub enum FocusPositions {
    Motif {
        // positions and rules for all motifs "stacked up"
        positions: FxHashMap<u32, StrandRule>,
        // mapping of position to IDs of motifs at that position
        positive_motif_ids: FxHashMap<u32, Vec<usize>>,
        // mapping of position to the motif ids on the negative stand
        negative_motif_ids: FxHashMap<u32, Vec<usize>>,
    },
    MotifCombineStrands {
        // positions and rules for all motifs "stacked up"
        positions: FxHashMap<u32, StrandRule>,
        // positions to a list of motifs present at that position and their ID,
        // only need positive because, at present, we just care about this data
        // for combining counts together.
        positive_motifs: BTreeMap<u32, Vec<(MotifInfo, usize)>>,
        // n.b. when combining strands, we don't need to know about
        // negative strand positions, but we need to know which motif
        // IDs are at each site when combining them together.
        negative_motif_ids: FxHashMap<u32, Vec<usize>>,
    },
    Regions {
        // positions from an extracted BED file
        pos_intervals: GenomeIntervals<()>,
        neg_intervals: GenomeIntervals<()>,
    },
    AllPositions,
}

impl FocusPositions {
    // semantics: return Some iff we keep that position, else None
    pub fn check_position(&self, pos: &u32) -> Option<StrandRule> {
        match &self {
            FocusPositions::Motif { positions, .. }
            | FocusPositions::MotifCombineStrands { positions, .. } => {
                positions.get(pos).map(|sr| *sr)
            }
            FocusPositions::Regions {
                pos_intervals: pos_lapper,
                neg_intervals: neg_lapper,
            } => {
                let pos = *pos as u64;
                let pos_hit = pos_lapper.find(pos, pos + 1).count() > 0;
                let neg_hit = neg_lapper.find(pos, pos + 1).count() > 0;
                match (pos_hit, neg_hit) {
                    (true, true) => Some(StrandRule::Both),
                    (true, false) => Some(StrandRule::Positive),
                    (false, true) => Some(StrandRule::Negative),
                    (false, false) => None,
                }
            }
            FocusPositions::AllPositions => Some(StrandRule::Both),
        }
    }

    pub fn get_positive_strand_motif_ids(
        &self,
        pos: &u32,
    ) -> Option<Vec<usize>> {
        match &self {
            FocusPositions::Motif { positive_motif_ids, .. } => {
                positive_motif_ids
                    .get(pos)
                    .map(|ids| ids.iter().copied().collect())
            }
            FocusPositions::MotifCombineStrands { positive_motifs, .. } => {
                positive_motifs.get(&pos).map(|motifs_at_position| {
                    motifs_at_position.iter().map(|(_, id)| *id).collect()
                })
            }
            _ => None,
        }
    }

    pub fn get_negative_strand_motif_ids(
        &self,
        pos: &u32,
    ) -> Option<Vec<usize>> {
        match &self {
            FocusPositions::Motif { negative_motif_ids, .. }
            | FocusPositions::MotifCombineStrands {
                negative_motif_ids, ..
            } => negative_motif_ids
                .get(pos)
                .map(|ids| ids.iter().copied().collect()),
            _ => None,
        }
    }
}

#[derive(new)]
pub(crate) struct ChromCoordinates {
    pub chrom_tid: u32,
    pub start_pos: u32,
    pub end_pos: u32,
    pub focus_positions: FocusPositions2,
    pub final_interval: bool,
}

impl ChromCoordinates {
    pub(crate) fn len(&self) -> u32 {
        self.end_pos.checked_sub(self.start_pos).unwrap_or(0u32)
    }

    pub(crate) fn merge(self, other: Self) -> Self {
        match (&self.focus_positions, &other.focus_positions) {
            (FocusPositions2::AllPositions, FocusPositions2::AllPositions) => {}
            _ => todo!("must be 'AllPositions' to merge"),
        }
        assert_eq!(self.chrom_tid, other.chrom_tid);
        Self {
            chrom_tid: self.chrom_tid,
            start_pos: std::cmp::min(self.start_pos, other.start_pos),
            end_pos: std::cmp::max(self.end_pos, other.end_pos),
            focus_positions: FocusPositions2::AllPositions,
            final_interval: self.final_interval || other.final_interval,
        }
    }
}

#[derive(new)]
pub(crate) struct MultiChromCoordinates(pub Vec<ChromCoordinates>);

impl MultiChromCoordinates {
    pub fn total_length(&self) -> u64 {
        self.0.iter().map(|cc| cc.len() as u64).sum::<u64>()
    }
}

pub(crate) trait TotalLength {
    fn total_length(&self) -> u64;
}

impl TotalLength for Vec<MultiChromCoordinates> {
    fn total_length(&self) -> u64 {
        self.iter().map(|x| x.total_length()).sum::<u64>()
    }
}

impl TotalLength for ReferenceIntervalBatchesFeeder {
    fn total_length(&self) -> u64 {
        self.contigs.iter().fold(
            self.curr_contig
                .as_ref()
                .map(|record| record.length as u64)
                .unwrap_or(0),
            |agg, next| agg.saturating_add(next.length as u64),
        )
    }
}

pub(crate) struct ReferenceIntervalBatchesFeeder {
    contigs: VecDeque<ReferenceRecord>,
    batch_size: usize,
    interval_size: u32,
    motifs: Option<MotifLocationsLookup>,
    position_filter: Option<StrandedPositionFilter<()>>,
    combine_strands: bool,
    // None denotes exhaustion, including an explicit zero-target terminal.
    curr_contig: Option<ReferenceRecord>,
    curr_position: u32,
}

impl ReferenceIntervalBatchesFeeder {
    pub fn new(
        reference_records: Vec<ReferenceRecord>,
        batch_size: usize,
        interval_size: u32,
        combine_strands: bool,
        multi_motif_locations: Option<MotifLocationsLookup>,
        position_filter: Option<StrandedPositionFilter<()>>,
    ) -> anyhow::Result<Self> {
        Self::new_inner(
            reference_records,
            batch_size,
            interval_size,
            combine_strands,
            multi_motif_locations,
            position_filter,
            false,
        )
    }

    /// Allow a zero-reference input to produce an already exhausted feeder.
    /// This is only appropriate when the caller separately processes unmapped
    /// records after the mapped interval phase.
    pub(crate) fn new_allowing_zero_reference_targets(
        reference_records: Vec<ReferenceRecord>,
        batch_size: usize,
        interval_size: u32,
        combine_strands: bool,
        multi_motif_locations: Option<MotifLocationsLookup>,
        position_filter: Option<StrandedPositionFilter<()>>,
    ) -> anyhow::Result<Self> {
        Self::new_inner(
            reference_records,
            batch_size,
            interval_size,
            combine_strands,
            multi_motif_locations,
            position_filter,
            true,
        )
    }

    fn new_inner(
        reference_records: Vec<ReferenceRecord>,
        batch_size: usize,
        interval_size: u32,
        combine_strands: bool,
        multi_motif_locations: Option<MotifLocationsLookup>,
        position_filter: Option<StrandedPositionFilter<()>>,
        allow_zero_reference_targets: bool,
    ) -> anyhow::Result<Self> {
        if combine_strands & !multi_motif_locations.is_some() {
            bail!("cannot combine strands without a motif")
        }
        let zero_reference_targets = reference_records.is_empty();
        let mut contigs =
            reference_records.into_iter().collect::<VecDeque<_>>();
        let n_contigs = contigs.iter().map(|r| r.tid).unique().count();

        if n_contigs == 1 {
            debug!(
                "there is a single contig to work on (in {} parts)",
                contigs.len()
            );
        } else {
            debug!(
                "there are {n_contigs} contig(s) to work on ({} parts)",
                contigs.len()
            );
        }
        let curr_contig = contigs.pop_front();
        if curr_contig.is_none()
            && !(allow_zero_reference_targets && zero_reference_targets)
        {
            return Err(anyhow!("should be at least 1 contig"));
        }
        let curr_position =
            curr_contig.as_ref().map(|record| record.start).unwrap_or(0);
        Ok(Self {
            contigs,
            batch_size,
            interval_size,
            motifs: multi_motif_locations,
            combine_strands,
            position_filter,
            curr_contig,
            curr_position,
        })
    }

    fn update_current(&mut self) {
        if let Some(reference_record) = self.contigs.pop_front() {
            self.curr_position = reference_record.start;
            self.curr_contig = Some(reference_record);
        } else {
            debug!("no more records to process");
            self.curr_contig = None;
        }
    }

    fn next_batch(
        &mut self,
    ) -> anyhow::Result<Option<Vec<MultiChromCoordinates>>> {
        let mut ret = Vec::new();

        let mut batch = Vec::new();
        let mut batch_length = 0u32;

        loop {
            if self.curr_contig.is_none() {
                break;
            } else if ret.len() >= self.batch_size {
                break;
            }
            let curr_contig = self.curr_contig.as_ref().expect(
                "non-terminal feeder must have a current reference record",
            );
            debug_assert!(self.curr_position < curr_contig.end());
            let start = self.curr_position;
            let tid = curr_contig.tid;
            // in the case where we're on a large chrom end will be < length,
            // but batch length will be equal to interval size
            let end =
                std::cmp::min(start + self.interval_size, curr_contig.end());
            // get the sequence here.
            let (focus_positions, end) =
                if let Some(lookup) = self.motifs.as_mut() {
                    // todo change everything to u64
                    let range = (start as u64)..(end as u64);
                    let (fps, end) = lookup.get_motif_positions(
                        &curr_contig.name,
                        tid,
                        curr_contig.end(),
                        range,
                        self.position_filter.as_ref(),
                        self.combine_strands,
                    )?;
                    (fps, end)
                } else if let Some(_pos_filt) = self.position_filter.as_ref() {
                    // TODO: Currently blocked by check that include-bed always
                    // has either a motif or modified bases,
                    // remove this eventually.
                    todo!()
                } else {
                    (FocusPositions2::AllPositions, end)
                };
            let end = std::cmp::min(end, curr_contig.end());
            // in the "short contig" case, chrom_coords.len() will be less than
            // interval size so batch length will be less than
            // interval size for a few rounds
            let chrom_coords = ChromCoordinates::new(
                tid,
                start,
                end,
                focus_positions,
                end == curr_contig.end(),
            );
            batch_length += chrom_coords.len();
            batch.push(chrom_coords);
            if batch_length >= self.interval_size {
                // in the "normal" case, the batch will have 1 element
                let finished_batch = std::mem::replace(&mut batch, Vec::new());
                ret.push(MultiChromCoordinates::new(finished_batch));
                batch_length = 0;
            } else {
                // we're going to accumulate another chrom_coords in this batch
            }
            // might need to update the pointers, check if we're at the end of
            // this contig
            if end >= curr_contig.end() {
                self.update_current();
            } else {
                self.curr_position = end;
            }
        }

        if !batch.is_empty() {
            ret.push(MultiChromCoordinates::new(batch));
        }

        if ret.is_empty() {
            Ok(None)
        } else {
            Ok(Some(ret))
        }
    }
}

impl Iterator for ReferenceIntervalBatchesFeeder {
    type Item = anyhow::Result<Vec<MultiChromCoordinates>>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_batch().transpose()
    }
}

impl TotalLength for ChromCoordinatesFeeder {
    fn total_length(&self) -> u64 {
        self.contigs.iter().fold(
            self.curr_contig
                .as_ref()
                .map(|record| record.length as u64)
                .unwrap_or(0),
            |agg, next| agg.saturating_add(next.length as u64),
        )
    }
}

#[derive(Clone)]
pub(crate) struct ChromCoordinatesFeeder {
    contigs: VecDeque<ReferenceRecord>,
    interval_size: u32,
    motifs: Option<MotifLocationsLookup>,
    // TODO: make this a borrow
    position_filter: Option<StrandedPositionFilter<()>>,
    combine_strands: bool,
    // None denotes exhaustion, including an explicit zero-target terminal.
    curr_contig: Option<ReferenceRecord>,
    curr_position: u32,
    zero_reference_target_terminal: bool,
}

impl ChromCoordinatesFeeder {
    pub(crate) fn new(
        reference_records: &[ReferenceRecord],
        interval_size: u32,
        motifs: Option<MotifLocationsLookup>,
        combine_strands: bool,
        position_filter: Option<StrandedPositionFilter<()>>,
    ) -> anyhow::Result<Self> {
        Self::new_inner(
            reference_records,
            interval_size,
            motifs,
            combine_strands,
            position_filter,
            false,
        )
    }

    /// Allow a zero-reference input to produce an already exhausted feeder.
    /// Inputs whose targets are removed by motif filtering remain errors.
    pub(crate) fn new_allowing_zero_reference_targets(
        reference_records: &[ReferenceRecord],
        interval_size: u32,
        motifs: Option<MotifLocationsLookup>,
        combine_strands: bool,
        position_filter: Option<StrandedPositionFilter<()>>,
    ) -> anyhow::Result<Self> {
        Self::new_inner(
            reference_records,
            interval_size,
            motifs,
            combine_strands,
            position_filter,
            true,
        )
    }

    fn new_inner(
        reference_records: &[ReferenceRecord],
        interval_size: u32,
        motifs: Option<MotifLocationsLookup>,
        combine_strands: bool,
        position_filter: Option<StrandedPositionFilter<()>>,
        allow_zero_reference_targets: bool,
    ) -> anyhow::Result<Self> {
        debug!("there is {} reference record(s)", reference_records.len());
        if combine_strands & !motifs.is_some() {
            bail!("cannot combine strands without a motif")
        }
        let zero_reference_targets = reference_records.is_empty();
        let mut contigs = reference_records
            .into_iter()
            .filter(|rr| {
                if let Some(motifs) = motifs.as_ref() {
                    motifs.has_sequence(&rr.name)
                } else {
                    true
                }
            })
            .cloned()
            .collect::<VecDeque<_>>();
        let n_contigs = contigs.iter().map(|r| r.tid).unique().count();

        if n_contigs == 1 {
            debug!(
                "there is a single contig to work on (in {} parts)",
                contigs.len()
            );
        } else {
            debug!(
                "there are {n_contigs} contig(s) to work on ({} parts)",
                contigs.len()
            );
        }
        let curr_contig = contigs.pop_front();
        if curr_contig.is_none()
            && !(allow_zero_reference_targets && zero_reference_targets)
        {
            return Err(anyhow!(
                "should be at least 1 contig, are all of the reads unmapped?"
            ));
        }
        let curr_position =
            curr_contig.as_ref().map(|record| record.start).unwrap_or(0);
        let zero_reference_target_terminal =
            allow_zero_reference_targets && zero_reference_targets;
        Ok(Self {
            contigs,
            interval_size,
            motifs,
            position_filter,
            combine_strands,
            curr_contig,
            curr_position,
            zero_reference_target_terminal,
        })
    }

    fn update_current(&mut self, end: u32) {
        let curr_contig_end = self
            .curr_contig
            .as_ref()
            .expect("non-terminal feeder must have a current reference record")
            .end();
        if end >= curr_contig_end {
            if let Some(rr) = self.contigs.pop_front() {
                self.curr_position = rr.start;
                self.curr_contig = Some(rr);
            } else {
                debug!("feeder is done");
                self.curr_contig = None;
            }
        } else {
            self.curr_position = end
        }
    }

    fn get_next(&mut self) -> MkResult<Option<ChromCoordinates>> {
        let Some(curr_contig) = self.curr_contig.as_ref() else {
            return Ok(None);
        };
        let start = self.curr_position;
        let tid = curr_contig.tid;
        let end = std::cmp::min(
            start.saturating_add(self.interval_size),
            curr_contig.end(),
        );

        let (focus_positions, end) = if let Some(lookup) = self.motifs.as_mut()
        {
            // todo change everything to u64
            let range = (start as u64)..(end as u64);
            let (fps, end) = lookup.get_motif_positions(
                &curr_contig.name,
                tid,
                curr_contig.end(),
                range,
                self.position_filter.as_ref(),
                self.combine_strands,
            )?;
            (fps, end)
        } else if let Some(pos_filt) = self.position_filter.as_ref() {
            let mut mask = bitvec![0; (end.saturating_sub(start) as usize) * 2];
            if !pos_filt.contains_chrom_id(&(tid as i64)) {
                (FocusPositions2::MaskedPositions { mask }, end)
            } else {
                if let Some(pos_strand_intervals) =
                    pos_filt.pos_positions.get(&tid)
                {
                    let mut cursor = 0usize;
                    for (i, pos) in (start..end).enumerate() {
                        if pos_strand_intervals
                            .seek(
                                pos as u64,
                                pos.saturating_add(1) as u64,
                                &mut cursor,
                            )
                            .count()
                            >= 1
                        {
                            mask.set(i * 2, true);
                        }
                    }
                }
                if let Some(negative_strand_intervals) =
                    pos_filt.neg_positions.get(&tid)
                {
                    let mut cursor = 0usize;
                    for (i, pos) in (start..end).enumerate() {
                        if negative_strand_intervals
                            .seek(
                                pos as u64,
                                pos.saturating_add(1) as u64,
                                &mut cursor,
                            )
                            .count()
                            >= 1
                        {
                            mask.set((i * 2) + 1, true);
                        }
                    }
                }
                (FocusPositions2::MaskedPositions { mask }, end)
            }
        } else {
            (FocusPositions2::AllPositions, end)
        };

        let end = std::cmp::min(end, curr_contig.end());
        let chrom_coords = ChromCoordinates::new(
            tid,
            start,
            end,
            focus_positions,
            end == curr_contig.end(),
        );
        self.update_current(end);

        Ok(Some(chrom_coords))
    }

    pub(crate) fn ref_motifs(&self) -> Option<&MotifLocationsLookup> {
        self.motifs.as_ref()
    }

    pub(crate) fn has_position_filter(&self) -> bool {
        self.position_filter.is_some()
    }

    pub(crate) fn is_zero_reference_target_terminal(&self) -> bool {
        self.zero_reference_target_terminal
    }

    pub(crate) fn get_motif_bases(&self) -> Option<[DnaBase; 4]> {
        if let Some(motif_lookup) = self.ref_motifs() {
            let motif_primary_bases = motif_lookup
                .motif_infos()
                .map(|motif_info| motif_info.primary_base)
                .collect::<BTreeSet<DnaBase>>();
            let mut motif_bases = [DnaBase::A; 4];
            for (i, b) in motif_primary_bases.into_iter().enumerate() {
                motif_bases[i] = b
            }
            Some(motif_bases)
        } else {
            None
        }
    }
}

impl Iterator for ChromCoordinatesFeeder {
    type Item = MkResult<ChromCoordinates>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next().transpose()
    }
}

#[derive(Copy, Clone)]
struct ChromEnd {
    chrom_id: u32,
    end_pos: u32,
}

#[derive(new)]
pub(crate) struct LinkedChromCoordinates {
    chrom_coordinates: ChromCoordinates,
    prev_chrom_end: Option<ChromEnd>,
}

impl LinkedChromCoordinates {
    pub(crate) fn chrom_tid(&self) -> u32 {
        self.chrom_coordinates.chrom_tid
    }

    pub(crate) fn start_pos(&self) -> u32 {
        self.chrom_coordinates.start_pos
    }

    pub(crate) fn end_pos(&self) -> u32 {
        self.chrom_coordinates.end_pos
    }

    pub(crate) fn prev_end(&self) -> Option<u32> {
        self.prev_chrom_end.map(|x| x.end_pos)
    }
}

pub(crate) struct MultiLinkedChromCoordinates(pub Vec<LinkedChromCoordinates>);

impl MultiLinkedChromCoordinates {
    pub(crate) fn total_length(&self) -> u64 {
        self.0.iter().map(|c| c.chrom_coordinates.len() as u64).sum::<u64>()
    }
}

pub(crate) struct LastEnd<
    I: Iterator<Item = Vec<MultiChromCoordinates>> + Sized,
> {
    iter: I,
    chrom_end: Option<ChromEnd>,
}

pub(crate) trait WithPrevEnd<I: Iterator<Item = Vec<MultiChromCoordinates>>> {
    fn with_prev_end(self) -> LastEnd<Self>
    where
        Self: Iterator<Item = Vec<MultiChromCoordinates>> + Sized,
    {
        LastEnd { iter: self, chrom_end: None }
    }
}

impl<I: Iterator<Item = Vec<MultiChromCoordinates>>> WithPrevEnd<I> for I {}

impl<I> Iterator for LastEnd<I>
where
    I: Iterator<Item = Vec<MultiChromCoordinates>>,
{
    type Item = Vec<MultiLinkedChromCoordinates>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(chrom_coords) = self.iter.next() {
            let (super_batch, last_batch_end) = chrom_coords.into_iter().fold(
                (Vec::new(), self.chrom_end),
                |(mut super_batch_acc, prev_batch_end), chrom_coords| {
                    let (linked_chrom_coords, batch_end) =
                        chrom_coords.0.into_iter().fold(
                            (Vec::new(), prev_batch_end),
                            |(mut acc, prev_end), next| {
                                // decide whether to keep the previous end:
                                let prev_end = prev_end.and_then(|pe| {
                                    // if this interval is on the same chrom:
                                    // keep
                                    if pe.chrom_id == next.chrom_tid {
                                        Some(pe)
                                    } else {
                                        // if we're on a new chrom, don't keep
                                        None
                                    }
                                });

                                let end = Some(ChromEnd {
                                    chrom_id: next.chrom_tid,
                                    end_pos: next.end_pos,
                                });
                                acc.push(LinkedChromCoordinates::new(
                                    next, prev_end,
                                ));
                                (acc, end)
                            },
                        );
                    super_batch_acc
                        .push(MultiLinkedChromCoordinates(linked_chrom_coords));
                    (super_batch_acc, batch_end)
                },
            );
            self.chrom_end = last_batch_end;
            Some(super_batch)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod interval_chunks_tests {
    use rust_htslib::faidx;

    use crate::interval_chunks::{
        slice_dna_sequence, ChromCoordinatesFeeder,
        ReferenceIntervalBatchesFeeder, TotalLength,
    };
    use crate::test_utils::load_test_sequence;
    use crate::util::ReferenceRecord;

    #[test]
    fn test_check_sequence_slicing_is_same_as_fetch() {
        let fasta_fp = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
        let fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let name = "oligo_1512_adapters";
        let dna = load_test_sequence(name);
        let start = 49;
        let end = 99;
        let slice_a = slice_dna_sequence(&dna, start, end);
        let slice_b = fasta_reader.fetch_seq_string(name, start, end).unwrap();
        assert_eq!(slice_a, slice_b);
    }

    #[test]
    fn strict_feeders_reject_zero_reference_targets() {
        assert!(ReferenceIntervalBatchesFeeder::new(
            Vec::new(),
            1,
            100,
            false,
            None,
            None,
        )
        .is_err());
        assert!(
            ChromCoordinatesFeeder::new(&[], 100, None, false, None).is_err()
        );
    }

    #[test]
    fn interval_batches_zero_reference_terminal_is_exhausted() {
        let mut feeder =
            ReferenceIntervalBatchesFeeder::new_allowing_zero_reference_targets(
                Vec::new(),
                1,
                100,
                false,
                None,
                None,
            )
            .unwrap();
        assert_eq!(feeder.total_length(), 0);
        assert!(feeder.next().is_none());
        assert!(feeder.next().is_none());
    }

    #[test]
    fn chrom_coordinates_zero_reference_terminal_and_clone_are_exhausted() {
        let reference_record =
            ReferenceRecord::new(0, 0, 1, "chr1".to_string());
        let active =
            ChromCoordinatesFeeder::new_allowing_zero_reference_targets(
                &[reference_record],
                100,
                None,
                false,
                None,
            )
            .unwrap();
        assert!(!active.is_zero_reference_target_terminal());

        let mut feeder =
            ChromCoordinatesFeeder::new_allowing_zero_reference_targets(
                &[],
                100,
                None,
                false,
                None,
            )
            .unwrap();
        let mut cloned = feeder.clone();
        assert!(feeder.is_zero_reference_target_terminal());
        assert!(cloned.is_zero_reference_target_terminal());
        assert_eq!(feeder.total_length(), 0);
        assert_eq!(cloned.total_length(), 0);
        assert!(feeder.next().is_none());
        assert!(cloned.next().is_none());
        assert!(feeder.next().is_none());
        assert!(cloned.next().is_none());
    }

    #[test]
    fn test_interval_chunks() {

        // let seq = "ABCDEF".chars().collect::<Vec<char>>();
        // let mut ic = IntervalChunks::new(0, seq.len() as u32, 3);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 0);
        // assert_eq!(e, 3);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['A', 'B', 'C']);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 2);
        // assert_eq!(e, 5);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['C', 'D', 'E']);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 4);
        // assert_eq!(e, 6);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['E', 'F']);
    }
}
