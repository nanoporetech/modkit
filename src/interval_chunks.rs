use std::collections::{BTreeMap, HashMap, VecDeque};

use anyhow::{anyhow, bail};
use derive_new::new;
use log::debug;
use rustc_hash::FxHashMap;

use crate::find_motifs::motif_bed::{
    MotifInfo, MotifLocations, MultipleMotifLocations, RegexMotif,
};
use crate::position_filter::{GenomeIntervals, Iv, StrandedPositionFilter};
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
        pos_lapper: GenomeIntervals<()>,
        neg_lapper: GenomeIntervals<()>,
    },
    AllPositions,
}

impl FocusPositions {
    fn new_motif(
        multiple_motif_locations: &MultipleMotifLocations,
        chrom_tid: u32,
        start: u32,
        end: u32,
    ) -> Self {
        let mut positions = FxHashMap::<u32, StrandRule>::default();
        let mut positive_motif_ids = FxHashMap::<u32, Vec<usize>>::default();
        let mut negative_motif_ids = FxHashMap::<u32, Vec<usize>>::default();

        let all_single_base = multiple_motif_locations
            .motif_locations
            .iter()
            .all(|ml| ml.motif_length() == 1);
        if multiple_motif_locations.motif_locations.len() == 1
            && all_single_base
        {
            // simple case, we just have 1 single base to worry about;
            let motif = &multiple_motif_locations.motif_locations[0];
            let positions_this_contig =
                motif.get_locations_unchecked(chrom_tid);
            for (position, strand_rule) in
                positions_this_contig.range(start..end)
            {
                positions.insert(*position, *strand_rule);
                match strand_rule {
                    StrandRule::Positive => {
                        positive_motif_ids.insert(*position, vec![0]);
                    }
                    StrandRule::Negative => {
                        negative_motif_ids.insert(*position, vec![0]);
                    }
                    StrandRule::Both => {
                        positive_motif_ids.insert(*position, vec![0]);
                        negative_motif_ids.insert(*position, vec![0]);
                    }
                }
            }
        } else if multiple_motif_locations.motif_locations.len() == 1 {
            // we need to potentially combine the positions together, e.g.
            // CGCG-2 and CG-0
            let motif = &multiple_motif_locations.motif_locations[0];
            let positions_this_contig =
                motif.get_locations_unchecked(chrom_tid);
            for (position, strand_rule) in
                positions_this_contig.range(start..end)
            {
                if let Some(rule) = positions.get_mut(position) {
                    *rule = rule.combine(*strand_rule);
                } else {
                    positions.insert(*position, *strand_rule);
                }
                match strand_rule {
                    StrandRule::Positive => {
                        positive_motif_ids.insert(*position, vec![0]);
                    }
                    StrandRule::Negative => {
                        negative_motif_ids.insert(*position, vec![0]);
                    }
                    StrandRule::Both => {
                        positive_motif_ids.insert(*position, vec![0]);
                        negative_motif_ids.insert(*position, vec![0]);
                    }
                }
            }
        } else if all_single_base {
            // we have more than 1 single base
            let single_base_motifs = multiple_motif_locations
                .motif_locations
                .iter()
                .enumerate()
                .map(|(id, ml)| (ml.motif().raw_motif.as_str(), (id, ml)))
                .collect::<HashMap<&str, (usize, &MotifLocations)>>();
            Self::add_single_base_motifs(
                &mut positions,
                &mut positive_motif_ids,
                &mut negative_motif_ids,
                "A",
                "T",
                &single_base_motifs,
                chrom_tid,
                start,
                end,
            );
            Self::add_single_base_motifs(
                &mut positions,
                &mut positive_motif_ids,
                &mut negative_motif_ids,
                "C",
                "G",
                &single_base_motifs,
                chrom_tid,
                start,
                end,
            );
        } else {
            // we have some mixture of the above cases.. this is the most
            // expensive.
            for (motif_id, motif) in
                multiple_motif_locations.motif_locations.iter().enumerate()
            {
                let positions_this_contig =
                    motif.get_locations_unchecked(chrom_tid);
                for (position, strand_rule) in
                    positions_this_contig.range(start..end)
                {
                    if let Some(rule) = positions.get_mut(position) {
                        *rule = rule.combine(*strand_rule);
                    } else {
                        positions.insert(*position, *strand_rule);
                    }
                    match strand_rule {
                        StrandRule::Positive => {
                            positive_motif_ids
                                .entry(*position)
                                .or_insert(Vec::new())
                                .push(motif_id);
                        }
                        StrandRule::Negative => {
                            negative_motif_ids
                                .entry(*position)
                                .or_insert(Vec::new())
                                .push(motif_id);
                        }
                        StrandRule::Both => {
                            positive_motif_ids
                                .entry(*position)
                                .or_insert(Vec::new())
                                .push(motif_id);
                            negative_motif_ids
                                .entry(*position)
                                .or_insert(Vec::new())
                                .push(motif_id);
                        }
                    }
                }
            }
        }

        Self::Motif { positions, positive_motif_ids, negative_motif_ids }
    }

    fn add_single_base_motifs(
        positions: &mut FxHashMap<u32, StrandRule>,
        positive_motif_ids: &mut FxHashMap<u32, Vec<usize>>,
        negative_motif_ids: &mut FxHashMap<u32, Vec<usize>>,
        top_base: &str,
        bottom_base: &str,
        single_base_motifs: &HashMap<&str, (usize, &MotifLocations)>,
        chrom_tid: u32,
        start: u32,
        end: u32,
    ) {
        // if we have A..
        if let Some((a_id, a_mls)) = single_base_motifs.get(top_base) {
            let positions_this_contig =
                a_mls.get_locations_unchecked(chrom_tid);
            // check if we have T, if so, we only need the A positions and
            // make the strand rule "Both"
            if let Some((t_id, _t_mls)) = single_base_motifs.get(bottom_base) {
                for (position, _strand_rule) in
                    positions_this_contig.range(start..end)
                {
                    positions.insert(*position, StrandRule::Both);
                    positive_motif_ids.insert(*position, vec![*a_id, *t_id]);
                    negative_motif_ids.insert(*position, vec![*a_id, *t_id]);
                }
            } else {
                // else, just add the A positions, and we know the strand
                // rule doesn't need to be combined
                for (position, strand_rule) in
                    positions_this_contig.range(start..end)
                {
                    positions.insert(*position, *strand_rule);
                    match strand_rule {
                        StrandRule::Positive => {
                            positive_motif_ids.insert(*position, vec![*a_id]);
                        }
                        StrandRule::Negative => {
                            negative_motif_ids.insert(*position, vec![*a_id]);
                        }
                        StrandRule::Both => {} // can't happen
                    }
                }
            }
        }
    }

    fn new_motif_combine_strands(
        motif_positions: &MultipleMotifLocations,
        chrom_tid: u32,
        start: u32,
        end: u32,
    ) -> Self {
        let mut positions = FxHashMap::<u32, StrandRule>::default();
        let mut positive_motifs = BTreeMap::new();
        let mut negative_motif_ids = FxHashMap::default();
        for (motif_id, motif) in
            motif_positions.motif_locations.iter().enumerate()
        {
            let positions_this_contig =
                motif.get_locations_unchecked(chrom_tid);
            for (position, strand_rule) in
                positions_this_contig.range(start..end)
            {
                if let Some(rule) = positions.get_mut(position) {
                    *rule = rule.combine(*strand_rule);
                } else {
                    positions.insert(*position, *strand_rule);
                }
                match strand_rule {
                    // todo make sure I want both here.. probably doesn't matter
                    //  since a motif can't really be both.
                    StrandRule::Positive | StrandRule::Both => {
                        let motif_info = motif.motif().motif_info;
                        positive_motifs
                            .entry(*position)
                            .or_insert(Vec::new())
                            .push((motif_info, motif_id));
                    }
                    StrandRule::Negative => {
                        negative_motif_ids
                            .entry(*position)
                            .or_insert(Vec::new())
                            .push(motif_id);
                    }
                }
            }
        }

        Self::MotifCombineStrands {
            positions,
            positive_motifs,
            negative_motif_ids,
        }
    }

    fn new_regions(
        stranded_position_filter: &StrandedPositionFilter<()>,
        chrom_id: u32,
        start: u32,
        end: u32,
    ) -> Self {
        let start = start as u64;
        let stop = end as u64;

        let pos_intervals = stranded_position_filter
            .pos_positions
            .get(&chrom_id)
            .map(|lp| {
                lp.find(start, stop)
                    .into_iter()
                    .map(|iv| Iv {
                        start: std::cmp::max(iv.start, start),
                        stop: std::cmp::min(iv.stop, stop),
                        val: (),
                    })
                    .collect::<Vec<Iv>>()
            })
            .unwrap_or(Vec::new());

        let neg_intervals = stranded_position_filter
            .neg_positions
            .get(&chrom_id)
            .map(|lp| {
                lp.find(start, stop)
                    .into_iter()
                    .map(|iv| Iv {
                        start: std::cmp::max(iv.start, start),
                        stop: std::cmp::min(iv.stop, stop),
                        val: (),
                    })
                    .collect::<Vec<Iv>>()
            })
            .unwrap_or(Vec::new());
        let mut pos_lapper = GenomeIntervals::new(pos_intervals);
        pos_lapper.merge_overlaps();
        let mut neg_lapper = GenomeIntervals::new(neg_intervals);
        neg_lapper.merge_overlaps();

        Self::Regions { pos_lapper, neg_lapper }
    }

    // semantics: return Some iff we keep that position, else None
    pub fn check_position(&self, pos: &u32) -> Option<StrandRule> {
        match &self {
            FocusPositions::Motif { positions, .. }
            | FocusPositions::MotifCombineStrands { positions, .. } => {
                positions.get(pos).map(|sr| *sr)
            }
            FocusPositions::Regions { pos_lapper, neg_lapper } => {
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

pub struct ChromCoordinates {
    pub chrom_tid: u32,
    pub start_pos: u32,
    pub end_pos: u32,
    pub focus_positions: FocusPositions,
}

impl ChromCoordinates {
    fn new(
        chrom_tid: u32,
        start_pos: u32,
        end_pos: u32,
        combine_strands: bool,
        motif_positions: Option<&MultipleMotifLocations>,
        position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> Self {
        // todo/warn currently the assumption is made that motifs, if given,
        // have  been pre-filtered so that the position filter can be
        // ignored..
        let focus_positions = match (motif_positions, position_filter) {
            (Some(motif), _) => {
                if combine_strands {
                    FocusPositions::new_motif_combine_strands(
                        motif, chrom_tid, start_pos, end_pos,
                    )
                } else {
                    FocusPositions::new_motif(
                        motif, chrom_tid, start_pos, end_pos,
                    )
                }
            }
            (_, Some(spf)) => {
                FocusPositions::new_regions(spf, chrom_tid, start_pos, end_pos)
            }
            (None, None) => FocusPositions::AllPositions,
        };

        Self { chrom_tid, start_pos, end_pos, focus_positions }
    }

    pub(crate) fn len(&self) -> u32 {
        self.end_pos.checked_sub(self.start_pos).unwrap_or(0u32)
    }

    pub(crate) fn merge(self, other: Self) -> Self {
        match (&self.focus_positions, &other.focus_positions) {
            (FocusPositions::AllPositions, FocusPositions::AllPositions) => {}
            _ => todo!("must be 'AllPositions' to merge"),
        }
        assert_eq!(self.chrom_tid, other.chrom_tid);
        Self {
            chrom_tid: self.chrom_tid,
            start_pos: std::cmp::min(self.start_pos, other.start_pos),
            end_pos: std::cmp::max(self.end_pos, other.end_pos),
            focus_positions: FocusPositions::AllPositions,
        }
    }
}

#[derive(new)]
pub struct MultiChromCoordinates(pub Vec<ChromCoordinates>);

impl MultiChromCoordinates {
    pub fn total_length(&self) -> u64 {
        self.0.iter().map(|cc| cc.len() as u64).sum::<u64>()
    }
}

pub struct ReferenceIntervalsFeeder {
    contigs: VecDeque<ReferenceRecord>,
    batch_size: usize,
    interval_size: u32,
    motifs: Option<MultipleMotifLocations>,
    position_filter: Option<StrandedPositionFilter<()>>,
    combine_strands: bool,
    curr_contig: ReferenceRecord,
    curr_position: u32,
    done: bool,
}

impl ReferenceIntervalsFeeder {
    pub fn new(
        reference_records: Vec<ReferenceRecord>,
        batch_size: usize,
        interval_size: u32,
        combine_strands: bool,
        multi_motif_locations: Option<MultipleMotifLocations>,
        position_filter: Option<StrandedPositionFilter<()>>,
    ) -> anyhow::Result<Self> {
        if combine_strands & !multi_motif_locations.is_some() {
            bail!("cannot combine strands without a motif")
        }
        let mut contigs = if let Some(position_filter) =
            position_filter.as_ref()
        {
            reference_records
                .into_iter()
                .filter_map(|contig| {
                    position_filter.contig_ends(&contig.tid).map(|(s, t)| {
                        let length = t.checked_sub(s).unwrap_or(0u64);
                        debug!(
                            "narrowing record {} to {s}-{t} ({length} bases)",
                            contig.name.as_str()
                        );
                        ReferenceRecord::new(
                            contig.tid,
                            s as u32,
                            length as u32,
                            contig.name,
                        )
                    })
                })
                .collect::<VecDeque<_>>()
        } else {
            reference_records.into_iter().collect::<VecDeque<_>>()
        };

        if contigs.len() == 1 {
            debug!("there is a single contig to work on");
        } else {
            debug!("there are {} contigs to work on", contigs.len());
        }
        let curr_contig = contigs
            .pop_front()
            .ok_or(anyhow!("should be at least 1 contig"))?;
        let curr_position = curr_contig.start;
        Ok(Self {
            contigs,
            batch_size,
            interval_size,
            motifs: multi_motif_locations,
            combine_strands,
            position_filter,
            curr_contig,
            curr_position,
            done: false,
        })
    }

    pub fn total_length(&self) -> usize {
        self.contigs
            .iter()
            .fold(self.curr_contig.length as usize, |agg, next| {
                agg.saturating_add(next.length as usize)
            })
    }

    fn find_end(&self, mut end: u32, tid: u32) -> u32 {
        if let Some(mls) = self.motifs.as_ref() {
            'creep_loop: loop {
                let motifs_at_position = mls
                    .motifs_at_position(tid, end.checked_sub(1).unwrap_or(end))
                    .into_iter()
                    .filter(|motif| motif.length() > 1)
                    .collect::<Vec<&RegexMotif>>();
                if motifs_at_position.is_empty() {
                    break 'creep_loop;
                }
                let creep_out = motifs_at_position
                    .iter()
                    .map(|motif| {
                        motif
                            .length()
                            .checked_sub(motif.forward_offset())
                            .unwrap_or(motif.length())
                    })
                    .max();
                if let Some(creep) = creep_out {
                    end += creep as u32;
                    if creep == 0 {
                        // just in case so that we don't get stuck
                        // forever
                        break 'creep_loop;
                    } else {
                        continue 'creep_loop;
                    }
                } else {
                    // shouldn't ever really hit this since we have
                    // a break when motifs_at_position is empty
                    break 'creep_loop;
                }
            }
        }
        end
    }

    fn update_current(&mut self) {
        if let Some(reference_record) = self.contigs.pop_front() {
            self.curr_position = reference_record.start;
            self.curr_contig = reference_record;
        } else {
            debug!("no more records to process");
            self.done = true;
        }
    }
}

impl Iterator for ReferenceIntervalsFeeder {
    type Item = Vec<MultiChromCoordinates>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut ret = Vec::new();

        let mut batch = Vec::new();
        let mut batch_length = 0u32;

        loop {
            if self.done {
                // debug!("done!");
                break;
            } else if ret.len() >= self.batch_size {
                break;
            }
            debug_assert!(self.curr_position < self.curr_contig.end());
            let start = self.curr_position;
            let tid = self.curr_contig.tid;
            // in the case where we're on a large chrom end will be < length,
            // but batch length will be equal to interval size
            let end = std::cmp::min(
                start + self.interval_size,
                self.curr_contig.end(),
            );
            // if self.combine_strands...
            let end = self.find_end(end, tid);
            let end = std::cmp::min(end, self.curr_contig.end());
            // in the "short contig" case, chrom_coords.len() will be less than
            // interval size so batch length will be less than
            // interval size for a few rounds
            let chrom_coords = ChromCoordinates::new(
                tid,
                start,
                end,
                self.combine_strands,
                self.motifs.as_ref(),
                self.position_filter.as_ref(),
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
            if end >= self.curr_contig.end() {
                self.update_current();
            } else {
                self.curr_position = end;
            }
        }

        if !batch.is_empty() {
            ret.push(MultiChromCoordinates::new(batch));
        }

        if ret.is_empty() {
            None
        } else {
            Some(ret)
        }
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

    use crate::interval_chunks::slice_dna_sequence;
    use crate::test_utils::load_test_sequence;

    #[test]
    fn test_check_sequence_slicing_is_same_as_fetch() {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
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
