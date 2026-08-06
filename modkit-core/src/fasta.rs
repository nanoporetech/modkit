use std::path::PathBuf;

use anyhow::{bail, Context};
use bitvec::bitvec;
use bitvec::vec::BitVec;
use itertools::Itertools;
use log::debug;
use rust_htslib::faidx;
use rustc_hash::FxHashMap;
use substring::Substring;

use crate::errs::{MkError, MkResult};
use crate::interval_chunks::FocusPositions2;
use crate::motifs::motif_bed::{find_motif_hits, MotifInfo, RegexMotif};
use crate::position_filter::StrandedPositionFilter;
use crate::util::Strand;

#[derive(Clone)]
struct HtsFastaHandle {
    fasta_fp: PathBuf,
    contigs: FxHashMap<String, u64>,
    preloaded: bool,
    sequences: FxHashMap<String, String>,
    #[cfg(test)]
    sequence_fetches: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    #[cfg(test)]
    sequence_bases_fetched: std::sync::Arc<std::sync::atomic::AtomicU64>,
}

impl HtsFastaHandle {
    fn from_file(fp: &PathBuf, preload: bool) -> anyhow::Result<Self> {
        let tmp_reader = faidx::Reader::from_path(fp)?;
        let contigs = (0..tmp_reader.n_seqs()).try_fold(
            FxHashMap::default(),
            |mut contigs, i| {
                let seq_name = tmp_reader
                    .seq_name(i as i32)
                    .context("failed to get name of {i}th sequence")?;
                let length = tmp_reader.fetch_seq_len(&seq_name);
                if let Some(_) = contigs.insert(seq_name.clone(), length) {
                    bail!("{seq_name} is in FASTA more than once")
                } else {
                    Ok(contigs)
                }
            },
        )?;

        let sequences = if preload {
            let fasta_reader = bio::io::fasta::Reader::from_file(fp)?;
            fasta_reader
                .records()
                .map_ok(|x| {
                    (
                        x.id().to_string(),
                        x.seq().iter().map(|b| *b as char).collect::<String>(),
                    )
                })
                .collect::<Result<FxHashMap<String, String>, _>>()?
        } else {
            FxHashMap::default()
        };

        if preload {
            debug!("preloaded {} sequences", sequences.len());
        }

        Ok(Self {
            fasta_fp: fp.to_owned(),
            contigs,
            sequences,
            preloaded: preload,
            #[cfg(test)]
            sequence_fetches: Default::default(),
            #[cfg(test)]
            sequence_bases_fetched: Default::default(),
        })
    }

    fn get_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> MkResult<String> {
        if let Some(length) = self.contigs.get(contig) {
            if start > end || end > *length {
                Err(MkError::InvalidReferenceCoordinates)
            } else if start == end {
                Ok(String::new())
            } else {
                #[cfg(test)]
                {
                    self.sequence_fetches
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    self.sequence_bases_fetched.fetch_add(
                        end - start,
                        std::sync::atomic::Ordering::Relaxed,
                    );
                }
                let seq = if self.preloaded
                    && self.sequences.get(contig).is_some()
                {
                    let s = self.sequences.get(contig).unwrap();
                    s.substring(start as usize, end as usize).to_string()
                } else {
                    let tmp_reader = faidx::Reader::from_path(&self.fasta_fp)
                        .map_err(|e| MkError::HtsLibError(e))?;
                    tmp_reader
                        .fetch_seq_string(
                            contig,
                            start as usize,
                            end.saturating_sub(1) as usize,
                        )
                        .map_err(|e| MkError::HtsLibError(e))?
                };
                Ok(seq)
            }
        } else {
            Err(MkError::ContigMissing(contig.to_string()))
        }
    }

    fn has_sequence(&self, seq_name: &str) -> bool {
        self.contigs.contains_key(seq_name)
    }

    #[cfg(test)]
    fn sequence_fetch_stats(&self) -> (usize, u64) {
        (
            self.sequence_fetches.load(std::sync::atomic::Ordering::Relaxed),
            self.sequence_bases_fetched
                .load(std::sync::atomic::Ordering::Relaxed),
        )
    }
}

pub struct HtsLibFastaRecords {
    curr_record_idx: i32,
    total_records: u64,
    reader: faidx::Reader,
}

impl Iterator for HtsLibFastaRecords {
    type Item = MkResult<(String, String)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair().transpose()
    }
}

impl HtsLibFastaRecords {
    pub fn from_file(file: &PathBuf) -> MkResult<Self> {
        let reader = faidx::Reader::from_path(file)
            .map_err(|e| MkError::HtsLibError(e))?;
        let total_records = reader.n_seqs();
        let curr_record_idx = 0;
        Ok(Self { curr_record_idx, total_records, reader })
    }

    fn next_pair(&mut self) -> MkResult<Option<(String, String)>> {
        if self.curr_record_idx as u64 >= self.total_records {
            return Ok(None);
        }
        let contig_name = self
            .reader
            .seq_name(self.curr_record_idx)
            .map_err(|e| MkError::HtsLibError(e))?;
        let length = self.reader.fetch_seq_len(&contig_name);
        let seq = self
            .reader
            .fetch_seq_string(&contig_name, 0usize, length as usize)
            .map_err(|e| MkError::HtsLibError(e))?;
        self.curr_record_idx = self.curr_record_idx.saturating_add(1);
        Ok(Some((contig_name, seq)))
    }
}

#[derive(Clone)]
pub struct MotifLocationsLookup {
    reader: HtsFastaHandle,
    mask: bool,
    motifs: Vec<RegexMotif>,
    longest_motif_length: u64,
}

impl MotifLocationsLookup {
    pub(crate) fn has_sequence(&self, seq_name: &str) -> bool {
        self.reader.has_sequence(seq_name)
    }

    pub fn from_paths(
        fasta_fp: &PathBuf,
        mask: bool,
        _index_fp: Option<&PathBuf>,
        motifs: Vec<RegexMotif>,
        preload: bool,
    ) -> anyhow::Result<Self> {
        if motifs.is_empty() {
            bail!("motifs is empty, are you sure you want to make a lookup?");
        }
        let reader = HtsFastaHandle::from_file(fasta_fp, preload)?;
        let longest_motif_length =
            motifs.iter().map(|m| m.length() as u64).max().unwrap();

        Ok(Self { reader, motifs, mask, longest_motif_length })
    }

    #[inline]
    fn get_motifs_for_owner(
        &self,
        contig: &str,
        owner: std::ops::Range<u64>,
        tid: u32,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> MkResult<BitVec> {
        debug_assert!(!self.motifs.is_empty());
        let num_motifs = self.motifs.len();
        assert!(num_motifs < u8::MAX as usize);
        let bits_per_pos = self.motifs.len() * 2; // two strands
        let ref_end = *self.reader.contigs.get(contig).ok_or_else(|| {
            MkError::ContigMissing(format!("{contig} not in FASTA index"))
        })?;
        if owner.start > owner.end || owner.end > ref_end {
            return Err(MkError::InvalidReferenceCoordinates);
        }
        let owner_len = (owner.end - owner.start) as usize;
        let mut mask = bitvec![0; owner_len * bits_per_pos];
        if owner.is_empty() {
            return Ok(mask);
        }

        // Motif anchors belong to exactly one half-open owner interval, but
        // discovering an anchor may require bases on either side of it. Fetch
        // enough clipped context for the longest motif, then translate every
        // local hit back to its global reference coordinate before filtering
        // and writing it into the owner-local mask.
        let context = self.longest_motif_length.saturating_sub(1);
        let fetch_start = owner.start.saturating_sub(context);
        let fetch_end = owner.end.saturating_add(context).min(ref_end);
        let seq = self.reader.get_sequence(contig, fetch_start, fetch_end)?;
        let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };

        for (offset, motif) in self.motifs.iter().enumerate() {
            for (local_pos, strand) in find_motif_hits(&seq, motif) {
                let global_pos = fetch_start.saturating_add(local_pos as u64);
                if !owner.contains(&global_pos)
                    || stranded_position_filter.is_some_and(|spf| {
                        !spf.contains(tid as i32, global_pos, strand)
                    })
                {
                    continue;
                }
                let owner_pos = (global_pos - owner.start) as usize;
                match strand {
                    Strand::Positive => {
                        mask.set((owner_pos * bits_per_pos) + offset, true)
                    }
                    Strand::Negative => mask.set(
                        ((owner_pos * bits_per_pos) + num_motifs) + offset,
                        true,
                    ),
                }
            }
        }
        Ok(mask)
    }

    /// Find complete palindromic-motif anchor pairs owned by their positive
    /// anchor. A combined-strand output row is written at that positive
    /// coordinate, so the region and any stranded position filter must select
    /// that coordinate; the paired negative anchor is then admitted
    /// atomically even when a strand-specific filter does not select it.
    fn get_combined_motif_pairs(
        &self,
        contig: &str,
        positive_owner: std::ops::Range<u64>,
        owner_limit: u64,
        tid: u32,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> MkResult<Vec<(u64, u64, usize)>> {
        let ref_end = *self.reader.contigs.get(contig).ok_or_else(|| {
            MkError::ContigMissing(format!("{contig} not in FASTA index"))
        })?;
        if positive_owner.start > positive_owner.end
            || positive_owner.end > owner_limit
            || owner_limit > ref_end
        {
            return Err(MkError::InvalidReferenceCoordinates);
        }
        if positive_owner.is_empty() {
            return Ok(Vec::new());
        }

        let context = self.longest_motif_length.saturating_sub(1);
        let fetch_start = positive_owner.start.saturating_sub(context);
        let fetch_end = positive_owner.end.saturating_add(context).min(ref_end);
        let seq = self.reader.get_sequence(contig, fetch_start, fetch_end)?;
        let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };

        let mut pairs = Vec::new();
        for (motif_idx, motif) in self.motifs.iter().enumerate() {
            debug_assert!(motif.is_palendrome());
            let forward_offset = motif.forward_offset() as u64;
            let reverse_offset = motif.reverse_offset() as u64;
            for (local_pos, strand) in find_motif_hits(&seq, motif) {
                if strand != Strand::Positive {
                    continue;
                }
                let Some(positive_position) =
                    fetch_start.checked_add(local_pos as u64)
                else {
                    continue;
                };
                if !positive_owner.contains(&positive_position)
                    || stranded_position_filter.is_some_and(|spf| {
                        !spf.contains(
                            tid as i32,
                            positive_position,
                            Strand::Positive,
                        )
                    })
                {
                    continue;
                }
                let negative_position = if reverse_offset >= forward_offset {
                    positive_position
                        .checked_add(reverse_offset - forward_offset)
                } else {
                    positive_position
                        .checked_sub(forward_offset - reverse_offset)
                };
                let Some(negative_position) = negative_position else {
                    continue;
                };

                // Both anchors must be in the same selected reference region
                // and owner. In particular, never expose a reverse anchor
                // whose positive owner is before the region start.
                if negative_position < positive_owner.start
                    || negative_position >= owner_limit
                {
                    continue;
                }
                pairs.push((positive_position, negative_position, motif_idx));
            }
        }
        Ok(pairs)
    }

    fn get_motif_positions_combine_strands(
        &mut self,
        contig: &str,
        tid: u32,
        ref_end: u64,
        range: std::ops::Range<u64>,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> MkResult<(FocusPositions2, u32)> {
        let contig_end = *self.reader.contigs.get(contig).ok_or_else(|| {
            MkError::ContigMissing(format!("{contig} not in FASTA index"))
        })?;
        let owner_limit = ref_end.min(contig_end);
        if range.start > range.end || range.end > owner_limit {
            return Err(MkError::InvalidReferenceCoordinates);
        }
        let num_motifs = self.motifs.len();
        let bits_per_pos = num_motifs * 2; // two strands
        if range.is_empty() {
            let mask = bitvec![0; 0];
            return Ok((
                FocusPositions2::MotifMask {
                    mask,
                    num_motifs: self.num_motifs() as u8,
                },
                range.end as u32,
            ));
        }

        let mut end = range.end;
        let (mask, end) = loop {
            let pairs = self.get_combined_motif_pairs(
                contig,
                range.start..end,
                owner_limit,
                tid,
                stranded_position_filter,
            )?;
            let required_end = pairs
                .iter()
                .filter_map(|(positive, negative, _)| {
                    positive.max(negative).checked_add(1)
                })
                .max()
                .unwrap_or(end);
            if required_end > end {
                debug!(
                    "extending combined motif owner from {end} to \
                     {required_end}, contig={contig}, range={range:?}"
                );
                end = required_end;
                continue;
            }

            let mut mask =
                bitvec![0; (end - range.start) as usize * bits_per_pos];
            for (positive, negative, motif_idx) in pairs {
                let positive_local = (positive - range.start) as usize;
                let negative_local = (negative - range.start) as usize;
                mask.set((positive_local * bits_per_pos) + motif_idx, true);
                mask.set(
                    (negative_local * bits_per_pos) + num_motifs + motif_idx,
                    true,
                );
            }
            break (mask, end);
        };
        Ok((
            FocusPositions2::MotifMask {
                mask,
                num_motifs: self.num_motifs() as u8,
            },
            end as u32,
        ))
    }

    pub(crate) fn get_motif_positions(
        &mut self,
        contig: &str,
        tid: u32,
        ref_length: u32,
        range: std::ops::Range<u64>,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
        combine_strands: bool,
    ) -> MkResult<(FocusPositions2, u32)> {
        if combine_strands && !self.motifs.is_empty() {
            self.get_motif_positions_combine_strands(
                contig,
                tid,
                ref_length as u64,
                range,
                stranded_position_filter,
            )
        } else {
            let mask = self.get_motifs_for_owner(
                contig,
                range.clone(),
                tid,
                stranded_position_filter,
            )?;
            let focus_positions = FocusPositions2::MotifMask {
                mask,
                num_motifs: self.num_motifs() as u8,
            };
            Ok((focus_positions, range.end as u32))
        }
    }

    pub(super) fn num_motifs(&self) -> usize {
        self.motifs.len()
    }

    pub(crate) fn motif_infos(&self) -> impl Iterator<Item = &MotifInfo> + '_ {
        self.motifs.iter().map(|x| &x.motif_info)
    }
}

#[cfg(test)]
mod fasta_mod_tests {
    use std::{
        fs::File,
        io::Write,
        ops::Range,
        path::{Path, PathBuf},
    };

    use rand::prelude::{SeedableRng, StdRng};
    use rust_lapper::Lapper;
    use rustc_hash::FxHashMap;
    use rv::prelude::Rv;

    use crate::{
        fasta::{HtsFastaHandle, MotifLocationsLookup},
        interval_chunks::FocusPositions2,
        motifs::motif_bed::RegexMotif,
        position_filter::{Iv, StrandedPositionFilter},
        util::Strand,
    };

    fn write_fasta(root: &Path, sequence: &str) -> PathBuf {
        let fasta_path = root.join("reference.fa");
        File::create(&fasta_path)
            .unwrap()
            .write_all(format!(">chr1\n{sequence}\n").as_bytes())
            .unwrap();
        File::create(root.join("reference.fa.fai"))
            .unwrap()
            .write_all(
                format!(
                    "chr1\t{}\t6\t{}\t{}\n",
                    sequence.len(),
                    sequence.len(),
                    sequence.len() + 1
                )
                .as_bytes(),
            )
            .unwrap();
        fasta_path
    }

    fn decode_motif_mask(
        focus_positions: FocusPositions2,
        range: Range<u64>,
        num_motifs: usize,
    ) -> Vec<(u64, Strand, usize)> {
        let FocusPositions2::MotifMask { mask, num_motifs: observed } =
            focus_positions
        else {
            panic!("expected motif mask")
        };
        assert_eq!(observed as usize, num_motifs);
        let bits_per_pos = num_motifs * 2;
        assert_eq!(
            mask.len(),
            (range.end - range.start) as usize * bits_per_pos,
            "motif mask must cover exactly its owner interval"
        );
        let mut hits = Vec::new();
        for local_pos in 0..(range.end - range.start) as usize {
            for motif_idx in 0..num_motifs {
                if mask[(local_pos * bits_per_pos) + motif_idx] {
                    hits.push((
                        range.start + local_pos as u64,
                        Strand::Positive,
                        motif_idx,
                    ));
                }
                if mask[(local_pos * bits_per_pos) + num_motifs + motif_idx] {
                    hits.push((
                        range.start + local_pos as u64,
                        Strand::Negative,
                        motif_idx,
                    ));
                }
            }
        }
        hits
    }

    fn collect_motif_hits(
        fasta_path: &PathBuf,
        motifs: Vec<RegexMotif>,
        preload: bool,
        interval_size: u64,
        combine_strands: bool,
        position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> (Vec<(u64, Strand, usize)>, Vec<Range<u64>>) {
        let mut lookup = MotifLocationsLookup::from_paths(
            fasta_path, false, None, motifs, preload,
        )
        .unwrap();
        let contig_end = lookup.reader.contigs["chr1"];
        let num_motifs = lookup.num_motifs();
        let mut start = 0;
        let mut hits = Vec::new();
        let mut owners = Vec::new();
        while start < contig_end {
            let requested_end =
                start.saturating_add(interval_size).min(contig_end);
            let (focus_positions, actual_end) = lookup
                .get_motif_positions(
                    "chr1",
                    0,
                    contig_end as u32,
                    start..requested_end,
                    position_filter,
                    combine_strands,
                )
                .unwrap();
            let owner = start..actual_end as u64;
            assert!(owner.end > owner.start);
            hits.extend(decode_motif_mask(
                focus_positions,
                owner.clone(),
                num_motifs,
            ));
            owners.push(owner.clone());
            start = owner.end;
        }
        hits.sort_by_key(|(position, strand, motif_idx)| {
            (*position, *strand, *motif_idx)
        });
        (hits, owners)
    }

    #[test]
    fn fasta_subsequences_are_half_open_with_or_without_preload() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(temp_dir.path(), "ACGT");

        for preload in [false, true] {
            let reader =
                HtsFastaHandle::from_file(&fasta_path, preload).unwrap();
            for (start, end, expected) in [
                (0, 0, ""),
                (0, 1, "A"),
                (1, 3, "CG"),
                (3, 4, "T"),
                (4, 4, ""),
                (0, 4, "ACGT"),
            ] {
                assert_eq!(
                    reader.get_sequence("chr1", start, end).unwrap(),
                    expected,
                    "preload={preload}, range={start}..{end}"
                );
            }
            assert!(reader.get_sequence("chr1", 3, 2).is_err());
            assert!(reader.get_sequence("chr1", 0, 5).is_err());
        }
    }

    #[test]
    fn motif_hits_are_interval_and_preload_invariant() {
        struct Case<'a> {
            sequence: &'a str,
            motif: &'a str,
            offset: usize,
            expected: Vec<(u64, Strand, usize)>,
        }

        let cases = [
            Case {
                sequence: "ACGTCG",
                motif: "CG",
                offset: 0,
                expected: vec![
                    (1, Strand::Positive, 0),
                    (2, Strand::Negative, 0),
                    (4, Strand::Positive, 0),
                    (5, Strand::Negative, 0),
                ],
            },
            Case {
                sequence: "GATCNNGATC",
                motif: "GATC",
                offset: 1,
                expected: vec![
                    (1, Strand::Positive, 0),
                    (2, Strand::Negative, 0),
                    (7, Strand::Positive, 0),
                    (8, Strand::Negative, 0),
                ],
            },
            Case {
                sequence: "CGTACG",
                motif: "CGT",
                offset: 0,
                expected: vec![
                    (0, Strand::Positive, 0),
                    (5, Strand::Negative, 0),
                ],
            },
        ];

        for case in cases {
            let temp_dir = tempfile::tempdir().unwrap();
            let fasta_path = write_fasta(temp_dir.path(), case.sequence);
            for preload in [false, true] {
                for interval_size in [1, 2, 3, case.sequence.len() as u64, 100]
                {
                    let motif =
                        RegexMotif::parse_string(case.motif, case.offset)
                            .unwrap();
                    let (observed, _) = collect_motif_hits(
                        &fasta_path,
                        vec![motif],
                        preload,
                        interval_size,
                        false,
                        None,
                    );
                    assert_eq!(
                        observed, case.expected,
                        "sequence={}, motif={} {}, preload={preload}, \
                         interval_size={interval_size}",
                        case.sequence, case.motif, case.offset
                    );
                }
            }
        }
    }

    #[test]
    fn motif_position_filter_uses_global_stranded_anchors() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(temp_dir.path(), "NNGATCNN");
        let mut pos_positions = FxHashMap::default();
        pos_positions
            .insert(0, Lapper::new(vec![Iv { start: 3, stop: 4, val: () }]));
        let mut neg_positions = FxHashMap::default();
        neg_positions
            .insert(0, Lapper::new(vec![Iv { start: 4, stop: 5, val: () }]));
        let position_filter =
            StrandedPositionFilter { pos_positions, neg_positions };
        let motif = RegexMotif::parse_string("GATC", 1).unwrap();

        for preload in [false, true] {
            let (observed, _) = collect_motif_hits(
                &fasta_path,
                vec![motif.clone()],
                preload,
                1,
                false,
                Some(&position_filter),
            );
            assert_eq!(
                observed,
                vec![(3, Strand::Positive, 0), (4, Strand::Negative, 0),]
            );
        }
    }

    #[test]
    fn combined_strand_lookup_preserves_boundary_extension() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(temp_dir.path(), "GATCNNNN");
        let motif = RegexMotif::parse_string("GATC", 1).unwrap();

        for preload in [false, true] {
            let (observed, owners) = collect_motif_hits(
                &fasta_path,
                vec![motif.clone()],
                preload,
                2,
                true,
                None,
            );
            assert_eq!(owners, vec![0..3, 3..5, 5..7, 7..8]);
            assert_eq!(
                observed,
                vec![(1, Strand::Positive, 0), (2, Strand::Negative, 0),]
            );
        }
    }

    #[test]
    fn combined_strand_lookup_extends_to_the_actual_cgcg_pair() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(temp_dir.path(), "CGCGNN");
        let motif = RegexMotif::parse_string("CGCG", 0).unwrap();

        for preload in [false, true] {
            for (interval_size, expected_owners) in [
                (1, vec![0..4, 4..5, 5..6]),
                (2, vec![0..4, 4..6]),
                (3, vec![0..4, 4..6]),
                (6, vec![0..6]),
            ] {
                let (observed, owners) = collect_motif_hits(
                    &fasta_path,
                    vec![motif.clone()],
                    preload,
                    interval_size,
                    true,
                    None,
                );
                assert_eq!(owners, expected_owners);
                assert_eq!(
                    observed,
                    vec![(0, Strand::Positive, 0), (3, Strand::Negative, 0),],
                    "preload={preload}, interval_size={interval_size}"
                );
            }
        }
    }

    #[test]
    fn combined_strand_lookup_closes_over_every_new_boundary_pair() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(temp_dir.path(), "CGCGCGNN");
        let motif = RegexMotif::parse_string("CGCG", 0).unwrap();

        for preload in [false, true] {
            let (observed, owners) = collect_motif_hits(
                &fasta_path,
                vec![motif.clone()],
                preload,
                1,
                true,
                None,
            );
            assert_eq!(owners, vec![0..6, 6..7, 7..8]);
            assert_eq!(
                observed,
                vec![
                    (0, Strand::Positive, 0),
                    (2, Strand::Positive, 0),
                    (3, Strand::Negative, 0),
                    (5, Strand::Negative, 0),
                ]
            );
        }
    }

    #[test]
    fn combined_strand_dense_overlap_scan_work_is_bounded() {
        let temp_dir = tempfile::tempdir().unwrap();
        let repeat_count = 256usize;
        let sequence = format!("{}NN", "CG".repeat(repeat_count));
        let fasta_path = write_fasta(temp_dir.path(), &sequence);
        let motif = RegexMotif::parse_string("CGCG", 0).unwrap();
        let mut expected = (0..repeat_count - 1)
            .flat_map(|i| {
                [
                    (2 * i as u64, Strand::Positive, 0),
                    (2 * i as u64 + 3, Strand::Negative, 0),
                ]
            })
            .collect::<Vec<_>>();
        expected.sort_by_key(|(position, strand, motif_idx)| {
            (*position, *strand, *motif_idx)
        });

        for preload in [false, true] {
            let mut lookup = MotifLocationsLookup::from_paths(
                &fasta_path,
                false,
                None,
                vec![motif.clone()],
                preload,
            )
            .unwrap();
            let (focus_positions, actual_end) = lookup
                .get_motif_positions(
                    "chr1",
                    0,
                    sequence.len() as u32,
                    0..1,
                    None,
                    true,
                )
                .unwrap();

            assert_eq!(actual_end, (repeat_count * 2) as u32);
            assert_eq!(
                decode_motif_mask(focus_positions, 0..actual_end as u64, 1,),
                expected,
                "preload={preload}"
            );

            let (fetches, bases_fetched) = lookup.reader.sequence_fetch_stats();
            assert!(
                fetches <= 12,
                "dense overlap required {fetches} sequence fetches with \
                 preload={preload}"
            );
            assert!(
                bases_fetched <= 3 * sequence.len() as u64,
                "dense overlap fetched {bases_fetched} bases for a {}-base \
                 reference with preload={preload}",
                sequence.len()
            );
        }
    }

    #[test]
    fn combined_strand_lookup_owns_and_filters_pairs_by_positive_anchor() {
        let temp_dir = tempfile::tempdir().unwrap();
        let sequence = format!("{}CGCG{}", "N".repeat(38), "N".repeat(28));
        let fasta_path = write_fasta(temp_dir.path(), &sequence);
        let motif = RegexMotif::parse_string("CGCG", 0).unwrap();

        for preload in [false, true] {
            let mut lookup = MotifLocationsLookup::from_paths(
                &fasta_path,
                false,
                None,
                vec![motif.clone()],
                preload,
            )
            .unwrap();

            // The reverse anchor at 41 is not an independently owned hit when
            // its positive pair anchor at 38 is before the selected region.
            let (focus_positions, actual_end) = lookup
                .get_motif_positions("chr1", 0, 62, 40..62, None, true)
                .unwrap();
            assert_eq!(actual_end, 62);
            assert!(decode_motif_mask(focus_positions, 40..62, 1).is_empty());

            let mut pos_positions = FxHashMap::default();
            pos_positions.insert(
                0,
                Lapper::new(vec![Iv { start: 38, stop: 39, val: () }]),
            );
            let positive_only = StrandedPositionFilter {
                pos_positions,
                neg_positions: FxHashMap::default(),
            };
            let (focus_positions, actual_end) = lookup
                .get_motif_positions(
                    "chr1",
                    0,
                    62,
                    38..62,
                    Some(&positive_only),
                    true,
                )
                .unwrap();
            assert_eq!(actual_end, 62);
            assert_eq!(
                decode_motif_mask(focus_positions, 38..62, 1),
                vec![(38, Strand::Positive, 0), (41, Strand::Negative, 0),]
            );

            let mut neg_positions = FxHashMap::default();
            neg_positions.insert(
                0,
                Lapper::new(vec![Iv { start: 41, stop: 42, val: () }]),
            );
            let negative_only = StrandedPositionFilter {
                pos_positions: FxHashMap::default(),
                neg_positions,
            };
            let (focus_positions, actual_end) = lookup
                .get_motif_positions(
                    "chr1",
                    0,
                    62,
                    38..62,
                    Some(&negative_only),
                    true,
                )
                .unwrap();
            assert_eq!(actual_end, 62);
            assert!(decode_motif_mask(focus_positions, 38..62, 1).is_empty());
        }
    }

    #[test]
    fn test_hts_fasta_reader() {
        let compressed_fp = std::path::Path::new(
            "../tests/resources/CGI_ladder_3.6kb_ref.fa.gz",
        )
        .to_path_buf();
        let fp =
            std::path::Path::new("../tests/resources/CGI_ladder_3.6kb_ref.fa")
                .to_path_buf();
        let compressed_reader =
            HtsFastaHandle::from_file(&compressed_fp, false).unwrap();
        let reader = HtsFastaHandle::from_file(&fp, false).unwrap();
        assert_eq!(&compressed_reader.contigs, &reader.contigs);
        let mut rng = StdRng::seed_from_u64(42);
        for (contig, len) in compressed_reader.contigs.iter() {
            let half = len / 2;
            let dist =
                rv::prelude::DiscreteUniform::<u64>::new(0, half).unwrap();
            let start: u64 = dist.draw(&mut rng);
            let dist =
                rv::prelude::DiscreteUniform::<u64>::new(half + 1, *len - 1)
                    .unwrap();
            let end: u64 = dist.draw(&mut rng);
            assert!(start < end, "{start} >= {end}");
            assert_ne!(start, end);
            let subseq_compressed =
                compressed_reader.get_sequence(contig, start, end).unwrap();
            let subset_uncompressed =
                reader.get_sequence(contig, start, end).unwrap();
            assert_eq!(subseq_compressed, subset_uncompressed);
        }
    }
}
