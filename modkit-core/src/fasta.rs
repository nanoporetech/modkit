use std::path::PathBuf;

use anyhow::{bail, Context};
use bitvec::bitvec;
use bitvec::vec::BitVec;
use itertools::Itertools;
use log::debug;
use rayon::prelude::*;
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
        })
    }

    fn get_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> MkResult<String> {
        if let Some(length) = self.contigs.get(contig) {
            if end > *length {
                Err(MkError::InvalidReferenceCoordinates)
            } else {
                let seq = if self.preloaded
                    && self.sequences.get(contig).is_some()
                {
                    let s = self.sequences.get(contig).unwrap();
                    s.substring(start as usize, end as usize).to_string()
                } else {
                    let tmp_reader = faidx::Reader::from_path(&self.fasta_fp)
                        .map_err(|e| MkError::HtsLibError(e))?;
                    tmp_reader
                        .fetch_seq_string(contig, start as usize, end as usize)
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
    fn get_motifs_on_seq(
        &self,
        seq: &str,
        start: u64,
        tid: u32,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> BitVec {
        debug_assert!(!self.motifs.is_empty());
        let num_motifs = self.motifs.len();
        assert!(num_motifs < u8::MAX as usize);
        let bits_per_pos = self.motifs.len() * 2; // two strands
        let mut mask = bitvec![0; seq.len() * bits_per_pos];

        for (offset, motif) in self.motifs.iter().enumerate() {
            let positions = if let Some(spf) = stranded_position_filter {
                find_motif_hits(seq, motif)
                    .into_par_iter()
                    .filter(|(pos, strand)| {
                        spf.contains(
                            tid as i32,
                            (*pos as u64).saturating_add(start),
                            *strand,
                        )
                    })
                    .collect::<Vec<(usize, Strand)>>()
            } else {
                find_motif_hits(seq, motif)
            };
            for (pos, strand) in positions {
                match strand {
                    Strand::Positive => {
                        mask.set((pos * bits_per_pos) + offset, true)
                    }
                    Strand::Negative => mask.set(
                        ((pos * bits_per_pos) + num_motifs) + offset,
                        true,
                    ),
                }
            }
        }
        mask
    }

    fn get_motif_positions_combine_strands(
        &mut self,
        contig: &str,
        tid: u32,
        _ref_end: u64,
        range: std::ops::Range<u64>,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> MkResult<(FocusPositions2, u32)> {
        let ref_end = *self.reader.contigs.get(contig).ok_or_else(|| {
            MkError::ContigMissing(format!("{contig} not in FASTA index"))
        })?;
        let buffer_size = self.longest_motif_length * 5;
        let num_motifs = self.motifs.len();
        let bits_per_pos = (num_motifs * 2) as u64; // two strands
        let mut end = range.end;
        let mut end_w_buffer = std::cmp::min(range.end + buffer_size, ref_end);
        let mask = 'fetch_loop: loop {
            let seq =
                self.reader.get_sequence(contig, range.start, end_w_buffer)?;
            let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };
            let mask = self.get_motifs_on_seq(
                &seq,
                range.start,
                tid,
                stranded_position_filter,
            );
            let mut end_idx = ((end - range.start - 1) * bits_per_pos) as usize;
            assert!(
                (end_idx + num_motifs) < mask.len(),
                "off the end {} {}",
                end_idx + num_motifs,
                mask.len()
            );

            while mask[end_idx..(end_idx + num_motifs)].any() {
                debug!("end_idx ({end_idx}) hits a motif.. end={end}",);
                end = std::cmp::min(
                    end.saturating_add(self.longest_motif_length),
                    ref_end,
                );
                end_idx = ((end - range.start - 1) * bits_per_pos) as usize;
                debug!(
                    "moved end_idx to {end_idx}, end={end} ref_end={ref_end}, \
                     contig={contig}, range={range:?}"
                );
                debug_assert!(end_idx + num_motifs < mask.len());
                if end >= end_w_buffer {
                    debug!("too close, re-fetching..");
                    end_w_buffer =
                        std::cmp::min(end.saturating_add(buffer_size), ref_end);
                    continue 'fetch_loop;
                }
                debug!("re-check..");
            }
            break mask;
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
            let seq =
                self.reader.get_sequence(contig, range.start, range.end)?;
            let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };
            let mask = self.get_motifs_on_seq(
                &seq,
                range.start,
                tid,
                stranded_position_filter,
            );
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
    use crate::fasta::HtsFastaHandle;
    use rand::prelude::{SeedableRng, StdRng};
    use rv::prelude::Rv;

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
