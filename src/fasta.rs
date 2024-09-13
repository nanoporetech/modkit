use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::{bail, Context};
use bio::io::fasta::IndexedReader as FastaReader;
use log::debug;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use rustc_hash::FxHashMap;

use crate::find_motifs::motif_bed::{
    find_motif_hits, MotifLocations, MultipleMotifLocations, RegexMotif,
};
use crate::position_filter::StrandedPositionFilter;
use crate::util::StrandRule;

pub struct MotifLocationsLookup {
    reader: FastaReader<std::fs::File>,
    mask: bool,
    motifs: Vec<RegexMotif>,
    longest_motif_length: u64,
}

impl MotifLocationsLookup {
    pub fn from_paths(
        fasta_fp: &PathBuf,
        mask: bool,
        _index_fp: Option<&PathBuf>,
        motifs: Vec<RegexMotif>,
    ) -> anyhow::Result<Self> {
        if motifs.is_empty() {
            bail!("motifs is empty, are you sure you want to make a lookup?");
        }
        let reader = FastaReader::from_file(fasta_fp)?;
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
    ) -> MultipleMotifLocations {
        let motif_locations = self
            .motifs
            .par_iter()
            .map(|motif| {
                let positions = find_motif_hits(seq, &motif)
                    .into_par_iter()
                    .map(|(pos, strand)| (pos as u64 + start, strand))
                    .filter_map(|(pos, strand)| {
                        if let Some(position_filter) = stranded_position_filter
                        {
                            if position_filter.contains(tid as i32, pos, strand)
                            {
                                Some((pos as u32, strand))
                            } else {
                                None
                            }
                        } else {
                            Some((pos as u32, strand))
                        }
                    })
                    .fold(
                        || BTreeMap::<u32, StrandRule>::new(),
                        |mut acc, (pos, strand)| {
                            if let Some(strand_rule) = acc.get_mut(&pos) {
                                *strand_rule = strand_rule.absorb(strand);
                            } else {
                                acc.insert(pos, strand.into());
                            }
                            acc
                        },
                    )
                    .reduce(
                        || BTreeMap::<u32, StrandRule>::new(),
                        |a, b| a.into_iter().chain(b).collect(),
                    );
                let tid_to_motif_positions =
                    FxHashMap::from_iter([(tid, positions)]);
                MotifLocations::new(tid_to_motif_positions, motif.clone())
            })
            .collect::<Vec<MotifLocations>>();
        MultipleMotifLocations::new(motif_locations)
    }

    fn get_motif_positions_combine_strands(
        &mut self,
        contig: &str,
        tid: u32,
        ref_end: u64,
        range: std::ops::Range<u64>,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
    ) -> anyhow::Result<(MultipleMotifLocations, u32)> {
        let buffer_size = self.longest_motif_length * 5;
        let mut end = range.end;
        let mut end_w_buffer = std::cmp::min(range.end + buffer_size, ref_end);
        let mut too_close =
            end_w_buffer.saturating_sub(self.longest_motif_length);
        'fetch_loop: loop {
            self.reader.fetch(contig, range.start, end_w_buffer)?;
            let l = end_w_buffer
                .checked_sub(range.start)
                .expect("end should be >= start") as usize;
            let mut buff = Vec::<u8>::with_capacity(l);
            self.reader.read(&mut buff)?;
            buff.shrink_to_fit();
            debug_assert_eq!(buff.len(), l);
            let seq = String::from_utf8(buff)
                .context("got illegal characters in sequence")?;
            let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };
            let motif_locations = self.get_motifs_on_seq(
                &seq,
                range.start,
                tid,
                stranded_position_filter,
            );

            let motif_ivs = motif_locations
                .motif_locations
                .iter()
                .flat_map(|mls| {
                    let adj = mls
                        .motif_length()
                        .checked_sub(mls.motif().forward_offset())
                        .unwrap_or(mls.motif_length())
                        as u64;

                    let locations = mls.get_locations_unchecked(tid);
                    locations.iter().map(move |(pos, _)| Interval {
                        start: *pos as u64,
                        stop: (*pos as u64) + adj,
                        val: (),
                    })
                })
                .collect::<Vec<_>>();
            let intervals = {
                let mut tmp = Lapper::new(motif_ivs);
                tmp.merge_overlaps();
                tmp.set_cov();
                tmp
            };
            let search_end = if let Some(iv) =
                intervals.find(end.saturating_sub(1), end).next()
            {
                iv.stop
            } else {
                end
            };

            if (search_end < too_close) || (end_w_buffer >= ref_end) {
                let motif_locations = motif_locations
                    .motif_locations
                    .into_iter()
                    .map(|mls| {
                        let locations = mls
                            .tid_to_motif_positions
                            .into_iter()
                            .map(|(tid, poss)| {
                                let filt = poss
                                    .into_iter()
                                    .filter(|(p, _)| *p as u64 <= search_end)
                                    .collect::<_>();
                                (tid, filt)
                            })
                            .collect();
                        MotifLocations::new(locations, mls.motif)
                    })
                    .collect();
                return Ok((
                    MultipleMotifLocations::new(motif_locations),
                    search_end as u32,
                ));
            } else {
                debug!("too close, re-fetching");
                end = end_w_buffer;
                end_w_buffer += buffer_size;
                too_close =
                    end_w_buffer.saturating_sub(self.longest_motif_length);
                continue 'fetch_loop;
            }
        }
    }

    pub fn get_motif_positions(
        &mut self,
        contig: &str,
        tid: u32,
        ref_length: u32,
        range: std::ops::Range<u64>,
        stranded_position_filter: Option<&StrandedPositionFilter<()>>,
        combine_strands: bool,
    ) -> anyhow::Result<(MultipleMotifLocations, u32)> {
        if combine_strands && !self.motifs.is_empty() {
            self.get_motif_positions_combine_strands(
                contig,
                tid,
                ref_length as u64,
                range,
                stranded_position_filter,
            )
        } else {
            self.reader.fetch(contig, range.start, range.end)?;
            let l = range
                .end
                .checked_sub(range.start)
                .expect("end should be >= start") as usize;
            let mut buff = Vec::<u8>::with_capacity(l);
            self.reader.read(&mut buff)?;
            buff.shrink_to_fit();
            debug_assert_eq!(buff.len(), l);
            let seq = String::from_utf8(buff)
                .context("got illegal characters in sequence")?;
            let seq = if self.mask { seq } else { seq.to_ascii_uppercase() };
            let multiple_motif_locations = self.get_motifs_on_seq(
                &seq,
                range.start,
                tid,
                stranded_position_filter,
            );
            Ok((multiple_motif_locations, range.end as u32))
        }
    }
}
