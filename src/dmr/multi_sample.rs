use std::collections::HashMap;
use std::path::PathBuf;

use crate::mod_base_code::DnaBase;
use anyhow::bail;
use bio::io::fasta::Reader as FastaReader;
use derive_new::new;
use indicatif::{MultiProgress, ProgressIterator};
use log::debug;
use rayon::prelude::*;
use rustc_hash::FxHashSet;

use crate::util::{get_ticker, Strand};

fn factorial(n: usize) -> anyhow::Result<usize> {
    if n < 1 {
        bail!("n must be > 1")
    } else {
        Ok((1..=n).product())
    }
}

pub(super) fn n_choose_2(n: usize) -> anyhow::Result<usize> {
    if n < 2 {
        bail!("n must be > 2")
    } else if n == 2 {
        Ok(1)
    } else {
        let numerator = factorial(n)?;
        let denom = 2 * factorial(n - 2)?;
        Ok(numerator / denom)
    }
}

#[derive(new, Debug)]
pub(super) struct DmrSample {
    pub(super) bedmethyl_fp: PathBuf,
    pub(super) index: PathBuf,
    pub(super) name: String,
}

pub(super) fn get_reference_modified_base_positions(
    ref_fasta: &PathBuf,
    modified_bases: &[DnaBase],
    mask: bool,
    mpb: MultiProgress,
) -> anyhow::Result<(
    HashMap<String, Vec<(DnaBase, u64)>>,
    HashMap<String, Vec<(DnaBase, u64)>>,
)> {
    // this should probably be a struct...
    let fasta_reader = FastaReader::from_file(ref_fasta)?;
    let reader_pb = mpb.add(get_ticker());
    reader_pb.set_message("sequences read");
    let positions_pb = mpb.add(get_ticker());
    positions_pb.set_message("positions found");

    let (snd, rcv) = crossbeam_channel::unbounded();
    std::thread::spawn(move || {
        fasta_reader
            .records()
            .progress_with(reader_pb)
            .filter_map(|r| r.ok())
            .filter_map(|record| {
                String::from_utf8(record.seq().to_vec())
                    .map(|s| if mask { s } else { s.to_ascii_uppercase() })
                    .ok()
                    .map(|s| (s, record.id().to_string()))
            })
            .for_each(|(seq, name)| match snd.send((seq, name)) {
                Ok(_) => {}
                Err(e) => {
                    debug!("failed to send seq on channel, {}", e.to_string());
                }
            });
    });
    let pos_bases = modified_bases
        .iter()
        .map(|b| b.char())
        .collect::<FxHashSet<char>>();
    let neg_bases = modified_bases
        .iter()
        .map(|b| b.complement().char())
        .collect::<FxHashSet<char>>();
    let (positive_hits, negative_hits) = rcv
        .iter()
        .par_bridge()
        .map(|(seq, name)| {
            let modified_positions = seq
                .chars()
                .enumerate()
                .filter_map(|(i, c)| {
                    if pos_bases.contains(&c) {
                        Some((i, Strand::Positive, c))
                    } else if neg_bases.contains(&c) {
                        Some((i, Strand::Negative, c))
                    } else {
                        None
                    }
                })
                .collect::<Vec<(usize, Strand, char)>>();
            positions_pb.inc(modified_positions.len() as u64);
            (name, modified_positions)
        })
        .fold(
            || (HashMap::new(), HashMap::new()),
            |(mut positive_positions, mut negative_positions),
             (name, positions)| {
                for (position, strand, raw_base) in positions {
                    let dna_base = DnaBase::parse(raw_base).unwrap();
                    match strand {
                        Strand::Positive => {
                            positive_positions
                                .entry(name.clone())
                                .or_insert(Vec::new())
                                .push((dna_base, position as u64));
                        }
                        Strand::Negative => {
                            negative_positions
                                .entry(name.clone())
                                .or_insert(Vec::new())
                                .push((dna_base, position as u64));
                        }
                    }
                }
                (positive_positions, negative_positions)
            },
        )
        .reduce(
            || (HashMap::new(), HashMap::new()),
            |(a_pos, a_neg), (b_pos, b_neg)| {
                let pos = a_pos.into_iter().chain(b_pos).collect();
                let neg = a_neg.into_iter().chain(b_neg).collect();
                (pos, neg)
            },
        );

    Ok((positive_hits, negative_hits))
}
