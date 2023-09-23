use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::bail;
use bio::io::fasta::Reader as FastaReader;
use derive_new::new;
use indicatif::{MultiProgress, ProgressIterator};
use log::debug;
use rayon::prelude::*;

use crate::motif_bed::{find_motif_hits, RegexMotif};
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
    modified_bases: &[RegexMotif],
    mask: bool,
    mpb: MultiProgress,
) -> anyhow::Result<(HashMap<String, Vec<u64>>, HashMap<String, Vec<u64>>)> {
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
    let (positive_hits, negative_hits) = rcv
        .iter()
        .par_bridge()
        .map(|(seq, name)| {
            let modified_positions = modified_bases
                .iter()
                .flat_map(|base| find_motif_hits(&seq, base))
                .collect::<Vec<(usize, Strand)>>();
            positions_pb.inc(modified_positions.len() as u64);
            (name, modified_positions)
        })
        .fold(
            || (HashMap::new(), HashMap::new()),
            |(mut positive_positions, mut negative_positions),
             (name, positions)| {
                for (position, strand) in positions {
                    match strand {
                        Strand::Positive => {
                            positive_positions
                                .entry(name.clone())
                                .or_insert(Vec::new())
                                .push(position as u64);
                        }
                        Strand::Negative => {
                            negative_positions
                                .entry(name.clone())
                                .or_insert(Vec::new())
                                .push(position as u64);
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
