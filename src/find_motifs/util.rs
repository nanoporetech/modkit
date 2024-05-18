use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::find_motifs::{EnrichedMotif, KmerRef, RawBase};
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;

// pub(super) fn string_kmer(kmer: &[u8]) -> String {
//     kmer.iter().map(|&x| x as char).collect::<String>()
// }
//
// pub(super) fn pretty_kmer_ref_table<T: Display>(
//     kmer_lookup: &FxHashMap<KmerRef, T>,
// ) -> Table {
//     let mut table = Table::new();
//     let kmer_iter = kmer_lookup
//         .iter()
//         .map(|(kmer, x)| {
//             let nt_kmer = kmer.iter().map(|&x| x as
// char).collect::<String>();             (nt_kmer, x)
//         })
//         .sorted_by(|(a, _), (b, _)| a.cmp(b));
//     for (nt_kmer, x) in kmer_iter {
//         table.add_row(row![nt_kmer, x]);
//     }
//
//     table
// }

#[inline]
pub(super) fn aggregate_base_counts_on_position(
    kmer_counts: &[(KmerRef, u32)],
) -> FxHashMap<RawBase, FxHashMap<usize, u32>> {
    kmer_counts
        .par_iter()
        .map(|(kmer, count)| {
            kmer.iter().enumerate().fold(
                FxHashMap::<u8, FxHashMap<usize, u32>>::default(),
                |mut counter, (pos, nt)| {
                    *counter
                        .entry(*nt)
                        .or_insert(FxHashMap::default())
                        .entry(pos)
                        .or_insert(0u32) += *count;
                    counter
                },
            )
        })
        .reduce(|| FxHashMap::default(), |a, b| a.op(b))
}

pub(super) fn log_odds(
    low_pos: u32,
    low_neg: u32,
    high_pos: u32,
    high_neg: u32,
) -> f32 {
    let numer = high_pos * low_neg;
    let denom = low_pos * high_neg;
    if denom == 0u32 {
        if numer == 0u32 {
            0f32
        } else {
            f32::INFINITY
        }
    } else if numer == 0u32 {
        f32::NEG_INFINITY
    } else {
        (numer as f32 / denom as f32).log2()
    }
}

pub(super) fn log_motifs(
    mod_code: ModCodeRepr,
    motifs: &[EnrichedMotif],
    kind: &str,
    duration: std::time::Duration,
) -> String {
    if motifs.is_empty() {
        format!(
            "({mod_code}) {kind} didn't find any motifs, took {}ms",
            duration.as_millis()
        )
    } else {
        let mut message = format!(
            "({mod_code}) {kind} ({}), took {}ms:\n",
            motifs.len(),
            duration.as_millis()
        );

        for motif in motifs.iter().take(motifs.len() - 1usize) {
            message.push_str(&format!("\t{motif}\n"));
        }
        message.push_str(&format!("\t{}", motifs.last().unwrap()));

        format!("{message}")
    }
}
