use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::find_motifs::{EnrichedMotif, KmerRef, RawBase};
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;

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

pub(super) fn log_odds<
    T: num_traits::Num + num_traits::cast::AsPrimitive<f32>,
>(
    low_pos: T,
    low_neg: T,
    high_pos: T,
    high_neg: T,
) -> f32 {
    let numer = high_pos * low_neg;
    let denom = low_pos * high_neg;
    if denom == T::zero() {
        if numer == T::zero() {
            0f32
        } else {
            f32::INFINITY
        }
    } else if numer == T::zero() {
        f32::NEG_INFINITY
    } else {
        let numer: f32 = numer.as_();
        let denom: f32 = denom.as_();
        (numer / denom).log2()
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
