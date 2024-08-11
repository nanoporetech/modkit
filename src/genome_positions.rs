use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::Range;
use std::path::PathBuf;

use bio::io::fasta::Reader as FastaReader;
use indicatif::{MultiProgress, ProgressIterator};
use log::debug;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::mod_base_code::DnaBase;
use crate::util::{get_ticker, Strand, StrandRule};

#[derive(Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub(crate) struct StrandedPosition<T>
where
    T: Copy + Debug + Ord + PartialOrd + Hash + Eq + PartialEq,
{
    pub(crate) position: u64,
    pub(crate) strand: Strand,
    pub(crate) value: T,
}

/// A struct to quickly check if a specific position on a strand should be used.
/// Assumes that a position can only be on either the (+) or (-) strand
pub(crate) struct GenomePositions {
    /// these are the bases we care to compare, i.e. C and A for 5mC and
    /// 6mA, respectively.
    pub positive_strand_bases: FxHashSet<char>,
    pub negative_strand_bases: FxHashSet<char>,
    /// this is the reference genome - we'll search it on the fly to
    /// reduce memory consumption.
    contigs: FxHashMap<String, Vec<char>>,
}

impl GenomePositions {
    pub(super) fn new_from_sequences(
        bases: &[DnaBase],
        fasta_fp: &PathBuf,
        mask: bool,
        all_contigs: &HashSet<String>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Self> {
        let fasta_reader = FastaReader::from_file(&fasta_fp)?;
        let reader_pb = multi_progress.add(get_ticker());
        reader_pb.set_message("sequences read");
        let pos_bases =
            bases.iter().map(|b| b.char()).collect::<FxHashSet<char>>();
        let neg_bases = bases
            .iter()
            .map(|b| b.complement().char())
            .collect::<FxHashSet<char>>();

        let contigs = fasta_reader
            .records()
            .progress_with(reader_pb)
            .inspect(|r| {
                if let Err(e) = r {
                    debug!("failed to parse FASTA sequence, {e}")
                }
            })
            .filter_map(|res| res.ok())
            .filter(|record| all_contigs.contains(record.id()))
            .map(|record| {
                let contig_name = record.id().to_string();
                let seq = record
                    .seq()
                    .iter()
                    .map(|b| {
                        let base = char::from(*b);
                        if mask {
                            base
                        } else {
                            base.to_ascii_uppercase()
                        }
                    })
                    .collect::<Vec<char>>();

                (contig_name, seq)
            })
            .collect::<FxHashMap<String, Vec<char>>>();

        Ok(Self {
            positive_strand_bases: pos_bases,
            negative_strand_bases: neg_bases,
            contigs,
        })
    }

    pub(crate) fn get_positions(
        &self,
        chrom_name: &str,
        dmr_interval: &Range<u64>,
        strand_rule: StrandRule,
    ) -> Option<Vec<StrandedPosition<DnaBase>>> {
        let interval =
            (dmr_interval.start as usize)..(dmr_interval.end as usize);
        self.contigs.get(chrom_name).map(|seq| {
            seq[interval.start..interval.end]
                .iter()
                .enumerate()
                .filter_map(|(i, base)| {
                    let position = i + interval.start;
                    if self.positive_strand_bases.contains(base)
                        && strand_rule.covers(Strand::Positive)
                    {
                        Some(StrandedPosition {
                            position: position as u64,
                            strand: Strand::Positive,
                            value: DnaBase::parse(*base).unwrap(),
                        })
                    } else if self.negative_strand_bases.contains(base)
                        && strand_rule.covers(Strand::Negative)
                    {
                        Some(StrandedPosition {
                            position: position as u64,
                            strand: Strand::Negative,
                            value: DnaBase::parse(*base).unwrap(),
                        })
                    } else {
                        None
                    }
                })
                .collect::<Vec<StrandedPosition<DnaBase>>>()
        })
    }

    pub(crate) fn contig_sizes(
        &self,
    ) -> impl Iterator<Item = (&String, usize)> {
        self.contigs.iter().map(|(name, seq)| (name, seq.len()))
    }
}
