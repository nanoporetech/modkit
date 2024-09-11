use std::marker::PhantomData;
use std::ops::Range;
use std::path::PathBuf;

use anyhow::{anyhow, Context};
use itertools::Itertools;
use rust_htslib::tbx::{Read, Reader as TbxReader};
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::util::StrandRule;

pub(crate) trait ParseBedLine {
    fn parse(l: &str) -> anyhow::Result<Self>
    where
        Self: Sized;
    fn overlaps(&self, strand_rule: StrandRule) -> bool;
}

impl ParseBedLine for BedMethylLine {
    fn parse(l: &str) -> anyhow::Result<Self> {
        BedMethylLine::parse(l)
    }

    fn overlaps(&self, strand_rule: StrandRule) -> bool {
        self.strand.overlaps(&strand_rule)
    }
}

pub(crate) struct HtsTabixHandler<T: ParseBedLine> {
    indexed_fp: PathBuf,
    contigs: FxHashMap<String, u64>,
    _t: PhantomData<T>,
}

impl<T: ParseBedLine> HtsTabixHandler<T> {
    pub(crate) fn from_path(path: &PathBuf) -> anyhow::Result<Self> {
        let reader = TbxReader::from_path(path)?;
        let contigs = reader
            .seqnames()
            .into_iter()
            .map(|n| reader.tid(&n).map(|tid| (n, tid)))
            .collect::<Result<FxHashMap<String, u64>, _>>()
            .context(
                "failed to collect contig IDs and names, invalid tabix header?",
            )?;
        Ok(Self { indexed_fp: path.to_owned(), contigs, _t: PhantomData })
    }

    pub(crate) fn has_contig(&self, contig: &str) -> bool {
        self.contigs.contains_key(contig)
    }

    pub(crate) fn fetch_region(
        &self,
        chrom: &str,
        range: &Range<u64>,
        strand_rule: StrandRule,
    ) -> anyhow::Result<Vec<T>> {
        let tid = *self
            .contigs
            .get(chrom)
            .ok_or(anyhow!("didn't find target-id for {chrom}"))?;
        let mut reader = TbxReader::from_path(&self.indexed_fp)?;
        reader.set_threads(4)?; // todo make param
        reader.fetch(tid, range.start, range.end).with_context(|| {
            format!("failed to fetch {chrom}:{}-{}", range.start, range.end)
        })?;

        reader
            .records()
            .map(|r| {
                r.map_err(|e| anyhow!("tbx failed to read, {e}"))
                    .and_then(|bs| {
                        String::from_utf8(bs).map_err(|e| {
                            anyhow!("failed to convert record to UTF8, {e}")
                        })
                    })
                    .and_then(|s| T::parse(&s))
            })
            .filter_ok(|t| t.overlaps(strand_rule))
            .collect::<anyhow::Result<Vec<T>>>()
    }
}
