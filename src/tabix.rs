use std::marker::PhantomData;
use std::ops::Range;
use std::path::PathBuf;

use anyhow::{anyhow, Context};
use itertools::Itertools;
use log_once::debug_once;
use rust_htslib::tbx::{Read, Reader as TbxReader};
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
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

    #[inline]
    fn fetch_region_it<'a>(
        &self,
        reader: &'a mut TbxReader,
        strand_rule: StrandRule,
    ) -> anyhow::Result<impl Iterator<Item = anyhow::Result<T>> + 'a> {
        Ok(reader
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
            .filter_ok(move |t| t.overlaps(strand_rule)))
    }

    fn get_reader(
        &self,
        chrom: &str,
        range: &Range<u64>,
    ) -> anyhow::Result<TbxReader> {
        let tid = *self
            .contigs
            .get(chrom)
            .ok_or(anyhow!("didn't find target-id for {chrom}"))?;
        let mut reader = TbxReader::from_path(&self.indexed_fp)?;
        reader.set_threads(4)?; // todo make param
        reader.fetch(tid, range.start, range.end).with_context(|| {
            format!("failed to fetch {chrom}:{}-{}", range.start, range.end)
        })?;
        Ok(reader)
    }

    pub(crate) fn fetch_region(
        &self,
        chrom: &str,
        range: &Range<u64>,
        strand_rule: StrandRule,
    ) -> anyhow::Result<Vec<T>> {
        let mut reader = self.get_reader(chrom, range)?;
        let it = self.fetch_region_it(&mut reader, strand_rule)?;
        it.collect()
    }

    pub fn get_contigs(&self) -> Vec<String> {
        self.contigs.keys().map(|x| x.to_owned()).collect()
    }
}
pub type BedMethylTbxIndex = HtsTabixHandler<BedMethylLine>;

impl HtsTabixHandler<BedMethylLine> {
    pub(crate) fn read_bedmethyl(
        &self,
        chrom: &str,
        range: &Range<u64>,
        min_coverage: u64,
        code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
    ) -> anyhow::Result<Vec<BedMethylLine>> {
        let mut reader = self.get_reader(chrom, range)?;
        let it = self.fetch_region_it(&mut reader, StrandRule::Both)?;
        it.filter_ok(|bml| bml.valid_coverage >= min_coverage)
            .filter_ok(|bml| {
                if code_lookup.contains_key(&bml.raw_mod_code) {
                    true
                } else {
                    debug_once!(
                        "encountered illegal mod-code for DMR, {}",
                        bml.raw_mod_code
                    );
                    false
                }
            })
            .collect()
    }
}
