use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use derive_new::new;
use itertools::Itertools;
use rustc_hash::FxHashMap;

use crate::mod_bam::{BaseModCall, BaseModProbs};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::read_ids_to_base_mod_probs::{
    ModProfile, ReadBaseModProfile, ReadsBaseModProfile,
};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util;
use crate::util::{
    create_out_directory, get_reference_mod_strand, Kmer, Strand,
};
use crate::writers::TsvWriter;

#[derive(new)]
pub(crate) struct PositionModCalls {
    query_position: usize,
    pub(crate) ref_position: Option<i64>,
    num_soft_clipped_start: usize,
    num_soft_clipped_end: usize,
    read_length: usize,
    pub(crate) base_mod_probs: BaseModProbs,
    q_base: u8,
    query_kmer: Kmer,
    pub(crate) mod_strand: Strand,
    pub(crate) alignment_strand: Option<Strand>,
    pub(crate) canonical_base: DnaBase,
}

impl PositionModCalls {
    fn header() -> String {
        let tab = '\t';
        format!(
            "\
            read_id{tab}\
            forward_read_position{tab}\
            ref_position{tab}\
            chrom{tab}\
            mod_strand{tab}\
            ref_strand{tab}\
            ref_mod_strand{tab}\
            fw_soft_clipped_start{tab}\
            fw_soft_clipped_end{tab}\
            read_length{tab}\
            call_prob{tab}\
            call_code{tab}\
            base_qual{tab}\
            ref_kmer{tab}\
            query_kmer{tab}\
            canonical_base{tab}\
            modified_primary_base{tab}\
            fail{tab}\
            inferred{tab}\
            within_alignment{tab}\
            flag"
        )
    }

    pub(crate) fn from_profile(
        read_base_mod_profile: &ReadBaseModProfile,
    ) -> Vec<Self> {
        type Key = (usize, Strand, DnaBase);
        let (grouped, mod_codes): (
            HashMap<Key, Vec<&ModProfile>>,
            HashSet<ModCodeRepr>,
        ) = read_base_mod_profile.iter_profiles().fold(
            (HashMap::new(), HashSet::new()),
            |(mut acc, mut codes), x| {
                let k = (x.query_position, x.mod_strand, x.canonical_base);
                acc.entry(k).or_insert(Vec::new()).push(x);
                codes.insert(x.raw_mod_code);
                (acc, codes)
            },
        );
        let mod_codes = mod_codes.into_iter().collect::<Vec<ModCodeRepr>>();

        grouped
            .into_iter()
            .fold(
                Vec::<Self>::new(),
                |mut acc, ((query_pos, strand, base), mod_profile)| {
                    let base_mod_probs = if mod_profile
                        .iter()
                        .any(|x| x.inferred)
                    {
                        BaseModProbs::new_inferred_canonical(mod_codes.iter())
                    } else {
                        let mut probs = mod_profile
                            .iter()
                            .map(|x| (x.raw_mod_code, x.q_mod))
                            .collect::<FxHashMap<ModCodeRepr, f32>>();
                        for code in mod_codes.iter() {
                            if !probs.contains_key(&code) {
                                probs.insert(*code, 0f32);
                            }
                        }

                        BaseModProbs::new(probs, false)
                    };
                    let template = &mod_profile[0];
                    let ref_position = template.ref_position;
                    let num_clip_start = template.num_soft_clipped_start;
                    let num_clip_end = template.num_soft_clipped_end;
                    let q_base = template.q_base;
                    let kmer = template.query_kmer;
                    let alignment_strand = template.alignment_strand;

                    let pos_mod_calls = PositionModCalls::new(
                        query_pos,
                        ref_position,
                        num_clip_start,
                        num_clip_end,
                        template.read_length,
                        base_mod_probs,
                        q_base,
                        kmer,
                        strand,
                        alignment_strand,
                        base,
                    );
                    acc.push(pos_mod_calls);

                    acc
                },
            )
            .into_iter()
            .sorted_by(|a, b| {
                if a.alignment_strand
                    .map(|s| s == Strand::Negative)
                    .unwrap_or(false)
                {
                    b.query_position.cmp(&a.query_position)
                } else {
                    a.query_position.cmp(&b.query_position)
                }
            })
            .collect()
    }

    fn within_alignment(&self) -> bool {
        util::within_alignment(
            self.query_position,
            self.num_soft_clipped_start,
            self.num_soft_clipped_end,
            self.read_length,
        )
    }

    pub(crate) fn to_row(
        &self,
        read_id: &str,
        chrom_name: Option<&String>,
        caller: &MultipleThresholdModCaller,
        reference_seqs: &HashMap<String, Vec<u8>>,
        flag: u16,
    ) -> String {
        let tab = '\t';
        let missing = ".".to_string();
        let chrom_name_label = chrom_name.unwrap_or(&missing).to_owned();
        let forward_read_position = self.query_position;
        let ref_position = self.ref_position.unwrap_or(-1);
        let mod_strand = self.mod_strand.to_char();
        let ref_strand =
            self.alignment_strand.map(|x| x.to_char()).unwrap_or('.');
        let ref_mod_strand = self
            .alignment_strand
            .map(|x| get_reference_mod_strand(self.mod_strand, x).to_char())
            .unwrap_or('.');
        let fw_soft_clipped_start = self.num_soft_clipped_start;
        let fw_soft_clipped_end = self.num_soft_clipped_end;
        let (mod_call_prob, mod_call_code) =
            match self.base_mod_probs.argmax_base_mod_call() {
                BaseModCall::Canonical(p) => (p, "-".to_string()),
                BaseModCall::Modified(p, code) => (p, code.to_string()),
                BaseModCall::Filtered => {
                    unreachable!("argmax should not output filtered calls")
                }
            };
        let read_length = self.read_length;
        let base_qual = self.q_base;
        let query_kmer = format!("{}", self.query_kmer);
        let ref_kmer = if let Some(ref_pos) = self.ref_position {
            if ref_pos < 0 {
                None
            } else {
                reference_seqs.get(&chrom_name_label).map(|s| {
                    Kmer::from_seq(s, ref_pos as usize, self.query_kmer.size)
                        .to_string()
                })
            }
        } else {
            None
        };
        let ref_kmer_rep = ref_kmer.as_ref().unwrap_or(&missing);
        let canonical_base = self.canonical_base.char();
        let modified_primary_base = if self.mod_strand == Strand::Negative {
            self.canonical_base.complement().char()
        } else {
            self.canonical_base.char()
        };
        let filtered = caller.call(&self.canonical_base, &self.base_mod_probs)
            == BaseModCall::Filtered;
        let inferred = self.base_mod_probs.inferred;
        let within_alignment = chrom_name.is_some() && self.within_alignment();

        format!(
            "\
            {read_id}{tab}\
            {forward_read_position}{tab}\
            {ref_position}{tab}\
            {chrom_name_label}{tab}\
            {mod_strand}{tab}\
            {ref_strand}{tab}\
            {ref_mod_strand}{tab}\
            {fw_soft_clipped_start}{tab}\
            {fw_soft_clipped_end}{tab}\
            {read_length}{tab}\
            {mod_call_prob}{tab}\
            {mod_call_code}{tab}\
            {base_qual}{tab}\
            {ref_kmer_rep}{tab}\
            {query_kmer}{tab}\
            {canonical_base}{tab}\
            {modified_primary_base}{tab}\
            {filtered}{tab}\
            {inferred}{tab}\
            {within_alignment}{tab}\
            {flag}\n"
        )
    }
}

pub trait OutwriterWithMemory<T> {
    fn write(&mut self, item: T, kmer_size: usize) -> anyhow::Result<u64>;
    fn num_reads(&self) -> usize;
}

pub struct TsvWriterWithContigNames<W: Write> {
    tsv_writer: TsvWriter<W>,
    tid_to_name: HashMap<u32, String>,
    name_to_seq: HashMap<String, Vec<u8>>,
    number_of_written_reads: usize,
    read_calls_writer: Option<TsvWriter<File>>,
    caller: MultipleThresholdModCaller,
}

impl<W: Write> TsvWriterWithContigNames<W> {
    pub(crate) fn new(
        output_writer: TsvWriter<W>,
        tid_to_name: HashMap<u32, String>,
        name_to_seq: HashMap<String, Vec<u8>>,
        read_calls_path: Option<&PathBuf>,
        caller: MultipleThresholdModCaller,
        force: bool,
    ) -> anyhow::Result<Self> {
        let read_calls_writer = read_calls_path
            .map(|fp| {
                create_out_directory(fp)?;
                TsvWriter::new_path(fp, force, Some(PositionModCalls::header()))
            })
            .transpose()?;
        Ok(Self {
            tsv_writer: output_writer,
            tid_to_name,
            name_to_seq,
            number_of_written_reads: 0,
            read_calls_writer,
            caller,
        })
    }
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W>
{
    fn write(
        &mut self,
        item: ReadsBaseModProfile,
        kmer_size: usize,
    ) -> anyhow::Result<u64> {
        let missing_chrom = ".".to_string();
        let mut rows_written = 0u64;
        for profile in item.profiles.iter() {
            let chrom_name = profile
                .chrom_id
                .and_then(|chrom_id| self.tid_to_name.get(&chrom_id));
            for mod_profile in profile.iter_profiles() {
                let row = mod_profile.to_row(
                    &profile.record_name,
                    chrom_name.unwrap_or(&missing_chrom),
                    &self.name_to_seq,
                    kmer_size,
                    profile.flag,
                );
                self.tsv_writer.write(row.as_bytes())?;
                rows_written += 1;
            }
            if let Some(read_calls_writer) = self.read_calls_writer.as_mut() {
                let position_calls = PositionModCalls::from_profile(&profile);
                for call in position_calls {
                    read_calls_writer.write(
                        call.to_row(
                            &profile.record_name,
                            chrom_name,
                            &self.caller,
                            &self.name_to_seq,
                            profile.flag,
                        )
                        .as_bytes(),
                    )?;
                }
            }
            self.number_of_written_reads += 1;
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.number_of_written_reads
    }
}
