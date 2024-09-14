use std::collections::HashMap;
use std::io::Write;

use crate::mod_bam::BaseModCall;
use crate::read_ids_to_base_mod_probs::{
    PositionModCalls, ReadsBaseModProfile,
};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{get_reference_mod_strand, Kmer, Strand};
use crate::writers::TsvWriter;

impl PositionModCalls {
    pub(super) fn header() -> String {
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

    pub(crate) fn to_row(
        &self,
        read_id: &str,
        chrom_name: Option<&String>,
        caller: &MultipleThresholdModCaller,
        reference_seqs: &HashMap<String, Vec<u8>>,
        flag: u16,
        pass_only: bool,
        skip_inferred: bool,
    ) -> Option<String> {
        let filtered = caller.call(&self.canonical_base, &self.base_mod_probs)
            == BaseModCall::Filtered;
        let inferred = self.base_mod_probs.inferred;
        if filtered && pass_only {
            return None;
        }
        if inferred && skip_inferred {
            return None;
        }

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
        let within_alignment = chrom_name.is_some() && self.within_alignment();

        Some(format!(
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
        ))
    }
}

pub trait OutwriterWithMemory<T> {
    fn write(&mut self, item: T, kmer_size: usize) -> anyhow::Result<u64>;
    fn num_reads(&self) -> usize;
}

pub struct TsvWriterWithContigNames<W: Write, C> {
    tsv_writer: TsvWriter<W>,
    tid_to_name: HashMap<u32, String>,
    name_to_seq: HashMap<String, Vec<u8>>,
    number_of_written_reads: usize,
    caller: C,
    pass_only: bool,
}

impl<W: Write> TsvWriterWithContigNames<W, ()> {
    pub(crate) fn new(
        output_writer: TsvWriter<W>,
        tid_to_name: HashMap<u32, String>,
        name_to_seq: HashMap<String, Vec<u8>>,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            tsv_writer: output_writer,
            tid_to_name,
            name_to_seq,
            number_of_written_reads: 0,
            caller: (),
            pass_only: false,
        })
    }
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W, ()>
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
            self.number_of_written_reads += 1;
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.number_of_written_reads
    }
}

impl<W: Write> TsvWriterWithContigNames<W, MultipleThresholdModCaller> {
    pub(crate) fn new_with_caller(
        output_writer: TsvWriter<W>,
        tid_to_name: HashMap<u32, String>,
        name_to_seq: HashMap<String, Vec<u8>>,
        caller: MultipleThresholdModCaller,
        pass_only: bool,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            tsv_writer: output_writer,
            tid_to_name,
            name_to_seq,
            number_of_written_reads: 0,
            caller,
            pass_only,
        })
    }
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W, MultipleThresholdModCaller>
{
    fn write(
        &mut self,
        item: ReadsBaseModProfile,
        _kmer_size: usize,
    ) -> anyhow::Result<u64> {
        let mut rows_written = 0u64;
        for profile in item.profiles.iter() {
            let chrom_name = profile
                .chrom_id
                .and_then(|chrom_id| self.tid_to_name.get(&chrom_id));
            let position_calls = PositionModCalls::from_profile(&profile);
            for call in position_calls {
                call.to_row(
                    &profile.record_name,
                    chrom_name,
                    &self.caller,
                    &self.name_to_seq,
                    profile.flag,
                    self.pass_only,
                    false,
                )
                .map(|s| self.tsv_writer.write(s.as_bytes()))
                .transpose()?;
                rows_written += 1;
            }
            self.number_of_written_reads += 1;
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.number_of_written_reads
    }
}
