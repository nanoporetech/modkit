use crate::extract::subcommand::PositionModCalls;
use crate::read_ids_to_base_mod_probs::ReadsBaseModProfile;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::writers::TsvWriter;
use derive_new::new;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

pub trait OutwriterWithMemory<T> {
    fn write(&mut self, item: T, kmer_size: usize) -> anyhow::Result<u64>;
    fn num_reads(&self) -> usize;
}

#[derive(new)]
pub struct TsvWriterWithContigNames<W: Write> {
    tsv_writer: TsvWriter<W>,
    tid_to_name: HashMap<u32, String>,
    name_to_seq: HashMap<String, Vec<u8>>,
    written_reads: HashSet<String>,
    read_calls_writer: Option<TsvWriter<File>>,
    read_pos_pileup_writer: Option<TsvWriter<File>>,
    caller: MultipleThresholdModCaller,
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
            if self.written_reads.contains(&profile.record_name) {
                continue;
            } else {
                let chrom_name = if let Some(chrom_id) = profile.chrom_id {
                    self.tid_to_name.get(&chrom_id)
                } else {
                    None
                };
                for mod_profile in profile.profile.iter() {
                    let row = mod_profile.to_row(
                        &profile.record_name,
                        chrom_name.unwrap_or(&missing_chrom),
                        &self.name_to_seq,
                        kmer_size,
                    );
                    self.tsv_writer.write(row.as_bytes())?;
                    rows_written += 1;
                }
                self.written_reads.insert(profile.record_name.to_owned());
                if let Some(read_calls_writer) = self.read_calls_writer.as_mut()
                {
                    let position_calls = PositionModCalls::from_profile(
                        profile.record_name.as_str(),
                        &profile.profile,
                    );
                    for call in position_calls {
                        read_calls_writer.write(
                            call.to_row(
                                &profile.record_name,
                                chrom_name,
                                &self.caller,
                                &self.name_to_seq,
                            )
                            .as_bytes(),
                        )?;
                    }
                }
            }
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.written_reads.len()
    }
}
