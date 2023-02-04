use crate::errs::RunError;
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::util::{record_is_secondary, Strand};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Debug)]
pub struct ModSummary {
    pub mod_called_bases: Vec<DnaBase>,
    pub reads_with_mod_calls: HashMap<DnaBase, u64>,
    pub mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub total_reads_used: usize,
}

pub fn summarize_modbam<T: AsRef<Path>>(
    bam_fp: T,
    threads: usize,
) -> Result<ModSummary, RunError> {
    let mut reader = bam::Reader::from_path(bam_fp)
        .map_err(|e| RunError::new_input_error(e.to_string()))?;
    reader.set_threads(threads).map_err(|e| {
        RunError::new_failed(format!(
            "failed to set threads on reader, {}",
            e.to_string()
        ))
    })?;

    let record_iter = reader
        .records()
        // skip records that fail to parse htslib
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        .filter_map(|record| ModBaseInfo::new_from_record(&record).ok());

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} [{elapsed_precise}] {pos} {msg}",
        )
        .unwrap()
        .tick_strings(&[
            "▹▹▹▹▹",
            "▸▹▹▹▹",
            "▹▸▹▹▹",
            "▹▹▸▹▹",
            "▹▹▹▸▹",
            "▹▹▹▹▸",
            "▪▪▪▪▪",
        ]),
    );
    spinner.set_message("records processed");

    let mut total_reads_used = 0;
    let mod_called_bases = HashSet::new();
    let mut reads_with_mod_calls = HashMap::new();
    let mut mod_call_counts = HashMap::new();
    for (i, modbase_info) in record_iter.enumerate() {
        if modbase_info.is_empty() {
            continue;
        }

        let (_converters, prob_iter) = modbase_info.into_iter_base_mod_probs();
        for (canonical_base, strand, seq_pos_mod_probs) in prob_iter {
            let canonical_base = match (DnaBase::parse(canonical_base), strand)
            {
                (Err(_), _) => continue,
                (Ok(dna_base), Strand::Positive) => dna_base,
                (Ok(dna_base), Strand::Negative) => dna_base.complement(),
            };
            let count = reads_with_mod_calls.entry(canonical_base).or_insert(0);
            *count += 1;
            let mod_counts = mod_call_counts
                .entry(canonical_base)
                .or_insert(HashMap::new());
            for (_position, base_mod_probs) in
                seq_pos_mod_probs.pos_to_base_mod_probs
            {
                match base_mod_probs.base_mod_call() {
                    BaseModCall::Canonical(_p) => {
                        let count = mod_counts
                            .entry(canonical_base.canonical_mod_code().unwrap())
                            .or_insert(0);
                        *count += 1;
                    }
                    BaseModCall::Modified(_p, mod_code) => {
                        let count = mod_counts.entry(mod_code).or_insert(0);
                        *count += 1;
                    }
                    BaseModCall::Filtered => {}
                }
            }
        }
        total_reads_used = i;
        spinner.inc(1);
    }

    let mut mod_called_bases = mod_called_bases.into_iter().collect::<Vec<_>>();
    mod_called_bases.sort();

    Ok(ModSummary {
        mod_called_bases,
        reads_with_mod_calls,
        mod_call_counts,
        total_reads_used: total_reads_used + 1,
    })
}
