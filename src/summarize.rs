use crate::errs::RunError;
use crate::mod_bam::BaseModCall;
use crate::mod_base_code::{DnaBase, ModCode};
use crate::thresholds::modbase_records;
use crate::util::{get_master_progress_bar, get_spinner, Strand};

use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug)]
pub struct ModSummary {
    pub reads_with_mod_calls: HashMap<DnaBase, u64>,
    pub mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub filtered_mod_calls: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub total_reads_used: usize,
}

pub fn summarize_modbam<T: AsRef<Path>>(
    bam_fp: T,
    threads: usize,
    threshold: f32,
    num_reads: Option<usize>,
) -> Result<ModSummary, RunError> {
    let mut reader = bam::Reader::from_path(bam_fp)
        .map_err(|e| RunError::new_input_error(e.to_string()))?;
    reader.set_threads(threads).map_err(|e| {
        RunError::new_failed(format!(
            "failed to set threads on reader, {}",
            e.to_string()
        ))
    })?;

    let record_iter = modbase_records(reader.records());
    let spinner = if let Some(n) = num_reads {
        get_master_progress_bar(n)
    } else {
        get_spinner()
    };

    spinner.set_message("records processed");

    let mut total_reads_used = 0;
    let mut reads_with_mod_calls = HashMap::new();
    let mut mod_call_counts = HashMap::new();
    let mut filtered_mod_calls = HashMap::new();
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
            let filtered_counts = filtered_mod_calls
                .entry(canonical_base)
                .or_insert(HashMap::new());
            for (_position, base_mod_probs) in
                seq_pos_mod_probs.pos_to_base_mod_probs
            {
                let count = match base_mod_probs.base_mod_call() {
                    BaseModCall::Canonical(p) => {
                        if p > threshold {
                            mod_counts
                                .entry(
                                    canonical_base
                                        .canonical_mod_code()
                                        .unwrap(),
                                )
                                .or_insert(0)
                        } else {
                            filtered_counts
                                .entry(
                                    canonical_base
                                        .canonical_mod_code()
                                        .unwrap(),
                                )
                                .or_insert(0)
                        }
                    }
                    BaseModCall::Modified(p, mod_code) => {
                        if p < threshold {
                            filtered_counts.entry(mod_code).or_insert(0)
                        } else {
                            mod_counts.entry(mod_code).or_insert(0)
                        }
                    }
                    BaseModCall::Filtered => {
                        unreachable!("should not encounter filtered calls")
                    }
                };
                *count += 1;
            }
        }
        total_reads_used = i;
        spinner.inc(1);
        let done = num_reads.map(|n| i >= n).unwrap_or(false);
        if done {
            break;
        }
    }
    spinner.finish_and_clear();

    Ok(ModSummary {
        reads_with_mod_calls,
        mod_call_counts,
        filtered_mod_calls,
        total_reads_used: total_reads_used + 1,
    })
}
