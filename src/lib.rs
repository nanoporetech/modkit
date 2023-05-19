pub mod commands;
pub mod errs;
// pub mod filter_thresholds;
pub mod adjust;
pub mod interval_chunks;
pub mod logging;
pub mod mod_bam;
pub mod mod_base_code;
pub mod mod_pileup;
pub mod monoid;
pub mod motif_bed;
pub mod summarize;
pub mod threshold_mod_caller;
pub mod thresholds;
pub mod writers;

mod read_cache;
mod read_ids_to_base_mod_probs;
/// Module contains functions for parallel processing
/// of individual reads and aggregating the results.
mod reads_sampler;
mod record_processor;
mod util;

#[cfg(test)]
pub mod test_utils {
    use rust_htslib::faidx;

    pub(crate) fn dna_complement(base: char) -> Option<char> {
        match base {
            'A' => Some('T'),
            'C' => Some('G'),
            'G' => Some('C'),
            'T' => Some('A'),
            _ => None,
        }
    }

    pub(crate) fn load_test_sequence(name: &str) -> String {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let dna = fasta_reader.fetch_seq_string(name, 0, 156).unwrap();
        dna
    }
}
