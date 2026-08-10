extern crate core;

pub mod adjust;
pub mod bedmethyl_util;
// pub mod commands;
pub mod command_utils;
pub mod dmr;
pub mod entropy;
pub mod errs;
pub mod extract;
pub mod interval_chunks;
pub mod localise;
// pub mod logging;
pub mod mod_bam;
pub mod mod_base_code;
pub mod modbam_util;
pub mod monoid;
pub mod motifs;
pub(crate) mod ordered_scheduler;
pub mod pileup;
pub mod position_filter;
pub mod read_ids_to_base_mod_probs;
pub mod reads_sampler;
pub mod record_processor;
pub mod repair_tags;
pub mod stats;
pub mod summarize;
pub mod threshold_mod_caller;
pub mod thresholds;
pub mod util;
pub mod validate;
pub mod writers;

pub(crate) mod genome_positions;
pub(crate) mod parsing_utils;
pub(crate) mod sample_probs;

mod fasta;
mod hmm;
mod read_cache;
mod tabix;

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
        let fasta_fp = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
        let fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let dna = fasta_reader.fetch_seq_string(name, 0, 156).unwrap();
        dna
    }
}
