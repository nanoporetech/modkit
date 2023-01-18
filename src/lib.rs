pub mod commands;
pub mod errs;
pub mod interval_chunks;
pub mod logging;
pub mod mod_bam;
pub mod mod_pileup;
pub mod motif_bed;
mod read_cache;
pub mod summarize;
pub mod thresholds;
mod util;
pub mod writers;

#[cfg(test)]
mod test_utils {
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
