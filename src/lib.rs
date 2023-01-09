pub mod commands;
pub mod errs;
pub mod interval_chunks;
pub mod mod_bam;
pub mod mod_pileup;
pub mod writers;
// pub mod mod_base_code;
// pub mod mod_pileup;
mod read_cache;
mod util;

#[cfg(test)]
mod test_utils {
    use rust_htslib::faidx;

    pub(crate) fn load_test_sequence(name: &str) -> String {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let dna = fasta_reader.fetch_seq_string(name, 0, 156).unwrap();
        dna
    }
}
