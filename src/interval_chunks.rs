use rayon::iter::plumbing::UnindexedConsumer;
use rayon::iter::ParallelIterator;

pub struct IntervalChunks {
    seq_len: u32,
    chunk_size: u32,
    overlap: u32,
    curr: u32,
}

impl IntervalChunks {
    pub fn new(seq_len: u32, chunk_size: u32, overlap: u32) -> Self {
        Self {
            seq_len,
            chunk_size,
            overlap,
            curr: 0,
        }
    }
}

impl Iterator for IntervalChunks {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.seq_len {
            None
        } else {
            let start = self.curr.checked_sub(self.overlap).unwrap_or(0);
            let end = std::cmp::min(start + self.chunk_size, self.seq_len);
            self.curr = end;
            Some((start, end))
        }
    }
}

pub fn slice_dna_sequence(str_seq: &str, start: usize, end: usize) -> String {
    str_seq
        .char_indices()
        .filter_map(|(pos, nt)| {
            if pos >= start && pos <= end {
                Some(nt)
            } else {
                None
            }
        })
        .collect::<String>()
}

#[cfg(test)]
mod interval_chunks_tests {
    use crate::interval_chunks::{slice_dna_sequence, IntervalChunks};
    use crate::test_utils::load_test_sequence;
    use rust_htslib::faidx;

    #[test]
    fn test_check_sequence_slicing_is_same_as_fetch() {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let name = "oligo_1512_adapters";
        let dna = load_test_sequence(name);
        let start = 49;
        let end = 99;
        let slice_a = slice_dna_sequence(&dna, start, end);
        let slice_b = fasta_reader.fetch_seq_string(name, start, end).unwrap();
        assert_eq!(slice_a, slice_b);
    }

    #[test]
    fn test_interval_chunks() {
        let seq = "ABCDEF".chars().collect::<Vec<char>>();
        let mut ic = IntervalChunks::new((seq.len() as u32), 3, 1);
        let (s, e) = ic.next().unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 3);
        let (s, e) = (s as usize, e as usize);
        assert_eq!(&seq[s..e], ['A', 'B', 'C']);
        let (s, e) = ic.next().unwrap();
        assert_eq!(s, 2);
        assert_eq!(e, 5);
        let (s, e) = (s as usize, e as usize);
        assert_eq!(&seq[s..e], ['C', 'D', 'E']);
        let (s, e) = ic.next().unwrap();
        assert_eq!(s, 4);
        assert_eq!(e, 6);
        let (s, e) = (s as usize, e as usize);
        assert_eq!(&seq[s..e], ['E', 'F']);
    }
}
