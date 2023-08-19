use rustc_hash::FxHashMap;

use crate::motif_bed::MultipleMotifLocations;

pub struct IntervalChunks {
    seq_len: u32,
    chunk_size: u32,
    curr: u32,
    motif_positions: FxHashMap<u32, u32>,
}

impl IntervalChunks {
    pub fn new_without_motifs(
        start: u32,
        seq_len: u32,
        chunk_size: u32,
        reference_id: u32,
    ) -> Self {
        Self::new(start, seq_len, chunk_size, reference_id, None)
    }

    pub fn new_with_multiple_motifs(
        start: u32,
        seq_len: u32,
        chunk_size: u32,
        reference_id: u32,
        motif_locations: Option<&MultipleMotifLocations>,
    ) -> Self {
        Self::new(start, seq_len, chunk_size, reference_id, motif_locations)
    }

    fn new(
        start: u32,
        seq_len: u32,
        chunk_size: u32,
        reference_id: u32,
        motif_locations: Option<&MultipleMotifLocations>,
    ) -> Self {
        let motif_positions = motif_locations
            .map(|locations| {
                locations
                    .motif_locations
                    .iter()
                    .map(|loc| {
                        loc.get_locations_unchecked(reference_id)
                            .keys()
                            .copied()
                            .map(|pos| (pos, loc.motif_length() as u32))
                            .collect::<Vec<(u32, u32)>>()
                    })
                    .fold(
                        FxHashMap::<u32, u32>::default(),
                        |mut acc, motif_positions| {
                            for (position, length) in motif_positions {
                                if let Some(l) = acc.get_mut(&position) {
                                    *l = std::cmp::max(*l, length);
                                } else {
                                    acc.insert(position, length);
                                }
                            }
                            acc
                        },
                    )
            })
            .unwrap_or_else(|| FxHashMap::default());

        Self {
            seq_len: start + seq_len,
            chunk_size,
            curr: start,
            motif_positions,
        }
    }
}

impl Iterator for IntervalChunks {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.seq_len {
            None
        } else {
            let start = self.curr;
            let mut end = std::cmp::min(start + self.chunk_size, self.seq_len);
            while let Some(&len) = self.motif_positions.get(&(end - 1)) {
                end += len;
            }
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
    use rust_htslib::faidx;

    use crate::interval_chunks::slice_dna_sequence;
    use crate::test_utils::load_test_sequence;

    #[test]
    fn test_check_sequence_slicing_is_same_as_fetch() {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
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

        // let seq = "ABCDEF".chars().collect::<Vec<char>>();
        // let mut ic = IntervalChunks::new(0, seq.len() as u32, 3);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 0);
        // assert_eq!(e, 3);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['A', 'B', 'C']);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 2);
        // assert_eq!(e, 5);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['C', 'D', 'E']);
        // let (s, e) = ic.next().unwrap();
        // assert_eq!(s, 4);
        // assert_eq!(e, 6);
        // let (s, e) = (s as usize, e as usize);
        // assert_eq!(&seq[s..e], ['E', 'F']);
    }
}
