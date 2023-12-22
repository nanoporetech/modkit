use crate::motif_bed::{MultipleMotifLocations, RegexMotif};

pub struct IntervalChunks {
    seq_len: u32,
    chunk_size: u32,
    curr: u32,
}

impl IntervalChunks {
    pub fn new(start: u32, seq_len: u32, chunk_size: u32) -> Self {
        Self { seq_len: start + seq_len, chunk_size, curr: start }
    }
}

impl Iterator for IntervalChunks {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.seq_len {
            None
        } else {
            let start = self.curr;
            let end = std::cmp::min(start + self.chunk_size, self.seq_len);
            self.curr = end;
            Some((start, end))
        }
    }
}

pub struct MotifAwareIntervalChunks<'a> {
    seq_len: u32,
    chunk_size: u32,
    curr: u32,
    reference_id: u32,
    motifs: Option<&'a MultipleMotifLocations>,
}

impl<'a> MotifAwareIntervalChunks<'a> {
    /// Make an interval iterator that will not make
    /// split points within motifs when the entire
    /// motif must be computed at once (as with
    /// `strand_combine`). E.g. if we have CpG
    /// motifs _and_ want to combine counts on the +
    /// and - strand calls, we need to make sure not to
    /// have one interval go up to the first C and leave
    /// the second C on the next interval. However, if
    /// we're not combining strands, they can be split
    /// since the counts are independent.
    pub fn new(
        start: u32,
        seq_len: u32,
        chunk_size: u32,
        reference_id: u32,
        strand_combine: bool,
        motif_locations: Option<&'a MultipleMotifLocations>,
    ) -> Self {
        Self {
            seq_len: start + seq_len,
            chunk_size,
            curr: start,
            motifs: if strand_combine { motif_locations } else { None },
            reference_id,
        }
    }
}

impl Iterator for MotifAwareIntervalChunks<'_> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.seq_len {
            None
        } else {
            let start = self.curr;
            let mut end = std::cmp::min(start + self.chunk_size, self.seq_len);
            if let Some(mls) = self.motifs {
                'creep_loop: loop {
                    let motifs_at_position = mls
                        .motifs_at_position(
                            self.reference_id,
                            end.checked_sub(1).unwrap_or(end),
                        )
                        .into_iter()
                        .filter(|motif| motif.length > 1)
                        .collect::<Vec<&RegexMotif>>();
                    if motifs_at_position.is_empty() {
                        break 'creep_loop;
                    }
                    let creep_out = motifs_at_position
                        .iter()
                        .map(|motif| {
                            motif
                                .length
                                .checked_sub(motif.forward_offset)
                                .unwrap_or(motif.length)
                        })
                        .max();
                    if let Some(creep) = creep_out {
                        end += creep as u32;
                        if creep == 0 {
                            // just in case so that we don't get stuck
                            // forever
                            break 'creep_loop;
                        } else {
                            continue 'creep_loop;
                        }
                    } else {
                        // shouldn't ever really hit this since we have
                        // a break when motifs_at_position is empty
                        break 'creep_loop;
                    }
                }
            }
            let end = std::cmp::min(end, self.seq_len);
            self.curr = end;
            Some((start, end))
        }
    }
}

pub fn slice_dna_sequence(str_seq: &str, start: usize, end: usize) -> String {
    str_seq
        .char_indices()
        .filter_map(
            |(pos, nt)| {
                if pos >= start && pos <= end {
                    Some(nt)
                } else {
                    None
                }
            },
        )
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
