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

#[cfg(test)]
mod interval_chunks_tests {
    use crate::interval_chunks::IntervalChunks;

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
