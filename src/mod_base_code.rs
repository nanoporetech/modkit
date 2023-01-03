use itertools::Itertools;
use std::collections::VecDeque;

use crate::util;
use rust_htslib::bgzf::CompressionLevel::Default;

pub trait ModBaseCode {
    fn num_mods(&self) -> usize;
    fn idx_for_mod_code(&self, raw_code: char) -> Option<usize>;
    fn raw_mod_codes(&self) -> &'static [char];
}

#[derive(Copy, Clone)]
pub struct CytosineModBaseCode;
const CYTOSINE_MODBASE_CODES: [char; 5] = ['C', 'c', 'f', 'h', 'm'];
impl ModBaseCode for CytosineModBaseCode {
    fn num_mods(&self) -> usize {
        5
    }

    fn idx_for_mod_code(&self, raw_code: char) -> Option<usize> {
        match raw_code {
            'C' => Some(0),
            'c' => Some(1),
            'f' => Some(2),
            'h' => Some(3),
            'm' => Some(4),
            _ => None,
        }
    }

    fn raw_mod_codes(&self) -> &'static [char] {
        &CYTOSINE_MODBASE_CODES
    }
}

#[derive(Copy, Clone)]
pub struct HydroxyMethylCytosineCode;
const HMC_MODBASE_CODE_SUBSET: [char; 2] = ['h', 'm'];
impl ModBaseCode for HydroxyMethylCytosineCode {
    fn num_mods(&self) -> usize {
        2
    }

    fn idx_for_mod_code(&self, raw_code: char) -> Option<usize> {
        match raw_code {
            'h' => Some(0),
            'm' => Some(1),
            _ => None,
        }
    }

    fn raw_mod_codes(&self) -> &'static [char] {
        &HMC_MODBASE_CODE_SUBSET
    }
}

pub trait ModificationMotif {
    // todo? would maybe be better API to have this return an iterator
    fn find_matches(&self, seq: &str) -> VecDeque<(usize, util::Strand)>;
    fn canonical_base(&self) -> char;
    fn len(&self) -> u32;
    fn required_overlap(&self) -> u32 {
        let base = self.len() / 2;
        base + (self.len() % 2)
    }
}

struct CustomModificationMotif {
    motif_seq: String,
    canonical_base: char,
}

pub struct CpGModificationMotif;
impl ModificationMotif for CpGModificationMotif {
    fn find_matches(&self, seq: &str) -> VecDeque<(usize, util::Strand)> {
        seq.match_indices("CG")
            .flat_map(|(pos, _)| {
                [
                    (pos, util::Strand::Positive),
                    (pos + 1, util::Strand::Negative),
                ]
            })
            .collect()
    }

    fn canonical_base(&self) -> char {
        'C'
    }
    fn len(&self) -> u32 {
        2
    }
}

pub struct CHHModificationMotif;
const CHH_MOTIF: [u8; 12] = [0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1];
impl ModificationMotif for CHHModificationMotif {
    fn find_matches(&self, seq: &str) -> VecDeque<(usize, util::Strand)> {
        let encoded_seq = util::encode_seq(seq);
        let fw_pos = util::find_matching_positions::<12, 4, 3>(&encoded_seq, &CHH_MOTIF)
            .map(|pos| (pos, util::Strand::Positive));
        let rc_seq = bio::alphabets::dna::revcomp(seq.as_bytes());
        let rc_seq = String::from_utf8(rc_seq).unwrap();
        let encoded_rc_seq = util::encode_seq(&rc_seq);
        let rv_pos = util::find_matching_positions::<12, 4, 3>(&encoded_rc_seq, &CHH_MOTIF)
            .map(|pos| (seq.len() - 1 - pos, util::Strand::Negative));
        let matches = fw_pos
            .chain(rv_pos)
            .sorted_by(|(x, _), (y, _)| x.cmp(y))
            .collect::<VecDeque<(usize, util::Strand)>>();
        matches
    }
    fn canonical_base(&self) -> char {
        'C'
    }
    fn len(&self) -> u32 {
        3
    }
}

#[cfg(test)]
mod mod_base_code_tests {
    use crate::interval_chunks::IntervalChunks;
    use lazy_regex::regex;
    use std::collections::VecDeque;

    use crate::mod_base_code::{CHHModificationMotif, CpGModificationMotif, ModificationMotif};
    use crate::util::Strand;

    #[test]
    fn test_cpg_modification_motif() {
        //               01     78 0
        //                +-    +- +-
        let seq = "CCGGTGACGGCGG";
        let motif = CpGModificationMotif;
        let matches = motif.find_matches(&seq);
        let expected = vec![
            (1, Strand::Positive),
            (2, Strand::Negative),
            (7, Strand::Positive),
            (8, Strand::Negative),
            (10, Strand::Positive),
            (11, Strand::Negative),
        ];
        assert_eq!(matches.into_iter().collect::<Vec<_>>(), expected);
    }

    #[test]
    fn test_chh_modification_motifs() {
        //                ++  -
        let seq = "ACCTAG";
        let motif = CHHModificationMotif;
        let matches = motif.find_matches(&seq);
        let expected = vec![
            (1, Strand::Positive),
            (2, Strand::Positive),
            (5, Strand::Negative),
        ];
        assert_eq!(matches.into_iter().collect::<Vec<_>>(), expected);
    }

    #[test]
    fn test_splitting_and_overlap() {
        let motif = CpGModificationMotif;
        let seq = "ACGCGTTT".chars().collect::<Vec<_>>();
        let chunk_size = 4;
        let mut intervals =
            IntervalChunks::new(seq.len() as u32, chunk_size, motif.required_overlap());
        let (start, end) = intervals.next().unwrap();
        let (start, end) = (start as usize, end as usize);
        let sub_seq = &seq[start..end].iter().collect::<String>();
        assert_eq!(sub_seq, "ACGC");
        let motif_positions = motif.find_matches(sub_seq).into_iter().collect::<Vec<_>>();
        assert_eq!(
            motif_positions,
            vec![(1, Strand::Positive), (2, Strand::Negative)]
        );

        let (start, end) = intervals.next().unwrap();
        let (start, end) = (start as usize, end as usize);
        let sub_seq = &seq[start..end].iter().collect::<String>();
        assert_eq!(sub_seq, "CGTT");
        let motif_positions = motif.find_matches(sub_seq).into_iter().collect::<Vec<_>>();
        assert_eq!(
            motif_positions,
            vec![(0, Strand::Positive), (1, Strand::Negative)]
        );

        let motif = CHHModificationMotif;
        assert_eq!(motif.required_overlap(), 2);
        let seq = "AACCATT".chars().collect::<Vec<_>>();
        let mut intervals =
            IntervalChunks::new(seq.len() as u32, chunk_size, motif.required_overlap());
        let (start, end) = intervals.next().unwrap();
        let (start, end) = (start as usize, end as usize);
        let sub_seq = &seq[start..end].iter().collect::<String>();
        assert_eq!(sub_seq, "AACC");
        let motif_positions = motif.find_matches(sub_seq).into_iter().collect::<Vec<_>>();
        assert_eq!(motif_positions, vec![]);

        let (start, end) = intervals.next().unwrap();
        let (start, end) = (start as usize, end as usize);
        let sub_seq = &seq[start..end].iter().collect::<String>();
        assert_eq!(sub_seq, "CCAT");
        let motif_positions = motif.find_matches(sub_seq).into_iter().collect::<Vec<_>>();
        assert_eq!(
            motif_positions,
            vec![(0, Strand::Positive), (1, Strand::Positive)]
        );
    }

    #[test]
    #[ignore = "not yet implemented"]
    fn test_chg_modification_motifs() {
        todo!("implement CHG motif and test")
    }

    #[test]
    #[ignore = "not yet implemented"]
    fn test_generative_against_regex() {
        todo!("write test that generates random DNA and tests against regex engine")
    }

    #[test]
    fn test_generic_motif() {}
}
