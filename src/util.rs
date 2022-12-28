use rust_htslib::bam::{self, ext::BamRecordExtensions, pileup::Pileup, record::Aux, Read};
use std::collections::VecDeque;

use crate::errs::{InputError, RunError};

struct ReferenceNameLookup {}

impl ReferenceNameLookup {
    pub(crate) fn get_reference_name(&self, tid: u32) -> String {
        unimplemented!()
    }
}

pub(crate) fn get_aligned_pairs_forward(
    record: &bam::Record,
) -> impl Iterator<Item = (usize, u64)> + '_ {
    let read_length = record.seq_len();
    record.aligned_pairs().map(move |pair| {
        assert_eq!(pair.len(), 2);
        let q_pos = pair[0] as usize;
        let q_pos = if record.is_reverse() {
            read_length
                .checked_sub(q_pos)
                .and_then(|x| x.checked_sub(1))
        } else {
            Some(q_pos)
        };
        assert!(q_pos.is_some(), "pair {:?} is invalid", pair);

        let r_pos = pair[1];
        assert!(r_pos >= 0);
        let r_pos = r_pos as u64;
        (q_pos.unwrap(), r_pos)
    })
}

#[inline]
pub(crate) fn get_forward_sequence(record: &bam::Record) -> Result<String, RunError> {
    let raw_seq = if record.is_reverse() {
        bio::alphabets::dna::revcomp(record.seq().as_bytes())
    } else {
        record.seq().as_bytes()
    };
    let seq = String::from_utf8(raw_seq).map_err(|e| {
        RunError::new_input_error(format!("failed to convert sequence to string, {}", e))
    })?;
    if seq.len() == 0 {
        return Err(RunError::new_failed("seq is empty"));
    }
    Ok(seq)
}

pub(crate) fn get_tag<T>(
    record: &bam::Record,
    tag_keys: &[&str; 2],
    parser: &dyn Fn(&Aux, &str) -> Result<T, RunError>,
) -> Option<Result<T, RunError>> {
    let tag_new = record.aux(tag_keys[0].as_bytes());
    let tag_old = record.aux(tag_keys[1].as_bytes());

    let tag = match (tag_new, tag_old) {
        (Ok(aux), _) => Some((aux, tag_keys[0])),
        (Err(_), Ok(aux)) => Some((aux, tag_keys[1])),
        _ => None,
    };

    tag.map(|(aux, t)| parser(&aux, t))
}

#[inline]
pub(crate) fn encode_seq(seq: &str) -> Vec<u8> {
    seq.chars()
        .flat_map(|c| match c {
            'A' => [1, 0, 0, 0],
            'C' => [0, 1, 0, 0],
            'G' => [0, 0, 1, 0],
            'T' => [0, 0, 0, 1],
            'H' => [1, 1, 0, 1],
            // todo, the rest of the IUPAC codes...
            _ => panic!("non-canonical base {c}"),
        })
        .collect()
}

#[inline]
fn array_mul(a: &[u8], b: &[u8]) -> u8 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| (*x) * (*y)).sum()
}

#[inline]
pub(crate) fn find_matching_positions<
    'a,
    const WINDOW_SIZE: usize,
    const STEP_SIZE: usize,
    const PATTERN_SIZE: u8,
>(
    encoded_seq: &'a [u8],
    encoded_pattern: &'a [u8],
) -> impl Iterator<Item = usize> + 'a {
    encoded_seq
        .windows(WINDOW_SIZE)
        .step_by(STEP_SIZE)
        .enumerate()
        .filter_map(|(i, chunk)| {
            let match_count = array_mul(chunk, encoded_pattern);
            if match_count == PATTERN_SIZE {
                Some(i)
            } else {
                None
            }
        })
}

#[inline]
pub(crate) fn dna_complement(base: char) -> Option<char> {
    match base {
        'A' => Some('T'),
        'C' => Some('G'),
        'G' => Some('C'),
        'T' => Some('A'),
        _ => None,
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash)]
pub enum Strand {
    Positive,
    Negative,
}

impl Strand {
    pub fn parse_char(x: char) -> Result<Self, InputError> {
        match x {
            '+' => Ok(Self::Positive),
            '-' => Ok(Self::Negative),
            _ => Err(format!("failed to parse strand {}", x).into()),
        }
    }
}
