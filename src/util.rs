use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::string::FromUtf8Error;

use crate::errs::{InputError, RunError};

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

pub(crate) fn get_query_name_string(
    record: &bam::Record,
) -> Result<String, FromUtf8Error> {
    String::from_utf8(record.qname().to_vec())
}

#[inline]
pub(crate) fn get_forward_sequence(
    record: &bam::Record,
) -> Result<String, RunError> {
    let raw_seq = if record.is_reverse() {
        bio::alphabets::dna::revcomp(record.seq().as_bytes())
    } else {
        record.seq().as_bytes()
    };
    let seq = String::from_utf8(raw_seq).map_err(|e| {
        RunError::new_input_error(format!(
            "failed to convert sequence to string, {}",
            e
        ))
    })?;
    if seq.len() == 0 {
        return Err(RunError::new_failed("seq is empty"));
    }
    Ok(seq)
}

pub(crate) fn get_tag<T>(
    record: &bam::Record,
    tag_keys: &[&'static str; 2],
    parser: &dyn Fn(&Aux, &str) -> Result<T, RunError>,
) -> Option<Result<(T, &'static str), RunError>> {
    let tag_new = record.aux(tag_keys[0].as_bytes());
    let tag_old = record.aux(tag_keys[1].as_bytes());

    let tag = match (tag_new, tag_old) {
        (Ok(aux), _) => Some((aux, tag_keys[0])),
        (Err(_), Ok(aux)) => Some((aux, tag_keys[1])),
        _ => None,
    };

    tag.map(|(aux, t)| {
        let parsed = parser(&aux, t);
        parsed.map(|res| (res, t))
    })
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
    pub fn to_char(&self) -> char {
        match self {
            Strand::Positive => '+',
            Strand::Negative => '-',
        }
    }

    pub fn opposite(&self) -> Self {
        match self {
            Strand::Positive => Strand::Negative,
            Strand::Negative =>Strand::Positive,
        }
    }
}

pub fn record_is_secondary(record: &bam::Record) -> bool {
    record.is_supplementary() || record.is_secondary() || record.is_duplicate()
}
