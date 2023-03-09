use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::string::FromUtf8Error;

use crate::errs::{InputError, RunError};
use derive_new::new;

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

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash, Default)]
pub enum Strand {
    #[default]
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
            Strand::Negative => Strand::Positive,
        }
    }
}

pub fn record_is_secondary(record: &bam::Record) -> bool {
    record.is_supplementary() || record.is_secondary() || record.is_duplicate()
}

#[derive(Debug, new)]
pub struct ReferenceRecord {
    pub tid: u32,
    pub start: u32,
    pub length: u32,
    pub name: String,
}

#[derive(Debug)]
pub struct Region {
    pub name: String,
    pub start: u32,
    pub end: u32,
}

impl Region {
    pub fn length(&self) -> u32 {
        self.end - self.start
    }

    pub fn parse_str(raw: &str) -> Result<Self, InputError> {
        let mut splitted = raw.split(':');
        let chrom_name = splitted
            .nth(0)
            .ok_or(InputError::new("failed to parse region {raw}"))?;
        let start_end = splitted.collect::<Vec<&str>>();
        if start_end.len() != 1 {
            return Err(InputError::new("failed to parse region {raw}"));
        } else {
            let start_end = start_end[0];
            let splitted = start_end
                .split('-')
                .map(|x| {
                    x.parse::<u32>()
                        .map_err(|e| InputError::new(&e.to_string()))
                })
                .collect::<Result<Vec<u32>, _>>()?;
            if splitted.len() != 2 {
                return Err(InputError::new("failed to parse region {raw}"));
            } else {
                let start = splitted[0];
                let end = splitted[1];
                if end <= start {
                    return Err(InputError::new(
                        "failed to parse region {raw}, end must be after start",
                    ));
                }
                Ok(Self {
                    name: chrom_name.to_owned(),
                    start,
                    end,
                })
            }
        }
    }
}

pub fn add_modkit_pg_records(header: &mut bam::Header) {
    let header_map = header.to_hashmap();
    let (id, pp) = if let Some(pg_tags) = header_map.get("PG") {
        let modkit_invocations = pg_tags.iter().filter_map(|tags| {
            tags.get("ID").and_then(|v| {
                if v.contains("modkit") {
                    let last_run = v.split('.').nth(1).unwrap_or("0");
                    last_run.parse::<usize>().ok()
                } else {
                    None
                }
            })
        });
        if let Some(latest_run_number) = modkit_invocations.max() {
            let pp = if latest_run_number > 0 {
                Some(format!("modkit.{}", latest_run_number))
            } else {
                Some(format!("modkit"))
            };
            (format!("modkit.{}", latest_run_number + 1), pp)
        } else {
            (format!("modkit"), None)
        }
    } else {
        (format!("modkit"), None)
    };

    let command_line = std::env::args().collect::<Vec<String>>();
    let command_line = command_line.join(" ");
    let version = env!("CARGO_PKG_VERSION");
    let mut modkit_header_record = HeaderRecord::new("PG".as_bytes());
    modkit_header_record.push_tag("ID".as_bytes(), &id);
    modkit_header_record.push_tag("PN".as_bytes(), &"modkit".to_owned());
    modkit_header_record.push_tag("VN".as_bytes(), &version.to_owned());
    if let Some(pp) = pp {
        modkit_header_record.push_tag("PP".as_bytes(), &pp);
    }
    modkit_header_record.push_tag("CL".as_bytes(), &command_line);

    header.push_record(&modkit_header_record);
}
