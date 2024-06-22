use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter};
use std::path::Path;
use std::str;
use std::string::FromUtf8Error;

use anyhow::Result as AnyhowResult;
use anyhow::{anyhow, bail};
use bio::alphabets::dna::complement;
use derive_new::new;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use linear_map::LinearMap;
use log::{debug, error, info};
use regex::Regex;
use rust_htslib::bam::{
    self, ext::BamRecordExtensions, header::HeaderRecord, record::Aux,
    HeaderView, Read,
};

use crate::errs::{InputError, RunError};

pub(crate) const TAB: char = '\t';

pub(crate) fn create_out_directory<T: AsRef<std::ffi::OsStr>>(
    raw_path: T,
) -> anyhow::Result<()> {
    if let Some(p) = Path::new(&raw_path).parent() {
        if !p.exists() && p != Path::new("") {
            info!("creating directory at {p:?}");
            std::fs::create_dir_all(p)?;
        }
    }
    Ok(())
}

pub(crate) fn get_ticker() -> ProgressBar {
    let ticker = ProgressBar::new_spinner();
    ticker.set_style(ProgressStyle::with_template("> {pos} {msg}").unwrap());
    ticker
}
pub(crate) fn get_spinner() -> ProgressBar {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} [{elapsed_precise}] {pos} {msg}",
        )
        .unwrap()
        .tick_strings(&[
            "▹▹▹▹▹",
            "▸▹▹▹▹",
            "▹▸▹▹▹",
            "▹▹▸▹▹",
            "▹▹▹▸▹",
            "▹▹▹▹▸",
            "▪▪▪▪▪",
        ]),
    );
    spinner
}

fn get_master_progress_bar_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.green/yellow} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-")
}

fn get_subroutine_progress_bar_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.blue/cyan} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-")
}

fn get_gauge_style() -> ProgressStyle {
    ProgressStyle::with_template("{bar:40.red/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("||-")
}

pub(crate) fn get_master_progress_bar(n: usize) -> ProgressBar {
    ProgressBar::new(n as u64).with_style(get_master_progress_bar_style())
}

pub(crate) fn get_subroutine_progress_bar(n: usize) -> ProgressBar {
    ProgressBar::new(n as u64).with_style(get_subroutine_progress_bar_style())
}

pub(crate) fn get_guage(n_max: usize) -> ProgressBar {
    ProgressBar::new(n_max as u64).with_style(get_gauge_style())
}

pub(crate) fn get_aligned_pairs_forward(
    record: &bam::Record,
) -> impl Iterator<Item = AnyhowResult<(usize, u64)>> + '_ {
    let read_length = record.seq_len();
    record.aligned_pairs().map(move |pair| {
        let q_pos = pair[0] as usize;
        let q_pos = if record.is_reverse() {
            read_length.checked_sub(q_pos).and_then(|x| x.checked_sub(1))
        } else {
            Some(q_pos)
        };
        if q_pos.is_none() || pair[1] < 0 {
            let read_id = get_query_name_string(&record)
                .unwrap_or("failed-to-parse-utf8".to_owned());
            debug!("record {read_id} has invalid aligned pair {:?}", pair);
            return Err(anyhow!("pair {:?} is invalid", pair));
        }

        let r_pos = pair[1];
        assert!(r_pos >= 0);
        let r_pos = r_pos as u64;
        Ok((q_pos.unwrap(), r_pos))
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
) -> Result<(T, &'static str), RunError> {
    let tag_new = record.aux(tag_keys[0].as_bytes());
    let tag_old = record.aux(tag_keys[1].as_bytes());

    let (tag, t) = match (tag_new, tag_old) {
        (Ok(aux), _) => Ok((aux, tag_keys[0])),
        (Err(_), Ok(aux)) => Ok((aux, tag_keys[1])),
        _ => Err(RunError::new_skipped("AUX data not found")),
    }?;
    parser(&tag, t).map(|v| (v, t))
}

pub(crate) fn parse_nm(record: &bam::Record) -> anyhow::Result<u32> {
    let nm_tag = record.aux("NM".as_bytes())?;
    match nm_tag {
        Aux::U8(x) => Ok(x as u32),
        Aux::U16(x) => Ok(x as u32),
        Aux::U32(x) => Ok(x),
        Aux::I8(x) => Ok(x as u32),
        Aux::I16(x) => Ok(x as u32),
        Aux::I32(x) => Ok(x as u32),
        _ => bail!("invalid NM tag {nm_tag:?}"),
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash, Default, PartialOrd, Ord)]
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

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum StrandRule {
    Positive,
    Negative,
    Both,
}

impl StrandRule {
    pub fn covers(&self, strand: Strand) -> bool {
        match &self {
            StrandRule::Positive => strand == Strand::Positive,
            StrandRule::Negative => strand == Strand::Negative,
            StrandRule::Both => true,
        }
    }

    pub fn same_as(&self, strand: Strand) -> bool {
        match &self {
            StrandRule::Positive => strand == Strand::Positive,
            StrandRule::Negative => strand == Strand::Negative,
            StrandRule::Both => false,
        }
    }

    pub fn absorb(self, strand: Strand) -> Self {
        if self.same_as(strand) {
            self
        } else {
            // self is either both or they are opposite strands
            // so that means to "absorb" the rule is now both
            StrandRule::Both
        }
    }

    pub fn combine(self, other: Self) -> Self {
        if self == other {
            self
        } else {
            Self::Both
        }
    }
}

impl From<Strand> for StrandRule {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Positive => Self::Positive,
            Strand::Negative => Self::Negative,
        }
    }
}

impl TryFrom<char> for StrandRule {
    type Error = anyhow::Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            '+' => Ok(Self::Positive),
            '-' => Ok(Self::Negative),
            '.' => Ok(Self::Both),
            _ => bail!("illegal strand rule {value}"),
        }
    }
}

impl Display for StrandRule {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let lab = match self {
            Self::Positive => '+',
            Self::Negative => '-',
            Self::Both => '.',
        };

        write!(f, "{}", lab)
    }
}

pub fn record_is_not_primary(record: &bam::Record) -> bool {
    record.is_supplementary() || record.is_secondary() || record.is_duplicate()
}

pub(crate) fn get_targets(
    header: &HeaderView,
    region: Option<&Region>,
) -> Vec<ReferenceRecord> {
    (0..header.target_count())
        .filter_map(|tid| {
            let chrom_name = String::from_utf8(header.tid2name(tid).to_vec())
                .unwrap_or("???".to_owned());
            if let Some(region) = &region {
                if chrom_name == region.name {
                    Some(ReferenceRecord::new(
                        tid,
                        region.start,
                        region.length(),
                        chrom_name,
                    ))
                } else {
                    None
                }
            } else {
                match header.target_len(tid) {
                    Some(size) => Some(ReferenceRecord::new(
                        tid,
                        0,
                        size as u32,
                        chrom_name,
                    )),
                    None => {
                        debug!(
                            "no size information for {chrom_name} (tid: {tid})"
                        );
                        None
                    }
                }
            }
        })
        .collect::<Vec<ReferenceRecord>>()
}

#[derive(Debug, new)]
pub struct ReferenceRecord {
    // todo make this usize and unify all of the "Genome types"
    pub tid: u32,
    pub start: u32,
    pub(crate) length: u32,
    pub name: String,
}

impl ReferenceRecord {
    pub fn end(&self) -> u32 {
        self.start + self.length
    }
}

#[derive(new, Debug, Eq, PartialEq)]
pub struct Region {
    pub name: String,
    pub start: u32,
    pub end: u32,
}

impl Region {
    pub fn length(&self) -> u32 {
        self.end - self.start
    }

    fn parse_raw_with_start_and_end(raw: &str) -> Result<Self, InputError> {
        let mut splitted = raw.split(':');
        let chrom_name = splitted
            .nth(0)
            .ok_or(InputError::new(&format!("failed to parse region {raw}")))?;
        let start_end = splitted.collect::<Vec<&str>>();
        if start_end.len() != 1 {
            return Err(InputError::new(&format!(
                "failed to parse region {raw}"
            )));
        } else {
            let start_end = start_end[0];
            let splitted = start_end
                .split('-')
                .map(|x| {
                    let cleaned = x.replace(",", "");
                    cleaned
                        .parse::<u32>()
                        .map_err(|e| InputError::new(&e.to_string()))
                })
                .collect::<Result<Vec<u32>, _>>()?;
            if splitted.len() != 2 {
                return Err(InputError::new(&format!(
                    "failed to parse region {raw}"
                )));
            } else {
                let start = splitted[0];
                let end = splitted[1];
                if end <= start {
                    return Err(InputError::new(&format!(
                        "failed to parse region {raw}, end must be after start"
                    )));
                }
                Ok(Self { name: chrom_name.to_owned(), start, end })
            }
        }
    }

    pub fn parse_str(
        raw: &str,
        header: &HeaderView,
    ) -> Result<Self, InputError> {
        if raw.contains(':') {
            Self::parse_raw_with_start_and_end(raw)
        } else {
            let target_id = (0..header.target_count()).find_map(|tid| {
                String::from_utf8(header.tid2name(tid).to_vec()).ok().and_then(
                    |contig| if &contig == raw { Some(tid) } else { None },
                )
            });
            let target_length =
                target_id.and_then(|tid| header.target_len(tid));
            if let Some(len) = target_length {
                Ok(Self { name: raw.to_owned(), start: 0, end: len as u32 })
            } else {
                Err(InputError::new(&format!(
                    "failed to find matching reference sequence for {raw} in \
                     BAM header"
                )))
            }
        }
    }

    pub fn get_fetch_definition(
        &self,
        header: &HeaderView,
    ) -> AnyhowResult<bam::FetchDefinition> {
        let tid = (0..header.target_count())
            .find_map(|tid| {
                String::from_utf8(header.tid2name(tid).to_vec()).ok().and_then(
                    |chrom| {
                        if &chrom == &self.name {
                            Some(tid)
                        } else {
                            None
                        }
                    },
                )
            })
            .ok_or(anyhow!(
                "failed to find target ID for chrom {}",
                self.name.as_str()
            ))?;
        let tid = tid as i32;
        Ok(bam::FetchDefinition::Region(
            tid,
            self.start as i64,
            self.end as i64,
        ))
    }

    pub(crate) fn to_string(&self) -> String {
        format!("{}:{}-{}", self.name, self.start, self.end)
    }
}

// shouldn't need this once it's fixed in rust-htslib or the repo moves to
// noodles..
fn header_to_hashmap(
    header: &bam::Header,
) -> anyhow::Result<HashMap<String, Vec<LinearMap<String, String>>>> {
    let mut header_map = HashMap::default();
    let record_type_regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
    let tag_regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]*)").unwrap();

    if let Ok(header_string) = String::from_utf8(header.to_bytes()) {
        for line in header_string.split('\n').filter(|x| !x.is_empty()) {
            let parts: Vec<_> =
                line.split('\t').filter(|x| !x.is_empty()).collect();
            if parts.is_empty() {
                continue;
            }
            let record_type = record_type_regex
                .captures(parts[0])
                .and_then(|captures| captures.get(1))
                .map(|m| m.as_str().to_owned());

            if let Some(record_type) = record_type {
                if record_type == "CO" {
                    continue;
                }
                let mut field = LinearMap::default();
                for part in parts.iter().skip(1) {
                    if let Some(cap) = tag_regex.captures(part) {
                        let tag = cap.get(1).unwrap().as_str().to_owned();
                        let value = cap.get(2).unwrap().as_str().to_owned();
                        field.insert(tag, value);
                    } else {
                        debug!("encounted illegal record line {line}");
                    }
                }
                header_map
                    .entry(record_type)
                    .or_insert_with(Vec::new)
                    .push(field);
            } else {
                debug!("encountered illegal record type in line {line}");
            }
        }
        Ok(header_map)
    } else {
        bail!("failed to parse header string")
    }
}

pub fn add_modkit_pg_records(header: &mut bam::Header) {
    let header_map = match header_to_hashmap(&header) {
        Ok(hm) => hm,
        Err(_) => {
            error!(
                "failed to parse input BAM header, not adding PG header \
                 record for modkit"
            );
            return;
        }
    };
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

#[derive(new, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone)]
pub struct SamTag {
    inner: [u8; 2],
}

#[cfg(test)]
impl SamTag {
    pub(crate) fn parse(chars: [char; 2]) -> Self {
        Self { inner: [chars[0] as u8, chars[1] as u8] }
    }
}

pub(crate) fn get_stringable_aux(
    record: &bam::Record,
    sam_tag: &SamTag,
) -> Option<String> {
    record.aux(&sam_tag.inner).ok().and_then(|aux| match aux {
        Aux::String(s) => Some(s.to_string()),
        Aux::Char(c) => Some(format!("{}", c)),
        Aux::Double(f) => Some(format!("{}", f)),
        Aux::Float(f) => Some(format!("{}", f)),
        Aux::HexByteArray(a) => Some(a.to_string()),
        Aux::I8(i) => Some(format!("{}", i)),
        Aux::I16(i) => Some(format!("{}", i)),
        Aux::I32(i) => Some(format!("{}", i)),
        Aux::U8(u) => Some(format!("{}", u)),
        Aux::U16(u) => Some(format!("{}", u)),
        Aux::U32(u) => Some(format!("{}", u)),
        _ => None,
    })
}

pub(crate) fn parse_partition_tags(
    raw_tags: &[String],
) -> anyhow::Result<Vec<SamTag>> {
    let mut tags_seen = HashSet::with_capacity(raw_tags.len());
    let mut tags = Vec::with_capacity(raw_tags.len());
    for raw_tag in raw_tags {
        if raw_tag.len() != 2 {
            bail!("illegal tag {raw_tag} should be length 2")
        }
        let raw_tag_parts = raw_tag.chars().collect::<Vec<char>>();
        assert_eq!(raw_tag_parts.len(), 2);
        let inner = [raw_tag_parts[0] as u8, raw_tag_parts[1] as u8];
        let tag = SamTag::new(inner);

        let inserted = tags_seen.insert(tag);
        if inserted {
            tags.push(tag);
        } else {
            bail!("cannot repeat partition-tags, got {raw_tag} twice")
        }
    }

    Ok(tags)
}

#[inline]
pub fn get_reference_mod_strand(
    read_mod_strand: Strand,
    alignment_strand: Strand,
) -> Strand {
    match (read_mod_strand, alignment_strand) {
        (Strand::Positive, Strand::Positive) => Strand::Positive,
        (Strand::Positive, Strand::Negative) => Strand::Negative,
        (Strand::Negative, Strand::Positive) => Strand::Negative,
        (Strand::Negative, Strand::Negative) => Strand::Positive,
    }
}

#[inline]
pub(crate) fn reader_is_bam(reader: &bam::IndexedReader) -> bool {
    unsafe {
        (*reader.htsfile()).format.format
            == rust_htslib::htslib::htsExactFormat_bam
    }
}

pub(crate) const KMER_SIZE: usize = 50;

#[derive(Copy, Clone)]
pub(crate) struct Kmer {
    inner: [u8; KMER_SIZE],
    pub(crate) size: usize,
}

impl Kmer {
    pub(crate) fn from_seq(seq: &[u8], pos: usize, kmer_size: usize) -> Kmer {
        Kmer::new(seq, pos, kmer_size)
    }

    // kinda risky, size needs to be < 12
    pub(crate) fn new(seq: &[u8], position: usize, size: usize) -> Self {
        if size > KMER_SIZE {
            debug!("kmers greater that size 12 will be corrupted");
        }
        let get_back_base_safe = |i| -> Option<u8> {
            position.checked_sub(i).and_then(|idx| seq.get(idx).map(|b| *b))
        };
        let before = if size % 2 == 0 { size / 2 - 1 } else { size / 2 };

        let after = size / 2;
        let mut buffer = [Some(45u8); KMER_SIZE];
        let mut i = 0;
        let mut assign = |b: Option<u8>| {
            buffer[i] = b;
            i = std::cmp::min(i + 1, KMER_SIZE - 1);
        };

        for offset in (1..=before).rev() {
            let b = get_back_base_safe(offset);
            assign(b);
        }
        assign(seq.get(position).map(|b| *b));
        for offset in 1..=after {
            assign(seq.get(position + offset).map(|b| *b))
        }
        let inner = buffer.map(|b| b.unwrap_or(45));
        Self { inner, size }
    }

    pub(crate) fn reverse_complement(self) -> Self {
        let mut inner = [45u8; KMER_SIZE];
        for (i, p) in (0..self.size).rev().enumerate() {
            let mut b = self.inner[p];
            if b != 45 {
                b = complement(b)
            }
            inner[i] = b
        }
        Self { inner, size: self.size }
    }
}

impl Debug for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = self.inner.iter().take(self.size).map(|b| *b as char).join("");
        write!(f, "{s}")
    }
}

impl Display for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[inline]
pub fn within_alignment(
    query_position: usize,
    num_soft_clipped_start: usize,
    num_soft_clipped_end: usize,
    read_length: usize,
) -> bool {
    read_length
        .checked_sub(num_soft_clipped_end)
        .map(|x| query_position >= num_soft_clipped_start && query_position < x)
        .unwrap_or_else(|| {
            debug!(
                "read_length ({read_length}) is less than \
                 num_soft_clipped_end ({num_soft_clipped_end})"
            );
            false
        })
}

pub fn format_int_with_commas(val: isize) -> String {
    let mut num = val
        .abs()
        .to_string()
        .as_bytes()
        .rchunks(3)
        .rev()
        .map(str::from_utf8)
        .collect::<Result<Vec<&str>, _>>()
        .unwrap()
        .join(",");
    if val < 0 {
        num = format!("-{num}")
    }
    num
}

#[cfg(test)]
mod utils_tests {
    use anyhow::Context;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;

    use crate::util::{
        get_query_name_string, get_stringable_aux, parse_partition_tags, SamTag,
    };

    #[test]
    fn test_util_get_stringable_tag() {
        let bam_fp = "tests/resources/bc_anchored_10_reads.sorted.bam";
        let mut reader = bam::Reader::from_path(bam_fp).unwrap();
        let mut checked = false;
        for record in reader.records().filter_map(|r| r.ok()) {
            let read_id = get_query_name_string(&record).unwrap();
            if read_id == "068ce426-129e-4870-bd34-16cd78edaa43".to_string() {
                let tag = SamTag {
                    inner: [102, 98], // 'fb' not a tag that's there
                };
                assert!(get_stringable_aux(&record, &tag).is_none());
                let expected_rg =
                    "5598049b1b3264566b162bf035344e7ec610d608_dna_r10.4.1_e8.\
                     2_400bps_hac@v3.5.2"
                        .to_string();
                let tag = SamTag { inner: [82, 71] };
                assert_eq!(
                    get_stringable_aux(&record, &tag),
                    Some(expected_rg)
                );
                let tag = SamTag { inner: [114, 110] };
                assert_eq!(
                    get_stringable_aux(&record, &tag),
                    Some("6335".to_string())
                );
                checked = true
            }
        }
        assert!(checked)
    }

    #[test]
    fn test_util_parse_partition_tags() {
        let raw_tags = ["HP".to_string(), "RG".to_string()];
        let parsed = parse_partition_tags(&raw_tags)
            .context("should have parsed raw tags")
            .unwrap();
        let expected =
            vec![SamTag::parse(['H', 'P']), SamTag::parse(['R', 'G'])];
        assert_eq!(parsed, expected);
    }
}
