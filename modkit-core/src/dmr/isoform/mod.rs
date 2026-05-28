use std::{
    collections::{HashMap, HashSet},
    fmt::{Debug, Write},
    fs::File,
    io::{BufRead, BufWriter, Write as IoWrite},
    ops::Neg,
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{
    dmr::isoform::{
        gtf_reader::open_gtf_reader,
        pair::run_gene_dmr,
        scoring::{dirichlet_multinomial_lrt, DirichletMultinomialTestResult},
        vis::{
            bar_color, build_compressed_layout, draw_intron_break, escape_xml,
            format_u64_with_commas, map_genomic_boundary_to_x,
            merged_gene_exon_blocks, methylation_color,
            modification_tick_position, GenomeSortedTranscriptModel,
        },
    },
    position_filter::Iv,
};
use crate::{
    mod_base_code::MOD_CODE_TO_DNA_BASE,
    util::{get_human_readable_table, TAB},
};
use anyhow::{anyhow, bail, Context};
use bitvec::bitvec;
use clap::Args;
use common_macros::hash_map;
use derive_new::new;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use itertools::Itertools;
use log::{debug, info, warn};
use log_once::warn_once;
use prettytable::row;
use rust_htslib::tbx::{Read, Reader as TbxReader};
use rustc_hash::FxHashMap;

use crate::{
    dmr::bedmethyl::BedMethylLine,
    errs::{MkError, MkResult},
    mod_base_code::ModCodeRepr,
    util::{get_ticker, Strand},
};
mod gtf_reader;
pub(super) mod pair;
mod scoring;
mod vis;
pub(super) use pair::PairTranscriptomeBedmethylHandler;

pub(crate) type TranscripId = String;
// pub(crate) type TranscriptKey = (TranscripId, u32);
pub(super) type GeneId = String;

pub(crate) trait GtfId {
    fn from_str(raw: &str) -> MkResult<Self>
    where
        Self: Sized;
}

#[derive(new, Debug, Clone, Hash, Eq, PartialEq, PartialOrd)]
pub(crate) struct GtfTranscript {
    tx_id: TranscripId,
    version: u32,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, PartialOrd)]
pub(crate) struct GtfGene {
    gene_id: GeneId,
    version: u32,
}

fn parse_gtf_id(raw: &str) -> MkResult<(String, u32)> {
    if raw.contains(".") {
        let Some((id, v)) = raw.split_once(".") else {
            return Err(MkError::InvalidGtfRecord);
        };
        let version =
            v.parse::<u32>().map_err(|_| MkError::InvalidGtfRecord)?;
        Ok((id.to_string(), version))
    } else {
        Ok((raw.to_string(), 0u32))
    }
}

impl GtfId for GtfTranscript {
    fn from_str(raw: &str) -> MkResult<Self> {
        let (tx_id, version) = parse_gtf_id(raw)?;
        Ok(Self { tx_id, version })
    }
}

impl GtfId for GtfGene {
    fn from_str(raw: &str) -> MkResult<Self> {
        let (gene_id, version) = parse_gtf_id(raw)?;
        Ok(Self { gene_id, version })
    }
}

impl std::fmt::Display for GtfTranscript {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}", self.tx_id, self.version)
    }
}

impl std::fmt::Display for GtfGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}", self.gene_id, self.version)
    }
}

#[derive(Debug, Clone)]
struct Exon {
    start0: u64, // 0-based inclusive
    end0: u64,   // 0-based exclusive
}

#[derive(Debug, Clone)]
struct Segment {
    start0: u64, // 0-based inclusive
    end0: u64,   // 0-based exclusive
}

#[derive(Debug, Clone)]
pub(crate) struct TranscriptModel {
    pub transcript_id: GtfTranscript,
    gene_id: GtfGene,
    gene_name: Option<String>,
    pub chrom: String,
    strand: char,              // '+' or '-'
    exons_tx_order: Vec<Exon>, // ordered in transcript 5'->3'
    pub transcript_len: u64,   // total spliced transcript length
}

#[derive(Debug, Clone)]
pub(super) struct GeneCommonCoord {
    pub gene_id: GtfGene,
    pub gene_name: Option<String>,
    pub chrom: String,
    strand: char,
    segments_tx_order: Vec<Segment>, // union-of-exons, ordered 5'->3'
    cumulative0: Vec<u64>,           /* 0-based start of each segment in
                                      * common coord */
    total_len: u64,
}

impl GeneCommonCoord {
    pub(super) fn empty() -> Self {
        Self {
            gene_id: GtfGene { gene_id: "".to_string(), version: 0u32 },
            gene_name: None,
            chrom: "".to_string(),
            strand: '+',
            segments_tx_order: Vec::new(),
            cumulative0: Vec::new(),
            total_len: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct GeneModel {
    gene_id: GtfGene,
    gene_name: Option<String>,
    chrom: String,
    strand: char,
    transcript_ids: Vec<GtfTranscript>,
    start0: u64,
    end0: u64,
}

// pub(super) struct GeneCommonCoord<'a> {
//     pub gene_id: &'a GeneId,
//     pub gene_name: &'a Option<String>,
//     pub chrom: &'a str,
//     strand: Strand,
//     segments_tx_order: Vec<Segment>, /* union of exonic intervals in
//                                          * transcript order */
//     cumulative: Vec<u64>, // 0-based segment starts in common coord
//     total_len: u64,
// }

fn parse_gtf_attributes(attr: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in attr.split(';') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        let mut it = p.splitn(2, ' ');
        if let (Some(key), Some(val)) = (it.next(), it.next()) {
            let cleaned = val.trim().trim_matches('"').to_string();
            map.insert(key.to_string(), cleaned);
        }
    }
    map
}

pub(crate) fn parse_gtf<P: AsRef<Path>>(
    gtf_path: P,
    multi_progress: &MultiProgress,
) -> anyhow::Result<IndexMap<GtfTranscript, TranscriptModel>> {
    let reader = open_gtf_reader(gtf_path)?;

    #[derive(Debug)]
    struct TmpTranscript {
        gene_id: GtfGene,
        gene_name: Option<String>,
        chrom: String,
        strand: char,
        exons: Vec<Exon>,
    }

    let mut tmp: IndexMap<GtfTranscript, TmpTranscript> = IndexMap::new();

    let pb = multi_progress.add(get_ticker());
    pb.set_message("parsed GTF records");
    let errs = multi_progress.add(get_ticker());
    errs.set_message("errored GTF records");
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chrom = fields[0].to_string();
        let feature = fields[2];
        if feature != "exon" {
            continue;
        }

        // GTF is 1-based inclusive
        let start1: u64 = fields[3].parse()?;
        let end1: u64 = fields[4].parse()?;

        if start1 == 0 || end1 < start1 {
            bail!("Invalid GTF exon coordinates: {}-{}", start1, end1)
        }

        // Convert to 0-based half-open
        let start0 = start1 - 1;
        let end0 = end1;

        let strand = fields[6]
            .chars()
            .next()
            .ok_or(anyhow!("Missing strand field in GTF"))?;

        let attrs = parse_gtf_attributes(fields[8]);
        let transcript_id = attrs
            .get("transcript_id")
            .ok_or(anyhow!("Missing transcript_id in GTF exon record"))?;

        let transcript_id = if let Some(transcript_version) =
            attrs.get("transcript_version").and_then(|x| x.parse::<u32>().ok())
        {
            GtfTranscript {
                tx_id: (*transcript_id).clone(),
                version: transcript_version,
            }
        } else {
            let Ok(transcript_id) = GtfTranscript::from_str(transcript_id)
            else {
                debug!("failed to parse transcript_id: {transcript_id}");
                errs.inc(1);
                continue;
            };
            transcript_id
        };
        let gene_id = attrs
            .get("gene_id")
            .ok_or(anyhow!("Missing gene_id in GTF exon record"))?;

        let gene_id = if let Some(gene_version) =
            attrs.get("gene_version").and_then(|x| x.parse::<u32>().ok())
        {
            GtfGene { gene_id: (*gene_id).clone(), version: gene_version }
        } else {
            let Ok(gene_id) = GtfGene::from_str(gene_id) else {
                debug!("failed to parse gene_id: {gene_id}");
                errs.inc(1);
                continue;
            };
            gene_id
        };
        let gene_name = attrs.get("gene_name").map(|x| x.to_owned());

        let entry = tmp.entry(transcript_id).or_insert(TmpTranscript {
            gene_id,
            gene_name,
            chrom,
            strand,
            exons: Vec::new(),
        });

        entry.exons.push(Exon { start0, end0 });
        pb.inc(1);
    }

    let mut out = IndexMap::new();

    for (transcript_id, mut t) in tmp {
        // Keep exons in transcript 5'->3' order
        if t.strand == '+' {
            t.exons.sort_by_key(|e| e.start0);
        } else {
            t.exons.sort_by(|a, b| b.start0.cmp(&a.start0));
        }

        let transcript_len: u64 =
            t.exons.iter().map(|e| e.end0 - e.start0).sum();

        out.insert(
            transcript_id.clone(),
            TranscriptModel {
                transcript_id,
                gene_id: t.gene_id,
                gene_name: t.gene_name,
                chrom: t.chrom,
                strand: t.strand,
                exons_tx_order: t.exons,
                transcript_len,
            },
        );
    }

    if out.is_empty() {
        bail!("failed to parse any transcript models")
    }

    Ok(out)
}

fn merge_intervals0(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut sorted = intervals.to_vec();
    sorted.sort_by_key(|x| x.0);

    let mut merged = Vec::new();
    let mut cur = sorted[0];

    for &(s, e) in &sorted[1..] {
        // For half-open intervals, merge if overlapping or abutting.
        if s <= cur.1 {
            cur.1 = cur.1.max(e);
        } else {
            merged.push(cur);
            cur = (s, e);
        }
    }

    merged.push(cur);
    merged
}

/// Build a common gene coordinate system as the union of all exons across
/// isoforms. All coordinates are 0-based half-open.
pub(super) fn build_gene_common_coords(
    transcripts: &IndexMap<GtfTranscript, TranscriptModel>,
) -> anyhow::Result<IndexMap<GtfGene, GeneCommonCoord>> {
    #[derive(Debug)]
    struct GeneTmp {
        chrom: String,
        gene_name: Option<String>,
        strand: char,
        intervals: Vec<(u64, u64)>,
    }

    let mut by_gene: IndexMap<GtfGene, GeneTmp> = IndexMap::new();

    for tx in transcripts.values() {
        let entry = by_gene.entry(tx.gene_id.clone()).or_insert(GeneTmp {
            chrom: tx.chrom.clone(),
            gene_name: tx.gene_name.to_owned(),
            strand: tx.strand,
            intervals: Vec::new(),
        });

        if entry.chrom != tx.chrom {
            bail!(
                "Gene {} spans multiple chromosomes: {} and {}",
                tx.gene_id,
                entry.chrom,
                tx.chrom
            )
        }
        if entry.strand != tx.strand {
            bail!("Gene {} has inconsistent strand assignments", tx.gene_id)
        }

        for exon in &tx.exons_tx_order {
            entry.intervals.push((exon.start0, exon.end0));
        }
    }

    let mut out = IndexMap::new();

    for (gene_id, g) in by_gene {
        let merged = merge_intervals0(&g.intervals);

        let segments_tx_order: Vec<Segment> = if g.strand == '+' {
            merged
                .into_iter()
                .map(|(s, e)| Segment { start0: s, end0: e })
                .collect()
        } else {
            let mut segs: Vec<Segment> = merged
                .into_iter()
                .map(|(s, e)| Segment { start0: s, end0: e })
                .collect();
            segs.sort_by(|a, b| b.start0.cmp(&a.start0));
            segs
        };

        let mut cumulative0 = Vec::with_capacity(segments_tx_order.len());
        let mut running = 0_u64;
        for seg in &segments_tx_order {
            cumulative0.push(running);
            running += seg.end0 - seg.start0;
        }

        out.insert(
            gene_id.clone(),
            GeneCommonCoord {
                gene_id,
                gene_name: g.gene_name,
                chrom: g.chrom,
                strand: g.strand,
                segments_tx_order,
                cumulative0,
                total_len: running,
            },
        );
    }

    Ok(out)
}

/// Convert a 0-based transcript position to a 0-based genomic position.
/// tx_pos0 must be in [0, transcript_len).
pub(crate) fn transcript_pos0_to_genomic0(
    tx: &TranscriptModel,
    tx_pos0: u64,
) -> MkResult<u64> {
    if tx_pos0 >= tx.transcript_len {
        debug!(
            "Transcript position {} out of bounds for transcript {} (len={})",
            tx_pos0, tx.transcript_id, tx.transcript_len
        );
        return Err(MkError::InvalidTranscriptPosition);
    }

    let mut cum0 = 0_u64;
    for exon in &tx.exons_tx_order {
        let exon_len = exon.end0 - exon.start0;
        if tx_pos0 < cum0 + exon_len {
            let offset = tx_pos0 - cum0; // 0-based within exon in transcript order

            let gpos0 = if tx.strand == '+' {
                exon.start0 + offset
            } else {
                (exon.end0 - 1) - offset
            };

            return Ok(gpos0);
        }
        cum0 += exon_len;
    }

    Err(MkError::InvalidTranscriptPosition)
}

/// Convert a 0-based genomic position to a 0-based common gene coordinate.
/// gpos0 must lie within the gene's union-of-exons.
fn genomic0_to_gene_common0(
    gene: &GeneCommonCoord,
    gpos0: u64,
) -> MkResult<u64> {
    for (i, seg) in gene.segments_tx_order.iter().enumerate() {
        if gpos0 >= seg.start0 && gpos0 < seg.end0 {
            let seg_common_start0 = gene.cumulative0[i];
            let offset = if gene.strand == '+' {
                gpos0 - seg.start0
            } else {
                (seg.end0 - 1) - gpos0
            };
            return Ok(seg_common_start0 + offset);
        }
    }

    debug!(
        "0-based genomic position {} not found in union-of-exons for gene {}",
        gpos0, gene.gene_id
    );

    Err(MkError::InvalidGenomicPosition)
}

/// Convert a 0-based common gene coordinate back to a 0-based genomic position.
/// common_pos0 must be in [0, total_len).
// fn gene_common0_to_genomic0(
//     gene: &GeneCommonCoord,
//     common_pos0: u64,
// ) -> MkResult<u64> {
//     if common_pos0 >= gene.total_len {
//         debug!(
//             "Gene common coordinate {} out of bounds for gene {} (len={})",
//             common_pos0, gene.gene_id, gene.total_len
//         );
//         return Err(MkError::InvalidGenomicPosition);
//     }

//     for (i, seg) in gene.segments_tx_order.iter().enumerate() {
//         let seg_common_start0 = gene.cumulative0[i];
//         let seg_len = seg.end0 - seg.start0;
//         let seg_common_end0 = seg_common_start0 + seg_len;

//         if common_pos0 >= seg_common_start0 && common_pos0 < seg_common_end0
// {             let offset = common_pos0 - seg_common_start0;

//             let gpos0 = if gene.strand == '+' {
//                 seg.start0 + offset
//             } else {
//                 (seg.end0 - 1) - offset
//             };

//             return Ok(gpos0);
//         }
//     }

//     debug!(
//         "Gene common coordinate {} not found in union-of-exons for gene {}",
//         common_pos0, gene.gene_id
//     );
//     Err(MkError::InvalidGenomicPosition)
// }

// pub(super) fn transcript_pos_to_genomic(
//     tx: &TranscriptModel,
//     tx_pos: u64,
// ) -> MkResult<u64> {
//     if tx_pos == 0 || tx_pos > tx.transcript_len {
//         debug!(
//             "Transcript position {} out of bounds for transcript {}
// (len={})",             tx_pos, tx.transcript_id, tx.transcript_len
//         );
//         return Err(MkError::InvalidTranscriptPosition);
//     }

//     let mut cum = 0_u64;
//     for (start, end) in &tx.exons_tx_order {
//         assert!(end > start);
//         let exon_len = *end - *start + 1;
//         if tx_pos <= cum + exon_len {
//             let offset = tx_pos - cum - 1; // 0-based within exon in
// transcript order

//             let gpos0 = if tx.strand == '+' {
//                 // exon.start is 1-based inclusive -> convert to 0-based
//                 (*start - 1) + offset
//             } else {
//                 // exon.end is 1-based inclusive, so the last base in 0-based
// is                 // exon.end - 1
//                 (*end - 1) - offset
//             };

//             return Ok(gpos0);
//         }
//         cum += exon_len;
//     }

//     debug!("Failed to map transcript position to genome");
//     Err(MkError::InvalidTranscriptPosition)
// }

// pub(super) fn genomic_to_gene_common(
//     gene: &GeneCommonCoord,
//     gpos: u64,
// ) -> MkResult<u64> {
//     for (i, (start, end)) in gene.segments_tx_order.iter().enumerate() {
//         if gpos >= *start && gpos <= *end {
//             let seg_common_start = gene.cumulative[i];
//             let offset = if gene.strand == Strand::Positive {
//                 gpos - *start
//             } else {
//                 *end - gpos
//             };
//             return Ok(seg_common_start + offset);
//         }
//     }

//     debug!(
//         "Genomic position {} not found in union-of-exons for gene {}",
//         gpos, gene.gene_id
//     );
//     return Err(MkError::InvalidGenomicPosition);
// }

pub(super) fn build_gene_models(
    transcripts: &FxHashMap<GtfTranscript, TranscriptModel>,
) -> anyhow::Result<HashMap<GtfGene, GeneModel>> {
    let mut by_gene: HashMap<GtfGene, GeneModel> = HashMap::new();

    for tx in transcripts.values() {
        let tx_start = tx
            .exons_tx_order
            .iter()
            .sorted_by_key(|e| e.start0)
            .map(|e| e.start0)
            .min()
            .unwrap_or(0);
        let tx_end = tx
            .exons_tx_order
            .iter()
            .sorted_by_key(|e| e.start0)
            .map(|e| e.end0)
            .max()
            .unwrap_or(0);

        let entry = by_gene.entry(tx.gene_id.clone()).or_insert(GeneModel {
            gene_id: tx.gene_id.clone(),
            gene_name: tx.gene_name.clone(),
            chrom: tx.chrom.clone(),
            strand: tx.strand,
            transcript_ids: Vec::new(),
            start0: tx_start,
            end0: tx_end,
        });

        if entry.chrom != tx.chrom {
            bail!("Gene {} spans multiple chromosomes", tx.gene_id)
        }
        if entry.strand != tx.strand {
            bail!("Gene {} has inconsistent strand assignments", tx.gene_id)
        }

        entry.transcript_ids.push(tx.transcript_id.clone());
        entry.start0 = entry.start0.min(tx_start);
        entry.end0 = entry.end0.max(tx_end);
    }

    Ok(by_gene)
}

#[derive(new, Copy, Clone)]
pub(super) struct TidAndLength {
    tid: u64,
    length: u64,
}

pub(super) struct TranscriptBedmethylHandler {
    tbx_reader: TbxReader,
    pub gene_id_to_transcript_ids: FxHashMap<GtfGene, Vec<TidAndLength>>,
    min_valid_coverage: u64,
    multi_progress: MultiProgress,
}

impl TranscriptBedmethylHandler {
    pub(super) fn from_path<R: AsRef<Path> + Debug + Clone>(
        fp: R,
        transcript_models: &FxHashMap<GtfTranscript, TranscriptModel>,
        min_valid_coverage: u64,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let tbx_reader = TbxReader::from_path(fp.clone())?;
        let mut missing_in_gtf = 0usize;
        let mut total_tx = 0usize;
        let gene_id_to_transcript_ids = tbx_reader
            .seqnames()
            .into_iter()
            .map(|n| tbx_reader.tid(&n).map(|tid| (n, tid)))
            .map_ok(|(contig, tid)| {
                let transcript_id = if contig.contains('|') {
                    let (x, _) = contig.split_once("|").unwrap();
                    x.to_string()
                } else {
                    contig
                };
                (transcript_id, tid)
            })
            .filter_map_ok(|(raw_transcript_id, tid)| {
                GtfTranscript::from_str(&raw_transcript_id)
                    .ok()
                    .map(|tx_id| (tx_id, tid))
            })
            .filter_map_ok(|(transcript_id, tid)| {
                if let Some(tm) = transcript_models.get(&transcript_id) {
                    total_tx = total_tx.saturating_add(1);
                    Some((tm.gene_id.clone(), tid, tm.transcript_len))
                } else {
                    missing_in_gtf = missing_in_gtf.saturating_add(1);
                    // debug!("transcript_id: {} not in GTF", transcript_id);
                    None
                }
            })
            .fold_ok(FxHashMap::default(), |mut acc, (gene_id, tid, length)| {
                acc.entry(gene_id)
                    .or_insert_with(Vec::new)
                    .push(TidAndLength::new(tid, length));
                acc
            })
            .context(
                "failed to collect transcript-ids, invalid tabix header?",
            )?;
        if gene_id_to_transcript_ids.is_empty() {
            bail!("zero transcripts in bedMethyl header could be found in GTF")
        }
        multi_progress.suspend(|| {
            if missing_in_gtf > 0 {
                warn!(
                    "{fp:?} has transcripts in tbx header missing from GTF, \
                     {missing_in_gtf} transcripts in bedMethyl do not have \
                     transcript models in GTF and cannot be used"
                );
            }
            let check = gene_id_to_transcript_ids
                .values()
                .map(|x| x.len())
                .sum::<usize>();
            assert_eq!(check, total_tx);
            info!(
                "{total_tx} transcripts parsed from bedmethyl header with \
                 corresponding transcript models in GTF"
            );
        });

        Ok(Self {
            tbx_reader,
            min_valid_coverage,
            gene_id_to_transcript_ids,
            multi_progress,
        })
    }

    pub(super) fn from_path_with_lookup_table<P: AsRef<Path>>(
        fp: P,
        gene_id_to_transcript_ids: FxHashMap<GtfGene, Vec<TidAndLength>>,
        min_valid_coverage: u64,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let tbx_reader = TbxReader::from_path(fp)?;
        Ok(Self {
            tbx_reader,
            min_valid_coverage,
            gene_id_to_transcript_ids,
            multi_progress,
        })
    }

    pub(super) fn get_read_bedmethyl_for_gene(
        &mut self,
        gene_id: &GtfGene,
        single_mod_code: Option<ModCodeRepr>,
    ) -> MkResult<Option<Vec<BedMethylLine>>> {
        let Some(tids) = self.gene_id_to_transcript_ids.get(gene_id) else {
            return Ok(None);
        };
        let mut out = Vec::new();
        let mut errs = 0usize;
        for tid_and_length in tids {
            self.tbx_reader.fetch(
                tid_and_length.tid,
                0u64,
                tid_and_length.length,
            )?;
            for record in self.tbx_reader.records() {
                let Ok(bs) = record else {
                    errs = errs.saturating_add(1);
                    continue;
                };
                let Ok(raw) = String::from_utf8(bs) else {
                    errs = errs.saturating_add(1);
                    continue;
                };
                if let Ok(bml) = BedMethylLine::parse(&raw) {
                    if bml.valid_coverage >= self.min_valid_coverage {
                        match single_mod_code {
                            Some(c) if bml.raw_mod_code == c => out.push(bml),
                            Some(_) => {}
                            None => out.push(bml),
                        }
                    }
                } else {
                    errs = errs.saturating_add(1);
                }
            }
        }
        if errs > 0usize {
            self.multi_progress.suspend(|| {
                warn_once!(
                    "some bedmethyl records failed to be read, check logs"
                );
            });
            debug!("{gene_id} had {errs} failures");
        }

        Ok(Some(out))
    }
}

#[derive(Debug)]
pub(super) struct GeneIsoformDmrScore {
    score: f64,
    p_value: f64,
    proportions: Vec<f64>,
    max_absolute_difference: f64,
    mod_codes: Vec<ModCodeRepr>,
    per_tx_proportions: Vec<f64>,
    per_tx_counts: Vec<u32>,
    tx_names: Vec<String>,
}

impl GeneIsoformDmrScore {
    fn from_score_result(
        res: DirichletMultinomialTestResult,
        tx_ids: &[&GtfTranscript],
        mod_codes: Vec<ModCodeRepr>,
        emit_full_results: bool,
    ) -> Self {
        let score = res.lrt;
        let p_value = res.p_value;
        let proportions =
            scoring::proportions_from_counts(&res.pooled_counts, 0.0);
        let mut mad = None;
        for (i, x_proportions) in res
            .isoform_counts
            .iter()
            .map(|x| scoring::proportions_from_counts(x, 0f64))
            .enumerate()
        {
            assert!(!x_proportions.is_empty());
            'y: for (j, y_proportions) in res
                .isoform_counts
                .iter()
                .map(|x| scoring::proportions_from_counts(x, 0f64))
                .enumerate()
            {
                if i == j {
                    continue 'y;
                }
                assert!(!y_proportions.is_empty());
                assert_eq!(x_proportions.len(), y_proportions.len());
                let abs_diff = x_proportions
                    .iter()
                    .zip(y_proportions.iter())
                    .map(|(x, y)| (*x - *y).abs())
                    .max_by(|x, y| {
                        x.partial_cmp(y)
                            .expect("should not have NaN in proportions")
                    })
                    .unwrap();
                if let Some(x) = mad.as_mut() {
                    if *x < abs_diff {
                        *x = abs_diff;
                    }
                } else {
                    mad = Some(abs_diff);
                }
            }
        }
        let max_absolute_difference = mad.unwrap();

        if !emit_full_results {
            Self {
                score,
                p_value,
                proportions,
                max_absolute_difference,
                mod_codes,
                per_tx_proportions: Vec::new(),
                per_tx_counts: Vec::new(),
                tx_names: Vec::new(),
            }
        } else {
            Self {
                score,
                p_value,
                proportions,
                max_absolute_difference,
                mod_codes,
                per_tx_proportions: res
                    .isoform_counts
                    .iter()
                    .map(|x| scoring::proportions_from_counts(x, 0f64))
                    .flatten()
                    .collect::<Vec<f64>>(),
                per_tx_counts: res
                    .isoform_counts
                    .into_iter()
                    .flatten()
                    .collect::<Vec<u32>>(),
                tx_names: tx_ids.iter().map(|x| x.to_string()).collect(),
            }
        }
    }
}

fn score_position(
    isoforms: &Vec<Vec<u32>>,
    tx_ids: &[&GtfTranscript],
    mod_codes: Vec<ModCodeRepr>,
    emit_full_results: bool,
) -> Option<GeneIsoformDmrScore> {
    if mod_codes.iter().filter_map(|x| MOD_CODE_TO_DNA_BASE.get(x)).count() == 1
    {
        dirichlet_multinomial_lrt(isoforms, 0.05, 0.5).map(|res| {
            GeneIsoformDmrScore::from_score_result(
                res,
                tx_ids,
                mod_codes,
                emit_full_results,
            )
        })
    } else {
        None
    }
}

#[derive(Debug)]
pub(super) struct GeneIsoformDmrRecord<S: Debug> {
    pub start: u64,
    end: u64,
    strand: Strand,
    n_transcripts: usize,
    score: S,
}

impl GeneIsoformDmrRecord<GeneIsoformDmrScore> {
    pub(super) fn header(full_results: bool) -> String {
        let mut h = vec![
            "#chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "p_value",
            "max_absolute_difference",
            "n_transcripts",
            "gene_id",
            "gene_name",
        ];
        if full_results {
            h.extend_from_slice(&[
                "pooled_proportions",
                "per_isoform_proportions",
                "per_isoform_counts",
            ]);
        }
        let mut h = h.join("\t");
        h.push('\n');
        h
    }
    pub(super) fn to_row(
        &self,
        chrom: &str,
        gene_id: &GtfGene,
        gene_name: Option<&String>,
        full_results: bool,
    ) -> String {
        let mut row = format!(
                        "\
                        {}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}",
                        chrom,
                        self.start,
                        self.end,
                        self.score.mod_codes.iter().join(","),
                        self.score.score,
                        self.strand.to_char(),
                        self.score.p_value,
                        self.score.max_absolute_difference,
                        self.n_transcripts,
                        gene_id,
                        gene_name.unwrap_or(&"-".to_string())
                    );
        if full_results {
            assert_eq!(
                self.score.mod_codes.len().saturating_add(1),
                self.score.proportions.len(),
                "offending {:?}, {row}",
                self.score
            );
            let mut proportions = hash_map! {
                "unmodified".to_string() => self.score.proportions[0],
            };
            for (i, code) in self.score.mod_codes.iter().enumerate() {
                let idx = i.saturating_add(1);
                proportions
                    .insert(format!("{code}"), self.score.proportions[idx]);
            }
            let proportions_ser =
                serde_json::to_string(&proportions).expect("should serialize");

            let mut per_tx_proportions = HashMap::new();
            let stride = self.score.mod_codes.len().saturating_add(1);
            for (i, tx) in self.score.tx_names.iter().enumerate() {
                let start = i * stride;
                let end = start.saturating_add(stride);
                let tx_proportions_sl =
                    &self.score.per_tx_proportions[start..end];
                assert_eq!(tx_proportions_sl.len(), stride);
                let mut tx_props = hash_map! {
                    "unmodified".to_string() => tx_proportions_sl[0],
                };
                for (i, code) in self.score.mod_codes.iter().enumerate() {
                    let idx = i.saturating_add(1);
                    tx_props.insert(format!("{code}"), tx_proportions_sl[idx]);
                }
                per_tx_proportions.insert(tx, tx_props);
            }
            let per_tx_props_ser = serde_json::to_string(&per_tx_proportions)
                .expect("should serialize per-tx proportions");

            let mut per_tx_counts = HashMap::new();
            for (i, tx) in self.score.tx_names.iter().enumerate() {
                let start = i * stride;
                let end = start.saturating_add(stride);
                let tx_counts_sl = &self.score.per_tx_counts[start..end];
                assert_eq!(tx_counts_sl.len(), stride);
                let mut tx_counts = hash_map! {
                    "unmodified".to_string() => tx_counts_sl[0],
                };
                for (i, code) in self.score.mod_codes.iter().enumerate() {
                    let idx = i.saturating_add(1);
                    tx_counts.insert(format!("{code}"), tx_counts_sl[idx]);
                }
                per_tx_counts.insert(tx, tx_counts);
            }
            let per_tx_counts_ser =
                serde_json::to_string(&per_tx_counts).unwrap();

            row.push_str(&format!(
                "{TAB}{proportions_ser}{TAB}{per_tx_props_ser}{TAB}{per_tx_counts_ser}"
            ));
        };
        row.push('\n');
        row
    }

    pub(super) fn score(&self, pvalue: bool) -> f64 {
        if pvalue {
            self.score.p_value
        } else {
            self.score.score
        }
    }

    pub(super) fn has_code(&self, mod_code: ModCodeRepr) -> bool {
        self.score.mod_codes.contains(&mod_code)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ModCodeCount {
    Empty,
    Count(ModCodeRepr, u32),
}
pub(super) struct TranscriptCounts<const NMODS: usize> {
    modified_counts: Vec<[ModCodeCount; NMODS]>,
    unmodified_counts: Vec<u32>,
}

pub(super) fn run_isoform_dmr_on_gene<'a, const NMODS: usize>(
    transcripts: &'a FxHashMap<GtfTranscript, TranscriptModel>,
    gene: &'a GeneCommonCoord,
    handler: &mut TranscriptBedmethylHandler,
    single_mod_code: Option<ModCodeRepr>,
    emit_full_results: bool,
) -> MkResult<Option<Vec<GeneIsoformDmrRecord<GeneIsoformDmrScore>>>> {
    let Some(bedmethyl_records) =
        handler.get_read_bedmethyl_for_gene(&gene.gene_id, single_mod_code)?
    else {
        return Ok(None);
    };

    let grouped_by_tx = group_bedmethyl_records_by_tx(bedmethyl_records);

    let mut valid_pos = bitvec![0; gene.total_len as usize];
    let mut genomic_positions = vec![0u64; gene.total_len as usize];
    let mut transcript_counts = HashMap::with_capacity(grouped_by_tx.len());

    'transcripts: for (tx_id, bedmethyl_records) in grouped_by_tx {
        assert!(!bedmethyl_records.is_empty());
        let Some(transcript_model) = transcripts.get(&tx_id) else {
            debug!("no transcript model for {tx_id}");
            continue 'transcripts;
        };
        let mut unmodified_counts = vec![0u32; gene.total_len as usize];
        let mut modified_counts =
            vec![[ModCodeCount::Empty; NMODS]; gene.total_len as usize];

        // base case
        let mut tx_start = bedmethyl_records[0].start();
        let mut genomic_start =
            transcript_pos0_to_genomic0(transcript_model, tx_start)?;
        let mut common_start =
            genomic0_to_gene_common0(gene, genomic_start)? as usize;
        if valid_pos[common_start] {
            assert_eq!(genomic_positions[common_start], genomic_start);
        } else {
            valid_pos.set(common_start, true);
            genomic_positions[common_start] = genomic_start
        }

        unmodified_counts[common_start] =
            bedmethyl_records[0].count_canonical as u32;

        for bml in bedmethyl_records.iter() {
            if bml.start() != tx_start {
                // update
                tx_start = bml.start();
                genomic_start =
                    transcript_pos0_to_genomic0(transcript_model, tx_start)?;
                common_start =
                    genomic0_to_gene_common0(gene, genomic_start)? as usize;
                valid_pos.set(common_start, true);
                genomic_positions[common_start] = genomic_start;
                unmodified_counts[common_start] = bml.count_canonical as u32;
            }

            let mod_counts = &mut modified_counts[common_start];
            if single_mod_code.is_some() {
                match mod_counts[0] {
                    ModCodeCount::Empty => {
                        let x = ModCodeCount::Count(
                            bml.raw_mod_code,
                            bml.count_methylated as u32,
                        );
                        mod_counts[0] = x;
                        unmodified_counts[common_start] = unmodified_counts
                            [common_start]
                            .saturating_add(bml.count_other as u32);
                    }
                    ModCodeCount::Count(_, _) => todo!(),
                }
            } else {
                'counts: for mod_count in mod_counts {
                    match mod_count {
                        ModCodeCount::Count(mod_code, count)
                            if *mod_code == bml.raw_mod_code =>
                        {
                            *count = count
                                .saturating_add(bml.count_methylated as u32);
                            break 'counts;
                        }
                        ModCodeCount::Count(_, _) => {
                            continue 'counts;
                        }
                        ModCodeCount::Empty => {
                            *mod_count = ModCodeCount::Count(
                                bml.raw_mod_code,
                                bml.count_methylated as u32,
                            );
                            break 'counts;
                        }
                    }
                }
            }
        }
        transcript_counts.insert(
            tx_id,
            TranscriptCounts { modified_counts, unmodified_counts },
        );
    }

    let gpos_iter =
        valid_pos.iter().zip(genomic_positions.iter()).enumerate().filter_map(
            |(i, (b, gpos))| if *b { Some((i, *gpos)) } else { None },
        );

    let mut out_records = Vec::new();
    for (i, gpos) in gpos_iter {
        let mut pos_mod_counts = Vec::new();
        let mut mod_codes_at_pos = Vec::new();
        let mut tx_ids = Vec::new();
        'tx: for (tx_id, counts) in transcript_counts.iter() {
            let modified_counts = &counts.modified_counts[i];
            if modified_counts[0] == ModCodeCount::Empty {
                continue 'tx;
            }
            let unmodified_count = counts.unmodified_counts[i];
            let mut tx_pos_mod_counts = vec![unmodified_count];
            'mod_counts: for mod_code_count in modified_counts {
                match mod_code_count {
                    ModCodeCount::Empty => break 'mod_counts,
                    ModCodeCount::Count(mod_code_repr, count) => {
                        mod_codes_at_pos.push(*mod_code_repr);
                        tx_pos_mod_counts.push(*count);
                    }
                }
            }
            tx_ids.push(tx_id);
            pos_mod_counts.push(tx_pos_mod_counts);
        }
        let mod_codes_at_pos =
            mod_codes_at_pos.into_iter().unique().collect::<Vec<ModCodeRepr>>();
        if let Some(gene_isoform_dmr_score) = score_position(
            &pos_mod_counts,
            &tx_ids,
            mod_codes_at_pos,
            emit_full_results,
        ) {
            out_records.push(GeneIsoformDmrRecord {
                start: gpos,
                end: gpos.saturating_add(1),
                strand: Strand::parse_char(gene.strand).unwrap(),
                n_transcripts: pos_mod_counts.len(),
                score: gene_isoform_dmr_score,
            });
        }
    }

    if gene.strand == '-' {
        Ok(Some(out_records.into_iter().rev().collect()))
    } else {
        Ok(Some(out_records))
    }
}

pub(super) struct GeneIsoformDmr {
    gene: GeneCommonCoord,
    records: Option<Vec<GeneIsoformDmrRecord<GeneIsoformDmrScore>>>,
}

impl GeneIsoformDmr {
    pub(super) fn empty() -> Self {
        Self { gene: GeneCommonCoord::empty(), records: None }
    }

    pub(super) fn clear(&mut self) {
        self.gene = GeneCommonCoord::empty();
        self.records = None;
    }

    pub(super) fn write(
        &self,
        writer: &mut Box<dyn IoWrite>,
        emit_full_results: bool,
    ) -> anyhow::Result<usize> {
        let mut records_written = 0usize;
        if let Some(records) = self.records.as_ref() {
            for record in records {
                let row = record.to_row(
                    &self.gene.chrom,
                    &self.gene.gene_id,
                    self.gene.gene_name.as_ref(),
                    emit_full_results,
                );
                writer.write(row.as_bytes())?;
                records_written = records_written.saturating_add(1);
            }
        }
        Ok(records_written)
    }
}

pub(super) struct GeneIsoformDmrWorker {
    transcript_bedmethyl_handler: TranscriptBedmethylHandler,
    transcript_models: Arc<FxHashMap<GtfTranscript, TranscriptModel>>,
    single_mod_code: Option<ModCodeRepr>,
    emit_full_results: bool,
}

impl GeneIsoformDmrWorker {
    pub(super) fn new(
        bedmethyl_fp: &PathBuf,
        transcript_models: Arc<FxHashMap<GtfTranscript, TranscriptModel>>,
        multi_progress: MultiProgress,
        gene_ids_to_transcript_ids: FxHashMap<GtfGene, Vec<TidAndLength>>,
        min_valid_coverage: u64,
        single_mod_code: Option<ModCodeRepr>,
        emit_full_results: bool,
    ) -> anyhow::Result<Self> {
        let transcript_bedmethyl_handler =
            TranscriptBedmethylHandler::from_path_with_lookup_table(
                bedmethyl_fp,
                gene_ids_to_transcript_ids,
                min_valid_coverage,
                multi_progress,
            )?;
        Ok(Self {
            transcript_bedmethyl_handler,
            transcript_models,
            single_mod_code,
            emit_full_results,
        })
    }
    pub(super) fn process_gene<const NMODS: usize>(
        &mut self,
        mut mem: GeneIsoformDmr,
        gene: GeneCommonCoord,
    ) -> MkResult<GeneIsoformDmr> {
        let records = run_isoform_dmr_on_gene::<NMODS>(
            &self.transcript_models,
            &gene,
            &mut self.transcript_bedmethyl_handler,
            self.single_mod_code,
            self.emit_full_results,
        )?;
        mem.gene = gene;
        mem.records = records;
        Ok(mem)
    }
}

#[derive(Args, Clone)]
pub(super) struct PlottingArgs {
    /// Comma-separated list of genomic positions to draw a vertical dotted
    /// line and a label.
    #[arg(
        long,
        requires = "plot",
        alias = "markers",
        num_args = 1..,
        value_delimiter = ','
    )]
    marker_positions: Option<Vec<u64>>,
    /// SVG width in px
    #[arg(
        long,
        default_value_t = 1600,
        requires = "plot",
        hide_short_help = true
    )]
    width: u32,

    /// Height per transcript track in px
    #[arg(
        long,
        default_value_t = 26,
        requires = "plot",
        hide_short_help = true
    )]
    track_height: u32,

    /// Gap between transcript tracks in px
    #[arg(
        long,
        default_value_t = 18,
        requires = "plot",
        hide_short_help = true
    )]
    track_gap: u32,

    /// Left margin in px
    #[arg(
        long,
        default_value_t = 220,
        requires = "plot",
        hide_short_help = true
    )]
    left_margin: u32,

    /// Right margin in px
    #[arg(
        long,
        default_value_t = 40,
        requires = "plot",
        hide_short_help = true
    )]
    right_margin: u32,

    /// Top margin in px
    #[arg(
        long,
        default_value_t = 210,
        requires = "plot",
        hide_short_help = true
    )]
    top_margin: u32,

    /// Bottom margin in px
    #[arg(
        long,
        default_value_t = 50,
        requires = "plot",
        hide_short_help = true
    )]
    bottom_margin: u32,

    /// Length, in px, of the compressed intron segments
    #[arg(
        long,
        default_value_t = 24.0,
        requires = "plot",
        hide_short_help = true
    )]
    intron_px: f64,

    /// Height in px of the top bar.
    #[arg(
        long,
        default_value_t = 90.0,
        requires = "plot",
        hide_short_help = true
    )]
    top_bar_height: f64,

    /// Only plot methylation markers with negative log p-value greater than
    /// this value.
    #[arg(long, requires = "plot", conflicts_with = "max_pval")]
    min_score: Option<f64>,

    /// Only plot methylation markers with p-value less than this value.
    #[arg(long, requires = "plot", conflicts_with = "min_score")]
    max_pval: Option<f64>,
}

pub(super) struct DmrMetric {
    pub start: u64,
    pub value: f64,
}

pub(super) fn plot_isoform_metrics(
    meth_by_tx: FxHashMap<GtfTranscript, Vec<BedMethylLine>>,
    dmr_metrics: &[DmrMetric],
    mod_code: ModCodeRepr,
    plot_dir: &PathBuf,
    gene_model: &GeneModel,
    gene_name: Option<&String>,
    all_transcripts: &FxHashMap<GtfTranscript, TranscriptModel>,
    plotting_args: &PlottingArgs,
    multi_progress: &MultiProgress,
) -> anyhow::Result<()> {
    let plot_fn = if let Some(gene_name) = gene_name {
        format!("{}_{gene_name}_{mod_code}.svg", &gene_model.gene_id)
    } else {
        format!("{}_{mod_code}.svg", &gene_model.gene_id)
    };
    let rel_ab_fn = if let Some(gene_name) = gene_name {
        format!("{}_{gene_name}.csv", &gene_model.gene_id)
    } else {
        format!("{}.csv", &gene_model.gene_id)
    };

    let tx_to_mean = meth_by_tx
        .iter()
        .map(|(tx_id, bm_records)| {
            let total = bm_records
                .iter()
                .map(|x| x.valid_coverage.saturating_add(x.count_fail))
                .sum::<u64>();
            (tx_id, total / bm_records.len() as u64)
        })
        .collect::<HashMap<_, _>>();

    let tot = tx_to_mean.values().sum::<u64>();
    let tx_to_relab = tx_to_mean
        .iter()
        .map(|(tx_id, mean)| ((*tx_id).clone(), *mean as f64 / tot as f64))
        .collect::<FxHashMap<GtfTranscript, f64>>();

    let plot_fp = plot_dir.join(plot_fn);
    let rel_ab_fp = plot_dir.join(rel_ab_fn);
    multi_progress.suspend(|| {
        info!("writing SVG plot to {plot_fp:?}");
        info!("writing relative abundances to {rel_ab_fp:?}");
    });

    let mut tab = get_human_readable_table();
    tab.set_titles(row![
        "transcript_id",
        "mean_total_coverage",
        "relative_abundance"
    ]);
    for (tx_id, mean) in tx_to_mean.iter() {
        tab.add_row(row![tx_id, mean, tx_to_relab.get(tx_id).unwrap()]);
    }
    let svg = render_gene_svg(
        gene_model,
        mod_code,
        &all_transcripts,
        meth_by_tx,
        tx_to_relab,
        dmr_metrics,
        plotting_args,
        multi_progress,
    )?;
    std::fs::write(plot_fp, svg)?;
    let tab_writer = BufWriter::new(File::create(rel_ab_fp)?);
    tab.to_csv(tab_writer)?;
    Ok(())
}

#[rustfmt::skip]
fn render_gene_svg(
    gene: &GeneModel,
    single_code: ModCodeRepr,
    all_transcripts: &FxHashMap<GtfTranscript, TranscriptModel>,
    meth_by_tx: FxHashMap<GtfTranscript, Vec<BedMethylLine>>,
    rel_abundances: FxHashMap<GtfTranscript, f64>,
    dmr_metrics: &[DmrMetric],
    args: &PlottingArgs,
    multi_progress: &MultiProgress,
) -> anyhow::Result<String> {
    let (min_score, has_score_thresh) = match (args.max_pval, args.min_score) {
        (None, None) => (0f64, false),
        (None, Some(min_score)) => {
            multi_progress.suspend(|| {
                info!("only plotting methylated positions with score >= {min_score:.4}");
            });
            (min_score, true)
        },
        (Some(max_pval), None) => {
            let min_score =  max_pval.ln().neg();
            multi_progress.suspend(|| {
                info!("only plotting methylated positions with score >= {min_score:.4} (max p-value {max_pval})");
            });
            (min_score, true)
        },
        (Some(_), Some(_)) => unreachable!(),
    };

    let meth_by_tx = meth_by_tx
        .into_iter()
        .map(|(tx_id, mut bmls)| {
            let tm = all_transcripts.get(&tx_id).unwrap();
            for bml in bmls.iter_mut() {
                let start = bml.start();
                // let stop = bml.stop();
                let genome_start = transcript_pos0_to_genomic0(tm, start).unwrap();
                let genome_stop = genome_start;
                // let genome_stop = transcript_pos0_to_genomic0(tm, stop).unwrap();
                bml.interval = Iv {start: genome_start, stop: genome_stop,  val: ()};
            }
            (tx_id, bmls)
        }).collect::<FxHashMap<GtfTranscript, Vec<BedMethylLine>>>();
    let mut tx_ids = gene.transcript_ids.clone()
        .into_iter()
        .filter_map(|x| {
            if meth_by_tx.contains_key(&x) {
                Some(x)
            } else {
                debug!("zero bedmethyl records for transcript-id: {x}");
                None
            }
        }).collect::<Vec<GtfTranscript>>();
    let transcripts = tx_ids
        .iter()
        .map(|tx_id| {
            (tx_id.to_owned(), all_transcripts.get(tx_id).unwrap().clone())
        })
        .map(|(tx_id, mut tm)| {
            tm.exons_tx_order.sort_by_key(|e| e.start0);
            let tm = GenomeSortedTranscriptModel::new(
                tm.strand,
                tm .exons_tx_order,
            );
            (tx_id, tm)
        })
        .collect::<HashMap<GtfTranscript, GenomeSortedTranscriptModel>>();

    tx_ids.sort_by(|a, b| {
        let ta = transcripts.get(a).unwrap();
        let tb = transcripts.get(b).unwrap();

        let a_key = ta.exons.first().map(|e| e.start0).unwrap_or(0);
        let b_key = tb.exons.first().map(|e| e.start0).unwrap_or(0);

        a_key.cmp(&b_key).then_with(|| a.to_string().cmp(&b.to_string()))
    });

    let inner_width = args
        .width
        .saturating_sub(args.left_margin + args.right_margin) as f64;
    let n_tracks = tx_ids.len() as u32;
    let height = args.top_margin
        + args.bottom_margin
        + n_tracks * args.track_height
        + n_tracks.saturating_sub(1) * args.track_gap
        + 40;

    let x0 = args.left_margin as f64;
    let x1 = (args.width - args.right_margin) as f64;
    let axis_y = args.top_margin as f64 - 25.0;

    let exon_blocks = merged_gene_exon_blocks(gene, &transcripts);
    let layout = build_compressed_layout(&exon_blocks, x0, inner_width, args.intron_px);

    let mut svg = String::new();
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"#,
        w = args.width,
        h = height
    )?;
    writeln!(svg, r#"<rect width="100%" height="100%" fill="white"/>"#)?;
    // Top genomic-position bar chart
    let gene_values: Vec<&DmrMetric> = dmr_metrics
        .iter()
        .filter(|r| r.start >= gene.start0 && r.start < gene.end0)
        .collect();

    if !gene_values.is_empty() {
        let bar_bottom = args.top_margin as f64 - 55.0;
        let bar_top = bar_bottom - args.top_bar_height;

        let min_value = gene_values
            .iter()
            .map(|r| r.value)
            .fold(f64::INFINITY, f64::min);

        let max_value = gene_values
            .iter()
            .map(|r| r.value)
            .fold(f64::NEG_INFINITY, f64::max);

        // Axis line
        writeln!(
            svg,
            r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#333" stroke-width="1"/>"##,
            x0,
            bar_bottom,
            x1,
            bar_bottom
        )?;

        // Label
        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="11" fill="#333">-log(p-value)</text>"##,
            x0 - 10.0,
            bar_bottom - args.top_bar_height / 2.0 + 4.0
        )?;

        // Min/max labels
        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="9" fill="#555">{:.3}</text>"##,
            x0 - 10.0,
            bar_bottom + 3.0,
            min_value
        )?;

        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="9" fill="#555">{:.3}</text>"##,
            x0 - 10.0,
            bar_top + 3.0,
            max_value
        )?;

        // Bars
        for r in gene_values {
            // Only draw if the coordinate maps into an exon block or compressed intron range.
            // For values defined at exonic positions, this will land exactly on the aligned exon.
            let px = map_genomic_boundary_to_x(&layout, r.start);

            let value_frac = if max_value > min_value {
                ((r.value - min_value) / (max_value - min_value)).clamp(0.0, 1.0)
            } else {
                0.5
            };

            let bar_h = value_frac * args.top_bar_height;
            let bar_y = bar_bottom - bar_h;
            let color = bar_color(r.value, min_value, max_value);

            let bar_w = 2.0;
            let bar_x = px - bar_w / 2.0;

            writeln!(
                svg,
                r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{}" stroke="none"/>"##,
                bar_x,
                bar_y,
                bar_w,
                bar_h.max(1.0),
                color
            )?;

            writeln!(
                svg,
                r##"<title>{}:{} value={:.4}</title>"##,
                escape_xml(&gene.chrom),
                r.start,
                r.value
            )?;
        }
    }

    let exon_h = (args.track_height as f64 * 0.55).max(8.0);
    // let heat_h = (args.track_height as f64 * 0.42).max(6.0);

    if let Some(gene_name) = gene.gene_name.as_ref() {
        writeln!(
            svg,
            r#"<text x="18" y="28" font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="bold">{}</text>"#,
            escape_xml(&format!(
                "{} ({gene_name}) ({}:{}-{}, strand {})",
                gene.gene_id,
                gene.chrom,
                gene.start0 + 1,
                gene.end0,
                gene.strand
            ))
        )?;
    } else {
        writeln!(
            svg,
            r#"<text x="18" y="28" font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="bold">{}</text>"#,
            escape_xml(&format!(
                "{}  ({}:{}-{}, strand {})",
                gene.gene_id,
                gene.chrom,
                gene.start0 + 1,
                gene.end0,
                gene.strand
            ))
        )?;
    }

    if has_score_thresh {
        writeln!(
            svg,
            r##"<text x="18" y="47" font-family="Arial, Helvetica, sans-serif" font-size="11" fill="#f96161">Only displaying sites with score >= {min_score:.4}</text>"##
        )?;
    }
    // Alternating background bands for merged exon blocks
    let plot_top = args.top_margin as f64 - 8.0;
    let plot_bottom = height as f64 - args.bottom_margin as f64 + 8.0;
    // Highlight specific genomic positions
    for pos0 in args.marker_positions.as_ref().unwrap_or(&Vec::new()) {
        if *pos0 < gene.start0 || *pos0 >= gene.end0 {
            continue;
        }

        let x = map_genomic_boundary_to_x(&layout, *pos0);

        writeln!(
            svg,
            r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#111" stroke-width="1.5" stroke-dasharray="4,3"/>"##,
            x,
            plot_top,
            x,
            plot_bottom
        )?;

        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#111">{}</text>"##,
            x - 3.0,
            plot_top + 4.0,
            format_u64_with_commas(*pos0 + 1)
        )?;
    }
    // for (i, block) in layout.blocks.iter().enumerate() {
    //     let bx1 = map_genomic_boundary_to_x(&layout, block.start0);
    //     let bx2 = map_genomic_boundary_to_x(&layout, block.end0);
    //     let fill = if i % 2 == 0 { "#fafafa" } else { "#f4f4f4" };
    //     writeln!(
    //         svg,
    //         r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{}" stroke="none"/>"##,
    //         bx1,
    //         plot_top,
    //         (bx2 - bx1).max(1.0),
    //         plot_bottom - plot_top,
    //         fill
    //     )?;
    // }

    // Vertical guide lines at merged exon boundaries
    // for block in &layout.blocks {
    //     let bx1 = map_genomic_boundary_to_x(&layout, block.start0);
    //     let bx2 = map_genomic_boundary_to_x(&layout, block.end0);
    //     writeln!(
    //         svg,
    //         r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#c9a5a5" stroke-width="1"/>"##,
    //         bx1, plot_top, bx1, plot_bottom
    //     )?;
    //     writeln!(
    //         svg,
    //         r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#c9a5a5" stroke-width="1"/>"##,
    //         bx2, plot_top, bx2, plot_bottom
    //     )?;
    // }
    // // Exon block genomic labels
    // if args.label_exons {
    //     for (i, block) in layout.blocks.iter().enumerate() {
    //         let bx1 = map_genomic_boundary_to_x(&layout, block.start0);
    //         let bx2 = map_genomic_boundary_to_x(&layout, block.end0);
    //         let bx_mid = (bx1 + bx2) / 2.0;

    //         let label = format_genomic_interval_label(&gene.chrom, block.start0, block.end0);

    //         // If the exon is very narrow, rotate the label to avoid overlap.
    //         let exon_px_width = bx2 - bx1;
    //         if exon_px_width >= 90.0 {
    //             writeln!(
    //                 svg,
    //                 r##"<text x="{:.2}" y="{:.2}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#555">exon block {}: {}</text>"##,
    //                 bx_mid,
    //                 plot_top - 8.0,
    //                 i + 1,
    //                 escape_xml(&label)
    //             )?;
    //         } else {
    //             writeln!(
    //                 svg,
    //                 r##"<text x="{:.2}" y="{:.2}" transform="rotate(-45 {:.2} {:.2})" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="9" fill="#555">E{} {}</text>"##,
    //                 bx_mid,
    //                 plot_top - 8.0,
    //                 bx_mid,
    //                 plot_top - 8.0,
    //                 i + 1,
    //                 escape_xml(&label)
    //             )?;
    //         }
    //     }
    // }

    // Top axis
    writeln!(
        svg,
        r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#333" stroke-width="1"/>"##,
        x0, axis_y, x1, axis_y
    )?;

    // Tick labels: gene start/end + merged exon block boundaries if not too many
    let mut tick_positions: Vec<u64> = vec![gene.start0, gene.end0];
    // if layout.blocks.len() <= 6 {
    //     for b in &layout.blocks {
    //         tick_positions.push(b.start0);
    //         tick_positions.push(b.end0);
    //     }
    // }
    tick_positions.sort_unstable();
    tick_positions.dedup();

    for pos0 in tick_positions {
        let px = map_genomic_boundary_to_x(&layout, pos0);
        writeln!(
            svg,
            r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#333" stroke-width="1"/>"##,
            px,
            axis_y,
            px,
            axis_y + 6.0
        )?;
        let label = if pos0 == gene.end0 {
            format_u64_with_commas(pos0)
        } else {
            format_u64_with_commas(pos0 + 1)
        };
        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="11" fill="#333">{}</text>"##,
            px,
            axis_y - 6.0,
            label
        )?;
    }

    // Legend
    // Methylation color legend
    let legend_w = 160.0;
    let legend_h = 10.0;
    let legend_x = x1 - legend_w;
    let legend_y = 16.0;
    let steps = 80usize;

    for i in 0..steps {
        let frac = i as f64 / (steps.saturating_sub(1)) as f64;
        let color = methylation_color(frac * 100.0);
        let lx = legend_x + frac * legend_w;

        writeln!(
            svg,
            r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{}" stroke="none"/>"##,
            lx,
            legend_y,
            legend_w / steps as f64 + 0.5,
            legend_h,
            color
        )?;
    }

    writeln!(
        svg,
        r##"<text x="{:.2}" y="{:.2}" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#333">percent modification ({})</text>"##,
        legend_x,
        legend_y - 3.0,
        single_code,
    )?;

    writeln!(
        svg,
        r##"<text x="{:.2}" y="{:.2}" text-anchor="start" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#333">0%</text>"##,
        legend_x,
        legend_y + legend_h + 12.0
    )?;

    writeln!(
        svg,
        r##"<text x="{:.2}" y="{:.2}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#333">50%</text>"##,
        legend_x + legend_w / 2.0,
        legend_y + legend_h + 12.0
    )?;

    writeln!(
        svg,
        r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#333">100%</text>"##,
        legend_x + legend_w,
        legend_y + legend_h + 12.0
    )?;

    for (row_idx, tx_id) in tx_ids.iter().enumerate() {
        let tx = transcripts.get(tx_id).unwrap();
        let y_top =
            args.top_margin as f64 + row_idx as f64 * (args.track_height + args.track_gap) as f64;
        let y_mid = y_top + args.track_height as f64 / 2.0;
        let exon_y = y_mid - exon_h / 2.0;
        // let heat_y = y_mid - heat_h / 2.0;

        // Transcript label
        writeln!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#111">{}({:.2}% rel)</text>"##,
            args.left_margin as f64 - 10.0,
            y_mid + 4.0,
            escape_xml(&tx_id.to_string()),
            rel_abundances.get(tx_id).copied().unwrap_or(f64::NAN) * 100f64
        )?;

        // Intron connectors only between consecutive exons
        for pair in tx.exons.windows(2) {
            let left_ex = &pair[0];
            let right_ex = &pair[1];

            let x_left = map_genomic_boundary_to_x(&layout, left_ex.end0);
            let x_right = map_genomic_boundary_to_x(&layout, right_ex.start0);

            writeln!(
                svg,
                r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="#666" stroke-width="1.2"/>"##,
                x_left, y_mid, x_right, y_mid
            )?;

            if x_right - x_left > 10.0 {
                draw_intron_break(&mut svg, (x_left + x_right) / 2.0, y_mid)?;
            }
        }

        // Exons
        for exon in &tx.exons {
            let ex1 = map_genomic_boundary_to_x(&layout, exon.start0);
            let ex2 = map_genomic_boundary_to_x(&layout, exon.end0);
            let ew = (ex2 - ex1).max(1.0);
            writeln!(
                svg,
                r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="#cfcfcf" stroke="#7a7a7a" stroke-width="0.6"/>"##,
                ex1, exon_y, ew, exon_h
            )?;
        }

        // Colored methylation ticks inside exon boxes
        if let Some(rows) = meth_by_tx.get(tx_id) {
            // debug!("got {} rows", rows.len());
            for r in rows {
                if r.stop() <= gene.start0 || r.start() >= gene.end0 {
                    // debug!("skipping {r:?}");
                    continue;
                }
                let metric_at_pos  = dmr_metrics.iter().find(|x| x.start == r.start());
                if metric_at_pos.map(|x| x.value < min_score).unwrap_or(true) {
                    continue;
                }
                // Use midpoint for longer intervals; for per-base BED rows this equals genomic_start.
                let pos0 = modification_tick_position(r);

                // Only draw ticks that fall inside this transcript's exons.
                let in_this_transcript = tx.exons.iter().any(|e| {
                    pos0 >= e.start0 && pos0 < e.end0
                });

                if !in_this_transcript {
                    continue;
                }

                let px = map_genomic_boundary_to_x(&layout, pos0);
                let color = methylation_color((r.frac_modified() * 100f32) as f64);

                // Make very narrow ticks still visible.
                let tick_w = 2.0_f64;
                let tick_x = px - tick_w / 2.0;

                writeln!(
                    svg,
                    r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{}" stroke="#222" stroke-width="0.25"/>"##,
                    tick_x,
                    exon_y,
                    tick_w,
                    exon_h,
                    color
                )?;

                // Tooltip in SVG viewers/browsers
                writeln!(
                    svg,
                    r##"<title>{}:{}-{} {}% methylation</title>"##,
                    escape_xml(&gene.chrom),
                    r.start()+ 1,
                    r.stop(),
                    r.frac_modified() * 100f32
                )?;
            }
        }
        // Strand hint
        if let (Some(first), Some(last)) = (tx.exons.first(), tx.exons.last()) {
            let arrow_size = 5.0;
            if tx.strand == '+' {
                let px = map_genomic_boundary_to_x(&layout, last.end0);
                writeln!(
                    svg,
                    r##"<path d="M {:.2} {:.2} l {:.2} {:.2} l {:.2} {:.2}" fill="none" stroke="#444" stroke-width="1"/>"##,
                    px - arrow_size,
                    y_mid - arrow_size / 2.0,
                    arrow_size,
                    arrow_size / 2.0,
                    -arrow_size,
                    arrow_size / 2.0
                )?;
            } else {
                let px = map_genomic_boundary_to_x(&layout, first.start0);
                writeln!(
                    svg,
                    r##"<path d="M {:.2} {:.2} l {:.2} {:.2} l {:.2} {:.2}" fill="none" stroke="#444" stroke-width="1"/>"##,
                    px + arrow_size,
                    y_mid - arrow_size / 2.0,
                    -arrow_size,
                    arrow_size / 2.0,
                    arrow_size,
                    arrow_size / 2.0
                )?;
            }
        }
    }

    writeln!(svg, "</svg>")?;
    Ok(svg)
}

pub(super) fn group_bedmethyl_records_by_tx(
    bedmethyl_records: Vec<BedMethylLine>,
) -> FxHashMap<GtfTranscript, Vec<BedMethylLine>> {
    bedmethyl_records.into_iter().fold(FxHashMap::default(), |mut agg, bml| {
        let txid = if let Some(tx_id) = bml.chrom.split("|").next() {
            &tx_id
        } else {
            bml.chrom.as_str()
        };
        if let Ok(gtf_tx) = GtfTranscript::from_str(&txid) {
            agg.entry(gtf_tx).or_insert_with(Vec::new).push(bml);
        }
        agg
    })
}

pub(super) struct GeneTxDmr {
    gene: GeneCommonCoord,
    records: Option<Vec<GeneIsoformDmrRecord<GeneDmrScore>>>,
}

impl GeneTxDmr {
    pub(super) fn empty() -> Self {
        Self { gene: GeneCommonCoord::empty(), records: None }
    }

    pub(super) fn clear(&mut self) {
        self.gene = GeneCommonCoord::empty();
        self.records = None;
    }

    pub(super) fn write(
        &self,
        writer: &mut Box<dyn IoWrite>,
        single_mod_code: Option<ModCodeRepr>,
        emit_full_results: bool,
    ) -> anyhow::Result<usize> {
        let mut records_written = 0usize;
        if let Some(records) = self.records.as_ref() {
            for record in records {
                let row = record.to_row(
                    &self.gene.chrom,
                    &self.gene.gene_id,
                    self.gene.gene_name.as_ref(),
                    single_mod_code,
                    emit_full_results,
                );
                writer.write(row.as_bytes())?;
                records_written = records_written.saturating_add(1);
            }
        }
        Ok(records_written)
    }

    pub(super) fn topk_records(
        &mut self,
        top_k: usize,
        min_effect_size: f64,
        labels: &HashSet<String>,
    ) -> Vec<PooledMethylationGenomePosition> {
        assert!(top_k > 0);
        if let Some(records) = self.records.as_mut() {
            let gene_string = self.gene.gene_id.to_string();
            let gene_name = self.gene.gene_name.clone();
            let label_point =
                gene_name.as_ref().is_some_and(|x| labels.contains(x));
            records.sort_by(|x, y| {
                let sx = x.score.p_value.ln().neg();
                let sy = y.score.p_value.ln().neg();
                sx.partial_cmp(&sy).unwrap_or(std::cmp::Ordering::Equal)
            });
            records
                .into_iter()
                .filter(|x| x.effect_size().abs() >= min_effect_size)
                .take(top_k)
                .map(|x| PooledMethylationGenomePosition {
                    genome_start: x.start,
                    gene: gene_string.clone(),
                    gene_name: gene_name.clone(),
                    log_fold_change: x.log_fc(),
                    neg_log_pvalue: x.score.p_value.ln().neg(),
                    effect_size: x.effect_size(),
                    label_point,
                })
                .collect()
        } else {
            Vec::new()
        }
    }
}

pub(super) struct GeneTxDmrWorker {
    pair_handler: PairTranscriptomeBedmethylHandler,
    transcript_models: Arc<FxHashMap<GtfTranscript, TranscriptModel>>,
    single_mod_code: Option<ModCodeRepr>,
}

impl GeneTxDmrWorker {
    pub(super) fn new<R: AsRef<Path>>(
        cond_a_path: R,
        cond_b_path: R,
        handler: &PairTranscriptomeBedmethylHandler,
        min_valid_coverage: u64,
        transcript_models: Arc<FxHashMap<GtfTranscript, TranscriptModel>>,
        single_mod_code: Option<ModCodeRepr>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let pair_handler = handler.get_copy(
            cond_a_path,
            cond_b_path,
            min_valid_coverage,
            multi_progress.clone(),
        )?;

        Ok(Self { pair_handler, transcript_models, single_mod_code })
    }
    pub(super) fn process_gene<const NMODS: usize>(
        &mut self,
        mut mem: GeneTxDmr,
        gene: GeneCommonCoord,
    ) -> MkResult<GeneTxDmr> {
        let records = run_gene_dmr::<NMODS>(
            &self.transcript_models,
            &gene,
            &mut self.pair_handler,
            self.single_mod_code,
        )?;

        mem.gene = gene;
        mem.records = records;
        Ok(mem)
    }
}

impl GeneIsoformDmrRecord<GeneDmrScore> {
    pub(super) fn to_row(
        &self,
        chrom: &str,
        gene_id: &GtfGene,
        gene_name: Option<&String>,
        single_mod_code: Option<ModCodeRepr>,
        emit_full_results: bool,
    ) -> String {
        let mut row = format!(
            "\
            {}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}{TAB}{}",
            chrom,
            self.start,
            self.end,
            self.score.mod_codes.iter().join(","),
            self.score.lrt,
            self.strand.to_char(),
            self.score.p_value,
            gene_id,
            gene_name.unwrap_or(&"-".to_string())
        );
        let log2_fc_field = if single_mod_code.is_some() {
            let log2_fc = self.log_fc();
            let effect_size = self.effect_size();
            format!("{TAB}{log2_fc}{TAB}{effect_size}")
        } else {
            "".to_string()
        };
        row.push_str(&log2_fc_field);

        if emit_full_results {
            assert_eq!(
                self.score.cond_a_proportions.len(),
                self.score.mod_codes.len() + 1
            );
            assert_eq!(
                self.score.cond_a_proportions.len(),
                self.score.cond_b_proportions.len()
            );
            assert_eq!(
                self.score.cond_a_proportions.len(),
                self.score.cond_a_counts.len()
            );
            assert_eq!(
                self.score.cond_a_proportions.len(),
                self.score.cond_b_counts.len()
            );

            let cond_a_unmodified = self.score.cond_a_proportions[0];
            let mut cond_a_props = hash_map! {
                "unmodified".to_string() => cond_a_unmodified
            };
            let cond_b_unmodified = self.score.cond_b_proportions[0];
            let mut cond_b_props = hash_map! {
                "unmodified".to_string() => cond_b_unmodified
            };

            for (i, mod_code) in self.score.mod_codes.iter().enumerate() {
                let idx = i + 1;
                let cond_a_prop = self.score.cond_a_proportions[idx];
                cond_a_props.insert(mod_code.to_string(), cond_a_prop);
                let cond_b_prop = self.score.cond_b_proportions[idx];
                cond_b_props.insert(mod_code.to_string(), cond_b_prop);
            }

            let cond_a_unmodified_counts = self.score.cond_a_counts[0];
            let mut cond_a_counts = hash_map! {
                "unmodified".to_string() => cond_a_unmodified_counts
            };
            let cond_b_unmodified_counts = self.score.cond_b_counts[0];
            let mut cond_b_counts = hash_map! {
                "unmodified".to_string() => cond_b_unmodified_counts
            };

            for (i, mod_code) in self.score.mod_codes.iter().enumerate() {
                let idx = i + 1;
                let cond_a_count = self.score.cond_a_counts[idx];
                cond_a_counts.insert(mod_code.to_string(), cond_a_count);
                let cond_b_count = self.score.cond_b_counts[idx];
                cond_b_counts.insert(mod_code.to_string(), cond_b_count);
            }
            let ser_a_props = serde_json::to_string(&cond_a_props).unwrap();
            let ser_b_props = serde_json::to_string(&cond_b_props).unwrap();
            let ser_a_counts = serde_json::to_string(&cond_a_counts).unwrap();
            let ser_b_counts = serde_json::to_string(&cond_b_counts).unwrap();
            let full = [ser_a_props, ser_b_props, ser_a_counts, ser_b_counts]
                .join("\t");
            row.push_str(&format!("{TAB}{full}"));
        }
        row.push('\n');
        row
    }
    pub(super) fn header(
        single_mod_code: Option<ModCodeRepr>,
        emit_full_results: bool,
    ) -> String {
        let mut h = vec![
            "#chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "p_value",
            "gene_id",
            "gene_name",
        ];
        if single_mod_code.is_some() {
            h.extend_from_slice(&["log2_fold_change", "effect_size"]);
        }
        if emit_full_results {
            h.extend_from_slice(&[
                "cond_a_proportions",
                "cond_b_proportions",
                "cond_a_counts",
                "cond_b_counts",
            ]);
        }
        let mut h = h.join("\t");
        h.push('\n');
        h
    }

    fn effect_size(&self) -> f64 {
        assert_eq!(self.score.cond_a_proportions.len(), 2);
        assert_eq!(self.score.cond_b_proportions.len(), 2);

        if self.score.cond_a_counts[1] == 0u32
            && self.score.cond_b_counts[1] == 0u32
        {
            return 0f64;
        }

        let m_a = self.score.cond_a_counts[1] as f64;
        let m_b = self.score.cond_b_counts[1] as f64;
        let u_a = self.score.cond_a_counts[0] as f64;
        let u_b = self.score.cond_b_counts[0] as f64;

        (m_a / (m_a + u_a)) - (m_b / (m_b + u_b))
    }

    fn log_fc(&self) -> f64 {
        assert_eq!(self.score.cond_a_proportions.len(), 2);
        assert_eq!(self.score.cond_b_proportions.len(), 2);

        if self.score.cond_a_counts[1] == 0u32
            && self.score.cond_b_counts[1] == 0u32
        {
            return 0f64;
        }

        let m_a = self.score.cond_a_counts[1] as f64;
        let m_b = self.score.cond_b_counts[1] as f64;
        let u_a = self.score.cond_a_counts[0] as f64;
        let u_b = self.score.cond_b_counts[0] as f64;

        if (m_a / (m_a + u_a)) - (m_b / (m_b + u_b)) == 0f64 {
            return 0f64;
        }

        let alpha = 0.5;
        let p_a = (m_a + alpha) / (m_a + u_a + 2.0 * alpha);
        let p_b = (m_b + alpha) / (m_b + u_b + 2.0 * alpha);
        (p_a / p_b).log2()
    }
}

#[derive(Debug)]
pub(super) struct GeneDmrScore {
    pub lrt: f64,
    pub p_value: f64,
    pub mod_codes: Vec<ModCodeRepr>,
    pub cond_a_proportions: Vec<f64>,
    pub cond_b_proportions: Vec<f64>,
    pub cond_a_counts: Vec<u32>,
    pub cond_b_counts: Vec<u32>,
}

impl GeneDmrScore {
    fn from_score_result(
        mut res: DirichletMultinomialTestResult,
        mod_codes: Vec<ModCodeRepr>,
    ) -> Self {
        assert_eq!(res.isoform_counts.len(), 2);
        let cond_a_counts = res.isoform_counts.pop().unwrap();
        let cond_b_counts = res.isoform_counts.pop().unwrap();
        let cond_a_proportions =
            scoring::proportions_from_counts(&cond_a_counts, 0f64);
        let cond_b_proportions =
            scoring::proportions_from_counts(&cond_b_counts, 0f64);
        Self {
            lrt: res.lrt,
            p_value: res.p_value,
            mod_codes,
            cond_a_proportions,
            cond_b_proportions,
            cond_a_counts,
            cond_b_counts,
        }
    }
}

pub(super) struct PooledMethylationGenomePosition {
    genome_start: u64,
    gene: String,
    gene_name: Option<String>,
    log_fold_change: f64,
    neg_log_pvalue: f64,
    effect_size: f64,
    label_point: bool,
}

#[rustfmt::skip]
pub(super) fn volcano_svg(points: &[PooledMethylationGenomePosition], title: Option<&String>) -> String {
    let width = 900.0;
    let height = 650.0;

    let margin_left = 80.0;
    let margin_right = 40.0;
    let margin_top = 40.0;
    let margin_bottom = 80.0;

    let plot_width = width - margin_left - margin_right;
    let plot_height = height - margin_top - margin_bottom;

    let valid_points: Vec<_> = points
        .iter()
        .filter(|p| {
            p.log_fold_change.is_finite() && p.neg_log_pvalue.is_finite()
        })
        .collect();

    if valid_points.is_empty() {
        return format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
  <text x="50%" y="50%" text-anchor="middle" font-family="sans-serif" font-size="20">
    No valid points to plot
  </text>
</svg>"#
        );
    }

    let max_abs_x = valid_points
        .iter()
        .map(|p| p.log_fold_change.abs())
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let min_x = -max_abs_x;
    let max_x = max_abs_x;

    let max_y = valid_points
        .iter()
        .map(|p| p.neg_log_pvalue)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let min_y = 0.0;

    let x_scale =
        |x: f64| margin_left + ((x - min_x) / (max_x - min_x)) * plot_width;

    let y_scale = |y: f64| {
        margin_top + plot_height - ((y - min_y) / (max_y - min_y)) * plot_height
    };

    let mut svg = String::new();

    if let Some(title) = title.as_ref() {
        svg.push_str(&format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
      <rect width="100%" height="100%" fill="white"/>

      <text x="{cx}" y="28" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="bold">
        {title}
      </text>
    "#,
            cx = width / 2.0,
            title = title,
        ));
    } else {
        svg.push_str(&format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
      <rect width="100%" height="100%" fill="white"/>
    "#,
        ));

    }

    // Axes
    // let x_axis_y = y_scale(0.0);
    let y_axis_x = x_scale(0.0);

    svg.push_str(&format!(
        r##"
  <line x1="{ml}" y1="{xb}" x2="{xr}" y2="{xb}" stroke="black" stroke-width="2"/>
  <line x1="{ml}" y1="{mt}" x2="{ml}" y2="{xb}" stroke="black" stroke-width="2"/>

  <line x1="{yx}" y1="{mt}" x2="{yx}" y2="{xb}" stroke="#999" stroke-width="1" stroke-dasharray="4 4"/>
"##,
        ml = margin_left,
        mt = margin_top,
        xr = margin_left + plot_width,
        xb = margin_top + plot_height,
        yx = y_axis_x
    ));

    // Axis labels
    svg.push_str(&format!(
        r#"
  <text x="{cx}" y="{y}" text-anchor="middle" font-family="sans-serif" font-size="16">
    log2 fold change
  </text>

  <text x="22" y="{cy}" text-anchor="middle" font-family="sans-serif" font-size="16"
        transform="rotate(-90 22 {cy})">
    -log p-value
  </text>
"#,
        cx = margin_left + plot_width / 2.0,
        y = height - 25.0,
        cy = margin_top + plot_height / 2.0
    ));

    // X ticks
    let tick_count = 6;
    for i in 0..=tick_count {
        let t = i as f64 / tick_count as f64;
        let value = min_x + t * (max_x - min_x);
        let x = x_scale(value);

        svg.push_str(&format!(
            r#"
  <line x1="{x:.2}" y1="{y1:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="black"/>
  <text x="{x:.2}" y="{ty:.2}" text-anchor="middle" font-family="sans-serif" font-size="12">
    {value:.2}
  </text>
"#,
            y1 = margin_top + plot_height,
            y2 = margin_top + plot_height + 6.0,
            ty = margin_top + plot_height + 22.0
        ));
    }

    // Y ticks
    for i in 0..=tick_count {
        let t = i as f64 / tick_count as f64;
        let value = min_y + t * (max_y - min_y);
        let y = y_scale(value);

        svg.push_str(&format!(
            r#"
  <line x1="{x1:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="black"/>
  <text x="{tx:.2}" y="{y_text:.2}" text-anchor="end" font-family="sans-serif" font-size="12">
    {value:.2}
  </text>
"#,
            x1 = margin_left - 6.0,
            x2 = margin_left,
            tx = margin_left - 10.0,
            y_text = y + 4.0
        ));
    }

    // Points and optional labels
    for p in valid_points {
        let x = x_scale(p.log_fold_change);
        let y = y_scale(p.neg_log_pvalue);

        let color = if p.log_fold_change.abs() >= 1.0 && p.neg_log_pvalue >= 2.0 {
            "#d62728"
        } else {
            "#4a90e2"
        };

        svg.push_str(&format!(
            r#"
      <circle cx="{x:.2}" cy="{y:.2}" r="4" fill="{color}" fill-opacity="0.75">
        <title>{gene} | {gene_name}genome_start={genome_start} | logFC={logfc:.4} | effect_size={effect_size} | -logP={neglogp:.4}</title>
      </circle>
    "#,
            gene = escape_xml(&p.gene),
            genome_start = p.genome_start,
            gene_name = p.gene_name.as_ref().map(|x| format!("{x} | ")).unwrap_or("".to_string()),
            logfc = p.log_fold_change,
            neglogp = p.neg_log_pvalue,
            effect_size = p.effect_size
        ));

        if p.label_point {
            svg.push_str(&format!(
                r##"
      <line x1="{x:.2}" y1="{y:.2}" x2="{lx:.2}" y2="{ly:.2}" stroke="#666" stroke-width="1"/>
      <text x="{tx:.2}" y="{ty:.2}" font-family="sans-serif" font-size="12" fill="black">
        {gene}
      </text>
    "##,
                lx = x + 8.0,
                ly = y - 8.0,
                tx = x + 10.0,
                ty = y - 10.0,
                gene = escape_xml(p.gene_name.as_ref().unwrap_or(&"-".to_string()))
            ));
        }
    }

    svg.push_str("</svg>\n");
    svg
}
