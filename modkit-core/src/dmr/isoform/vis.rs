use std::collections::HashMap;

use derive_new::new;
use std::fmt::Write;

use crate::dmr::{
    bedmethyl::BedMethylLine,
    isoform::{Exon, GeneModel, GtfTranscript},
};

#[derive(new)]
pub(super) struct GenomeSortedTranscriptModel {
    pub strand: char,     // '+' or '-'
    pub exons: Vec<Exon>, // ordered in transcript 5'->3'
}

pub(super) fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

pub(super) fn format_u64_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    let len = s.len();
    for (i, ch) in s.chars().enumerate() {
        out.push(ch);
        let rem = len - i - 1;
        if rem > 0 && rem % 3 == 0 {
            out.push(',');
        }
    }
    out
}
#[derive(Debug, Clone)]
pub(super) struct CompressedLayout {
    pub blocks: Vec<Exon>, // merged union-of-exons blocks, genomic order
    pub block_x_starts: Vec<f64>, // x start for each block
    pub exon_scale: f64,   // px per genomic bp inside exons
    pub intron_px: f64,    // fixed px width for each intron gap
    pub left: f64,
    pub inner_width: f64,
}

pub(super) fn merged_gene_exon_blocks(
    gene: &GeneModel,
    transcripts: &HashMap<GtfTranscript, GenomeSortedTranscriptModel>,
) -> Vec<Exon> {
    let mut all: Vec<Exon> = Vec::new();

    for tx_id in &gene.transcript_ids {
        if let Some(tx) = transcripts.get(tx_id) {
            all.extend(tx.exons.iter().cloned());
        }
    }

    merge_half_open_exons(&all)
}

pub(super) fn merge_half_open_exons(exons: &[Exon]) -> Vec<Exon> {
    if exons.is_empty() {
        return Vec::new();
    }

    let mut v = exons.to_vec();
    v.sort_by_key(|e| e.start0);

    let mut merged = Vec::new();
    let mut cur = v[0].clone();

    for e in v.into_iter().skip(1) {
        // half-open intervals: merge overlaps and abutting segments
        if e.start0 <= cur.end0 {
            cur.end0 = cur.end0.max(e.end0);
        } else {
            merged.push(cur);
            cur = e;
        }
    }
    merged.push(cur);
    merged
}

pub(super) fn build_compressed_layout(
    blocks: &[Exon],
    left: f64,
    inner_width: f64,
    intron_px_pref: f64,
) -> CompressedLayout {
    let n_introns = blocks.len().saturating_sub(1);
    let total_exonic_bp: u64 = blocks.iter().map(|b| b.end0 - b.start0).sum();

    // Prevent introns from consuming too much of the figure.
    let intron_px = if n_introns == 0 {
        0.0
    } else {
        intron_px_pref.min(inner_width * 0.40 / n_introns as f64)
    };

    let exon_px_total = inner_width - intron_px * n_introns as f64;
    let exon_scale = if total_exonic_bp == 0 {
        1.0
    } else {
        exon_px_total / total_exonic_bp as f64
    };

    let mut block_x_starts = Vec::with_capacity(blocks.len());
    let mut x = left;
    for (i, b) in blocks.iter().enumerate() {
        block_x_starts.push(x);
        x += (b.end0 - b.start0) as f64 * exon_scale;
        if i + 1 < blocks.len() {
            x += intron_px;
        }
    }

    CompressedLayout {
        blocks: blocks.to_vec(),
        block_x_starts,
        exon_scale,
        intron_px,
        left,
        inner_width,
    }
}

/// Maps a 0-based genomic boundary position into SVG x.
/// This is piecewise:
/// - linear within exon blocks
/// - linearly interpolated across compressed introns
pub(super) fn map_genomic_boundary_to_x(
    layout: &CompressedLayout,
    pos0: u64,
) -> f64 {
    if layout.blocks.is_empty() {
        return layout.left;
    }

    let first = &layout.blocks[0];
    let last = &layout.blocks[layout.blocks.len() - 1];

    if pos0 <= first.start0 {
        return layout.left;
    }
    if pos0 >= last.end0 {
        return layout.left + layout.inner_width;
    }

    for i in 0..layout.blocks.len() {
        let b = &layout.blocks[i];
        let x0 = layout.block_x_starts[i];
        let x1 = x0 + (b.end0 - b.start0) as f64 * layout.exon_scale;

        // inside block, including right boundary
        if pos0 >= b.start0 && pos0 <= b.end0 {
            return x0
                + (pos0.saturating_sub(b.start0)) as f64 * layout.exon_scale;
        }

        // inside intron between this block and the next
        if i + 1 < layout.blocks.len() {
            let next = &layout.blocks[i + 1];
            if pos0 > b.end0 && pos0 < next.start0 {
                let intron_len = next.start0 - b.end0;
                if intron_len == 0 {
                    return x1;
                }
                let frac = (pos0 - b.end0) as f64 / intron_len as f64;
                return x1 + frac * layout.intron_px;
            }
        }
    }

    layout.left + layout.inner_width
}

pub(super) fn draw_intron_break(
    svg: &mut String,
    x_mid: f64,
    y: f64,
) -> anyhow::Result<()> {
    let dx = 4.0;
    let dy = 5.0;
    writeln!(
        svg,
        r##"<path d="M {:.2} {:.2} l {:.2} {:.2} M {:.2} {:.2} l {:.2} {:.2}" fill="\#555" stroke="#555" stroke-width="5"/>"##,
        x_mid - dx,
        y - dy / 2.0,
        dx,
        dy,
        x_mid + 1.0,
        y - dy / 2.0,
        dx,
        dy
    )?;
    Ok(())
}

// pub(super) fn midpoint_u64(a: u64, b: u64) -> u64 {
//     a + (b - a) / 2
// }

// pub(super) fn format_genomic_interval_label(
//     chrom: &str,
//     start0: u64,
//     end0: u64,
// ) -> String {
//     // Display as 1-based inclusive genomic interval
//     format!(
//         "{}:{}-{}",
//         chrom,
//         format_u64_with_commas(start0 + 1),
//         format_u64_with_commas(end0)
//     )
// }

pub(super) fn modification_tick_position(row: &BedMethylLine) -> u64 {
    // For a BED-style interval, use the midpoint.
    // For per-base rows, genomic_end = genomic_start + 1, so this returns
    // genomic_start. midpoint_u64(row.start(), row.stop())
    row.start()
}

pub(super) fn methylation_color(percent: f64) -> String {
    // Clamp to 0..100
    let p = percent.clamp(0.0, 100.0) / 100.0;

    // Simple blue -> white -> red gradient
    // 0%   = blue
    // 50%  = white
    // 100% = red
    let (r, g, b) = if p <= 0.5 {
        let t = p / 0.5;
        let r = (255.0 * t).round() as u8;
        let g = (255.0 * t).round() as u8;
        let b = 255_u8;
        (r, g, b)
    } else {
        let t = (p - 0.5) / 0.5;
        let r = 255_u8;
        let g = (255.0 * (1.0 - t)).round() as u8;
        let b = (255.0 * (1.0 - t)).round() as u8;
        (r, g, b)
    };

    format!("#{:02x}{:02x}{:02x}", r, g, b)
}

pub(super) fn bar_color(value: f64, min_value: f64, max_value: f64) -> String {
    let t = if max_value > min_value {
        ((value - min_value) / (max_value - min_value)).clamp(0.0, 1.0)
    } else {
        0.5
    };

    // light gray -> dark purple
    let r = (230.0 + t * (84.0 - 230.0)).round() as u8;
    let g = (230.0 + t * (39.0 - 230.0)).round() as u8;
    let b = (230.0 + t * (143.0 - 230.0)).round() as u8;

    format!("#{:02x}{:02x}{:02x}", r, g, b)
}
