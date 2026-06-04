use std::{collections::HashSet, path::Path};

use anyhow::bail;
use bitvec::bitvec;
use indexmap::IndexMap;
use indicatif::MultiProgress;
use itertools::Itertools;
use log::debug;
use rustc_hash::FxHashMap;

use crate::{
    dmr::{
        bedmethyl::BedMethylLine,
        isoform::{
            dirichlet_multinomial_lrt, genomic0_to_gene_common0,
            group_bedmethyl_records_by_tx, transcript_pos0_to_genomic0,
            GeneCommonCoord, GeneDmrScore, GeneIsoformDmrRecord, GtfGene,
            GtfTranscript, ModCodeCount, TranscriptBedmethylHandler,
            TranscriptModel,
        },
    },
    errs::MkResult,
    mod_base_code::{ModCodeRepr, MOD_CODE_TO_DNA_BASE},
    util::Strand,
};

type PairBedmethylRecords = (Vec<BedMethylLine>, Vec<BedMethylLine>);

pub(crate) struct PairTranscriptomeBedmethylHandler {
    cond_a: TranscriptBedmethylHandler,
    cond_b: TranscriptBedmethylHandler,
}

impl PairTranscriptomeBedmethylHandler {
    pub(crate) fn from_handlers(
        mut cond_a: TranscriptBedmethylHandler,
        mut cond_b: TranscriptBedmethylHandler,
        genes_from_gtf: &mut IndexMap<GtfGene, GeneCommonCoord>,
    ) -> anyhow::Result<Self> {
        let common_genes = cond_a
            .gene_id_to_transcript_ids
            .iter()
            .filter(|(x, _)| cond_b.gene_id_to_transcript_ids.contains_key(*x))
            .map(|(x, _)| x)
            .cloned()
            .collect::<HashSet<GtfGene>>();
        if common_genes.is_empty() {
            bail!("no genes in common in two conditions")
        }
        let common_genes = common_genes
            .into_iter()
            .filter(|x| genes_from_gtf.contains_key(x))
            .collect::<HashSet<GtfGene>>();
        if common_genes.is_empty() {
            bail!("no genes in common in two conditions and GTF")
        }

        cond_a
            .gene_id_to_transcript_ids
            .retain(|x, _| common_genes.contains(x));
        cond_b
            .gene_id_to_transcript_ids
            .retain(|x, _| common_genes.contains(x));
        genes_from_gtf.retain(|x, _| common_genes.contains(x));

        Ok(Self { cond_a, cond_b })
    }

    fn get_bedmethyl_records_for_gene(
        &mut self,
        gene_id: &GtfGene,
        single_mod_code: Option<ModCodeRepr>,
    ) -> MkResult<Option<PairBedmethylRecords>> {
        let records_a = self
            .cond_a
            .get_read_bedmethyl_for_gene(gene_id, single_mod_code)?;
        let records_b = self
            .cond_b
            .get_read_bedmethyl_for_gene(gene_id, single_mod_code)?;
        Ok(records_a.and_then(|x| records_b.map(|y| (x, y))))
    }

    pub(super) fn get_copy<R: AsRef<Path>>(
        &self,
        cond_a_path: R,
        cond_b_path: R,
        min_valid_coverage: u64,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let cond_a = TranscriptBedmethylHandler::from_path_with_lookup_table(
            cond_a_path,
            self.cond_a.gene_id_to_transcript_ids.clone(),
            min_valid_coverage,
            multi_progress.clone(),
        )?;
        let cond_b = TranscriptBedmethylHandler::from_path_with_lookup_table(
            cond_b_path,
            self.cond_b.gene_id_to_transcript_ids.clone(),
            min_valid_coverage,
            multi_progress,
        )?;
        Ok(Self { cond_a, cond_b })
    }
}

fn add_to_modcount(mod_counts: &mut [ModCodeCount], bml: &BedMethylLine) {
    'counts: for mod_count in mod_counts {
        match mod_count {
            ModCodeCount::Count(mod_code, count)
                if *mod_code == bml.raw_mod_code =>
            {
                *count = count.saturating_add(bml.count_methylated as u32);
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

fn add_to_modcount_single_code(
    mod_counts: &mut [ModCodeCount],
    unmodified_counts: &mut [u32],
    common_start: usize,
    bml: &BedMethylLine,
) {
    match &mut mod_counts[0] {
        ModCodeCount::Empty => {
            let x = ModCodeCount::Count(
                bml.raw_mod_code,
                bml.count_methylated as u32,
            );
            mod_counts[0] = x;
            unmodified_counts[common_start] = unmodified_counts[common_start]
                .saturating_add(bml.count_other as u32);
        }
        ModCodeCount::Count(mod_code, count) => {
            assert_eq!(*mod_code, bml.raw_mod_code);
            *count = count.saturating_add(bml.count_methylated as u32);
            unmodified_counts[common_start] = unmodified_counts[common_start]
                .saturating_add(bml.count_other as u32);
        }
    }
}

pub(super) fn run_gene_dmr<const NMODS: usize>(
    transcripts: &FxHashMap<GtfTranscript, TranscriptModel>,
    gene: &GeneCommonCoord,
    handler: &mut PairTranscriptomeBedmethylHandler,
    single_mod_code: Option<ModCodeRepr>,
) -> MkResult<Option<Vec<GeneIsoformDmrRecord<GeneDmrScore>>>> {
    let Some((cond_a_records, cond_b_records)) = handler
        .get_bedmethyl_records_for_gene(&gene.gene_id, single_mod_code)?
    else {
        debug!(
            "None bedmethyl records for {}",
            gene.gene_name.as_ref().unwrap()
        );
        return Ok(None);
    };
    let cond_a_records = group_bedmethyl_records_by_tx(cond_a_records);
    let cond_b_records = group_bedmethyl_records_by_tx(cond_b_records);

    let tx_iter = cond_a_records.iter().filter_map(|(tx, records)| {
        match (transcripts.get(tx), cond_b_records.get(tx)) {
            (Some(tm), Some(b_recs)) => Some((tm, records, b_recs)),
            _ => None,
        }
    });

    let mut valid_pos = bitvec![0; gene.total_len as usize];
    let mut genomic_positions = vec![0u64; gene.total_len as usize];
    let mut unmodified_a = vec![0u32; gene.total_len as usize];
    let mut unmodified_b = vec![0u32; gene.total_len as usize];
    let mut modified_counts_a =
        vec![[ModCodeCount::Empty; NMODS]; gene.total_len as usize];
    let mut modified_counts_b =
        vec![[ModCodeCount::Empty; NMODS]; gene.total_len as usize];

    'transcripts: for (tm, a_records, b_records) in tx_iter {
        if a_records.is_empty() || b_records.is_empty() {
            continue;
        }

        let mut joined = JoinByPos::new(a_records, b_records).peekable();
        let Some((first_a, first_b)) = joined.peek() else {
            continue 'transcripts;
        };
        assert_eq!(first_a.start(), first_b.start());
        let mut tx_start = first_a.start();
        let mut genomic_start = transcript_pos0_to_genomic0(tm, tx_start)?;
        let mut common_start =
            genomic0_to_gene_common0(gene, genomic_start)? as usize;
        if valid_pos[common_start] {
            assert_eq!(genomic_positions[common_start], genomic_start)
        } else {
            valid_pos.set(common_start, true);
            genomic_positions[common_start] = genomic_start
        }
        unmodified_a[common_start] = unmodified_a[common_start]
            .saturating_add(first_a.count_canonical as u32);
        unmodified_b[common_start] = unmodified_b[common_start]
            .saturating_add(first_b.count_canonical as u32);

        for (a, b) in joined {
            assert_eq!(a.start(), b.start());
            if a.start() != tx_start {
                tx_start = a.start();
                genomic_start = transcript_pos0_to_genomic0(tm, tx_start)?;
                common_start =
                    genomic0_to_gene_common0(gene, genomic_start)? as usize;
                valid_pos.set(common_start, true);
                genomic_positions[common_start] = genomic_start;
                unmodified_a[common_start] = unmodified_a[common_start]
                    .saturating_add(a.count_canonical as u32);
                unmodified_b[common_start] = unmodified_b[common_start]
                    .saturating_add(b.count_canonical as u32);
            }

            let mod_counts_a = &mut modified_counts_a[common_start];
            let mod_counts_b = &mut modified_counts_b[common_start];

            if single_mod_code.is_some() {
                add_to_modcount_single_code(
                    mod_counts_a,
                    &mut unmodified_a,
                    common_start,
                    a,
                );
                add_to_modcount_single_code(
                    mod_counts_b,
                    &mut unmodified_b,
                    common_start,
                    b,
                );
            } else {
                add_to_modcount(mod_counts_a, a);
                add_to_modcount(mod_counts_b, b);
            }
        }
    }
    let gpos_iter =
        valid_pos.iter().zip(genomic_positions.iter()).enumerate().filter_map(
            |(i, (b, gpos))| if *b { Some((i, *gpos)) } else { None },
        );
    let mut out_records = Vec::new();
    'gpos: for (i, gpos) in gpos_iter {
        let mut pos_counts_a = vec![unmodified_a[i]];
        let mut pos_counts_b = vec![unmodified_b[i]];
        let mut mod_codes_a = Vec::new();
        let mut mod_codes_b = Vec::new();

        let pos_mod_counts_a = modified_counts_a[i];
        'a_mods: for x in pos_mod_counts_a {
            match x {
                ModCodeCount::Empty => break 'a_mods,
                ModCodeCount::Count(mod_code_repr, count) => {
                    mod_codes_a.push(mod_code_repr);
                    pos_counts_a.push(count);
                }
            }
        }

        let pos_mod_counts_b = modified_counts_b[i];
        'b_mods: for x in pos_mod_counts_b {
            match x {
                ModCodeCount::Empty => break 'b_mods,
                ModCodeCount::Count(mod_code_repr, count) => {
                    mod_codes_b.push(mod_code_repr);
                    pos_counts_b.push(count);
                }
            }
        }

        if mod_codes_a == mod_codes_b
            && mod_codes_a
                .iter()
                .map(|x| MOD_CODE_TO_DNA_BASE.get(x))
                .unique()
                .count()
                == 1
        {
            let Some(test_result) = dirichlet_multinomial_lrt(
                &[pos_counts_a, pos_counts_b],
                0.05,
                0.5,
            ) else {
                continue 'gpos;
            };
            let score =
                GeneDmrScore::from_score_result(test_result, mod_codes_a);
            out_records.push(GeneIsoformDmrRecord {
                start: gpos,
                end: gpos.saturating_add(1),
                strand: Strand::parse_char(gene.strand).unwrap(),
                n_transcripts: 2,
                score: score,
            });
        } else {
            continue 'gpos;
        }
    }

    if gene.strand == '-' {
        Ok(Some(out_records.into_iter().rev().collect()))
    } else {
        Ok(Some(out_records))
    }
}

struct JoinByPos<'a> {
    left: &'a [BedMethylLine],
    right: &'a [BedMethylLine],
    i: usize,
    j: usize,
}

impl<'a> JoinByPos<'a> {
    fn new(left: &'a [BedMethylLine], right: &'a [BedMethylLine]) -> Self {
        Self { left, right, i: 0, j: 0 }
    }
}

impl<'a> Iterator for JoinByPos<'a> {
    type Item = (&'a BedMethylLine, &'a BedMethylLine);

    fn next(&mut self) -> Option<Self::Item> {
        while self.i < self.left.len() && self.j < self.right.len() {
            let l = &self.left[self.i];
            let r = &self.right[self.j];

            match l.start().cmp(&r.start()) {
                std::cmp::Ordering::Less => {
                    self.i += 1;
                }
                std::cmp::Ordering::Greater => {
                    self.j += 1;
                }
                std::cmp::Ordering::Equal => {
                    self.i += 1;
                    self.j += 1;

                    return Some((l, r));
                }
            }
        }

        None
    }
}
