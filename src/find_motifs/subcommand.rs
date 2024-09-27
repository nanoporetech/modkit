use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::RwLock;

use anyhow::{anyhow, bail, Context};
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ParallelProgressIterator};
use itertools::Itertools;
use log::{debug, info};
use prettytable::{row, Table};
use rayon::{prelude::*, ThreadPool};
use rustc_hash::FxHashMap;

use crate::find_motifs::motif_bed::motif_bed;
use crate::find_motifs::{
    find_motifs_for_mod, load_bedmethyl_and_references, make_tables,
    merge_motifs, parse_known_motifs, parse_motifs_from_table,
    parse_raw_known_motifs, EnrichedMotif, EnrichedMotifData, KmerMask,
    KmerModificationDb, MotifRelationship,
};
use crate::logging::init_logging;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::util::get_subroutine_progress_bar;

use super::args::{InputArgs, KnownMotifsArgs, RefineArgs};

#[derive(Subcommand)]
pub enum EntryMotifs {
    /// Search for modification-enriched subsequences in a reference genome
    Search(EntryFindMotifs),
    /// Use a previously defined list of motif sequences and further refine
    /// them with a bedMethyl table.
    Refine(EntryRefineMotifs),
    /// Calculate enrichment statistics on a set of motifs from a bedMethyl
    /// table.
    Evaluate(EntrySearchMotifs), // todo rename
    /// Create BED file with all locations of a sequence motif.
    /// Example: modkit motif bed CG 0
    Bed(EntryMotifBed),
}

impl EntryMotifs {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryMotifs::Search(x) => x.run(),
            EntryMotifs::Evaluate(x) => x.run(),
            EntryMotifs::Refine(x) => x.run(),
            EntryMotifs::Bed(x) => x.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryFindMotifs {
    #[clap(flatten)]
    input_args: InputArgs,
    /// Optionally output a machine-parsable TSV (human-readable table will
    /// always be output to the log).
    #[arg(short = 'o', long)]
    out_table: Option<PathBuf>,
    /// Optionally output machine parsable table with known motif
    /// modification frequencies that were not found during search.
    #[arg(long = "eval-motifs-table")]
    out_known_table: Option<PathBuf>,
    #[clap(flatten)]
    refine_args: RefineArgs,
    /// Initial "fixed" seed window size in base pairs around the modified
    /// base. Example: --init-context-size 2 2
    #[arg(long, num_args=2, default_values_t=vec![2, 2])]
    init_context_size: Vec<usize>,
    /// Format should be <sequence> <offset> <mod_code>.
    #[arg(long="known-motif", num_args = 3, action = clap::ArgAction::Append)]
    pub known_motifs: Option<Vec<String>>,
    /// Path to known motifs in tabular format. Tab-separated values:
    /// <mod_code>\t<motif_seq>\t<offset>. May have the same header as the
    /// output table from this command.
    #[arg(long = "known-motifs-table")]
    pub known_motifs_table: Option<PathBuf>,
    /// Specify which modification codes to process, default will process all
    /// modification codes found in the input bedMethyl file
    #[arg(long = "mod-code")]
    mod_codes: Option<Vec<String>>,
    /// Force override SAM specification of association of modification codes
    /// to primary sequence bases.
    #[arg(long = "force-override-spec", default_value_t = false)]
    override_spec: bool,
}

impl EntryFindMotifs {
    fn get_context(&self) -> [u64; 2] {
        [self.refine_args.context_size[0], self.refine_args.context_size[1]]
    }

    fn load_mod_db(
        &self,
        multi_progress: &MultiProgress,
        pool: &ThreadPool,
    ) -> anyhow::Result<KmerModificationDb> {
        load_bedmethyl_and_references(
            &self.input_args.reference_fasta,
            &self.input_args.in_bedmethyl,
            self.input_args.contig.clone(),
            self.refine_args.min_coverage,
            self.get_context(),
            self.refine_args.low_threshold,
            self.refine_args.high_threshold,
            multi_progress,
            self.input_args.io_threads,
            pool,
        )
    }

    fn parse_known_motifs(
        &self,
        mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
    ) -> anyhow::Result<Option<Vec<EnrichedMotif>>> {
        let context_size = [
            self.refine_args.context_size[0] as usize,
            self.refine_args.context_size[1] as usize,
        ];
        let mut known_motifs = Vec::new();
        let cli_motifs = self
            .known_motifs
            .as_ref()
            .map(|raw_motifs| {
                parse_raw_known_motifs(
                    raw_motifs,
                    context_size,
                    mod_code_lookup,
                )
            })
            .transpose()?;
        let table_motifs = self
            .known_motifs_table
            .as_ref()
            .map(|tab_fp| {
                parse_motifs_from_table(tab_fp, context_size, mod_code_lookup)
            })
            .transpose()?;
        if let Some(cli_motifs) = cli_motifs {
            known_motifs.extend(cli_motifs);
        }

        if let Some(table_motifs) = table_motifs {
            known_motifs.extend(table_motifs);
        }

        if known_motifs.is_empty() {
            Ok(None)
        } else {
            Ok(Some(known_motifs))
        }
    }

    fn parse_input_mod_codes(
        &self,
    ) -> anyhow::Result<Option<Vec<ModCodeRepr>>> {
        self.mod_codes
            .as_ref()
            .map(|raw_codes| {
                raw_codes
                    .iter()
                    .map(|raw_code| ModCodeRepr::parse(raw_code))
                    .collect::<anyhow::Result<Vec<ModCodeRepr>>>()
                    .map(|parsed_codes| {
                        parsed_codes.into_iter().sorted().collect()
                    })
            })
            .transpose()
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.input_args.log_filepath.as_ref());
        if self.refine_args.context_size.len() != 2 {
            bail!("context-size must be 2 elements")
        }
        if self.init_context_size.len() != 2 {
            bail!("init-context-size must be 2 elements")
        }
        if self.out_known_table.is_some() {
            if self.known_motifs.is_none() && self.known_motifs_table.is_none()
            {
                bail!(
                    "--eval-motifs-table requires input known motifs with \
                     --known-motif and/or --known-motifs-table"
                )
            }
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.input_args.threads)
            .build()?;

        let mpb = MultiProgress::new();
        if self.input_args.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let mod_db = self.load_mod_db(&mpb, &pool)?;
        let input_mod_codes = self.parse_input_mod_codes()?;
        let inferred_mod_codes =
            mod_db.get_inferred_mod_code_associations(!self.override_spec)?;
        let known_motifs = self.parse_known_motifs(&inferred_mod_codes)?;
        if let Some(km) = known_motifs.as_ref() {
            info!(
                "parsed {} known motifs {}",
                km.len(),
                km.iter().map(|x| x.to_string()).join(",")
            );
        }

        let mod_codes = if let Some(codes) = input_mod_codes {
            let codes = codes
                .into_iter()
                .filter_map(|c| inferred_mod_codes.get_key_value(&c))
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
                .map(|(c, b)| (*c, *b))
                .collect::<Vec<(ModCodeRepr, DnaBase)>>();
            if codes.is_empty() {
                bail!(
                    "zero modification codes in common with requested and \
                     bedMethyl"
                )
            } else {
                codes
            }
        } else {
            inferred_mod_codes
                .iter()
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
                .map(|(c, b)| (*c, *b))
                .collect::<Vec<(ModCodeRepr, DnaBase)>>()
        };

        let results = {
            let mut results = pool.install(|| {
                mod_codes
                    .par_iter()
                    .flat_map(|(mod_code, canonical_base)| {
                        find_motifs_for_mod(
                            *canonical_base,
                            *mod_code,
                            &mod_db,
                            &self.init_context_size,
                            self.refine_args.min_log_odds,
                            self.refine_args.min_sites,
                            self.refine_args.frac_sites_thresh,
                            self.refine_args.skip_search,
                            self.refine_args.exhaustive_seed_len,
                            self.refine_args.exhaustive_seed_min_log_odds,
                            &mpb,
                        )
                    })
                    .collect::<Vec<EnrichedMotifData>>()
            });
            results.par_sort_by(|a, b| {
                b.frac_modified()
                    .partial_cmp(&a.frac_modified())
                    .expect("should compare")
            });
            results
        };

        let motifs_to_score = if let Some(known_motifs) = known_motifs.as_ref()
        {
            let mut motifs_to_score = Vec::new();
            let context_size = [
                self.refine_args.context_size[0] as usize,
                self.refine_args.context_size[1] as usize,
            ];
            let grouped_by_base =
                results.iter().fold(HashMap::new(), |mut agg, next_motif| {
                    let canonical_base = next_motif.motif.canonical_base;
                    agg.entry(canonical_base)
                        .or_insert(Vec::new())
                        .push(next_motif);
                    agg
                });
            let mut found = 0usize;
            for known_motif in known_motifs.iter() {
                let found_it = grouped_by_base
                    .get(&known_motif.canonical_base)
                    .and_then(|discovered| {
                        discovered.iter().find(|mot| {
                            match mot.motif.compare(known_motif, context_size) {
                                MotifRelationship::Equal => true,
                                _ => false,
                            }
                        })
                    })
                    .is_some();
                if found_it {
                    found += 1;
                } else {
                    motifs_to_score.push(known_motif)
                }
            }
            info!(
                "found {found} of {} known motifs, {} were not found and will \
                 be scored",
                known_motifs.len(),
                motifs_to_score.len()
            );
            let pb =
                mpb.add(get_subroutine_progress_bar(motifs_to_score.len()));
            pb.set_message("scoring known motifs");
            Some(
                motifs_to_score
                    .into_par_iter()
                    .progress_with(pb)
                    .map(|motif| mod_db.get_enriched_motif_data(motif))
                    .collect::<Vec<EnrichedMotifData>>(),
            )
        } else {
            None
        };

        let known_motifs_lookup = known_motifs.as_ref().map(|motifs| {
            motifs.iter().fold(HashMap::new(), |mut agg, next| {
                let primary_base = next.canonical_base;
                agg.entry(primary_base).or_insert(Vec::new()).push(next);
                agg
            })
        });
        let n_motifs = results.len();

        let results_table = self.format_human_readable_table(
            &results,
            known_motifs_lookup.as_ref(),
        );
        info!("Found {n_motifs} motifs:\n{results_table}");

        if let Some(out_fp) = self.out_table.as_ref() {
            let mach_table = self.format_machine_readable_table(
                &results,
                known_motifs_lookup.as_ref(),
            );
            let writer = csv::WriterBuilder::new()
                .has_headers(true)
                .delimiter('\t' as u8)
                .from_path(out_fp)?;
            mach_table.to_csv_writer(writer)?;
        }

        match motifs_to_score {
            Some(unfound_motifs) if !unfound_motifs.is_empty() => {
                let human_table = self
                    .format_unfound_motifs_human_readable_table(
                        &unfound_motifs,
                        &results,
                    );
                info!("Known motifs that were not found:\n{human_table}");
                if let Some(fp) = self.out_known_table.as_ref() {
                    let mach_table = self
                        .format_unfound_motifs_machine_readable_table(
                            &unfound_motifs,
                            &results,
                        );
                    let writer = csv::WriterBuilder::new()
                        .has_headers(true)
                        .delimiter('\t' as u8)
                        .from_path(fp)?;
                    mach_table.to_csv_writer(writer)?;
                }
            }
            Some(_) => {
                info!("All known motifs found.")
            }
            _ => {}
        }

        info!("finished");
        Ok(())
    }

    #[inline]
    fn get_closest_motif(
        &self,
        motif: &EnrichedMotif,
        others: &HashMap<DnaBase, Vec<&EnrichedMotif>>,
    ) -> (String, String) {
        let context_size = [
            self.refine_args.context_size[0] as usize,
            self.refine_args.context_size[1] as usize,
        ];
        if let Some(motifs_for_base) = others.get(&motif.canonical_base) {
            motifs_for_base
                .iter()
                .map(|m| motif.compare(m, context_size))
                .enumerate()
                .min_by(|(_, a), (_, b)| a.cmp(b))
                .map(|(idx, rel)| {
                    (motifs_for_base[idx].to_string(), rel.to_string())
                })
                .unwrap_or(("-".to_string(), "-".to_string()))
        } else {
            ("-".to_string(), "-".to_string())
        }
    }

    fn format_unfound_motifs_machine_readable_table(
        &self,
        unfound_motifs: &[EnrichedMotifData],
        discovered_motifs: &[EnrichedMotifData],
    ) -> Table {
        let organized_by_can_base =
            discovered_motifs.iter().fold(HashMap::new(), |mut agg, next| {
                agg.entry(next.motif.canonical_base)
                    .or_insert(Vec::new())
                    .push(&next.motif);
                agg
            });

        let mut tab = Table::new();
        tab.set_titles(row![
            "mod_code",
            "motif",
            "offset",
            "frac_mod",
            "high_count",
            "low_count",
            "mid_count",
            "status",
            "closest_found_motif",
        ]);
        for result in unfound_motifs {
            let mod_code = result.motif.multi_sequence.mod_code;
            let motif = result.motif.format_seq();
            let offset = result.motif.multi_sequence.get_offset();
            let frac_mod = result.frac_modified();
            let high_count = result.total_high_count;
            let low_count = result.total_low_count;
            let mid_count = result.total_mid_count;
            let (closest, relationship) =
                self.get_closest_motif(&result.motif, &organized_by_can_base);
            let row = row![
                mod_code,
                motif,
                offset,
                frac_mod,
                high_count,
                low_count,
                mid_count,
                relationship,
                closest
            ];
            tab.add_row(row);
        }
        tab
    }

    fn format_unfound_motifs_human_readable_table(
        &self,
        unfound_motifs: &[EnrichedMotifData],
        discovered_motifs: &[EnrichedMotifData],
    ) -> Table {
        let organized_by_can_base =
            discovered_motifs.iter().fold(HashMap::new(), |mut agg, next| {
                agg.entry(next.motif.canonical_base)
                    .or_insert(Vec::new())
                    .push(&next.motif);
                agg
            });

        let mut tab = Table::new();
        tab.set_format(
            *prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE,
        );
        tab.set_titles(row![
            "motif",
            "frac_mod",
            "high_count",
            "low_count",
            "mid_count",
            "status",
            "closest_found_motif",
        ]);
        for result in unfound_motifs {
            let motif_repr = result.motif.to_string();
            let frac_mod = result.frac_modified();
            let high_count = result.total_high_count;
            let low_count = result.total_low_count;
            let mid_count = result.total_mid_count;
            let (closest, relationship) =
                self.get_closest_motif(&result.motif, &organized_by_can_base);
            let row = row![
                motif_repr,
                frac_mod,
                high_count,
                low_count,
                mid_count,
                relationship,
                closest
            ];
            tab.add_row(row);
        }
        tab
    }

    fn format_machine_readable_table(
        &self,
        results: &[EnrichedMotifData],
        known_motifs: Option<&HashMap<DnaBase, Vec<&EnrichedMotif>>>,
    ) -> Table {
        let mut tab = Table::new();
        if known_motifs.is_some() {
            tab.set_titles(row![
                "mod_code",
                "motif",
                "offset",
                "frac_mod",
                "high_count",
                "low_count",
                "mid_count",
                "status",
                "closest_known_motif",
            ]);
        } else {
            tab.set_titles(row![
                "mod_code",
                "motif",
                "offset",
                "frac_mod",
                "high_count",
                "low_count",
                "mid_count",
            ]);
        }

        for result in results.iter() {
            let mod_code = result.motif.multi_sequence.mod_code;
            let motif = result.motif.format_seq();
            let offset = result.motif.multi_sequence.get_offset();
            let frac_mod = result.frac_modified();
            let high_count = result.total_high_count;
            let low_count = result.total_low_count;
            let mid_count = result.total_mid_count;
            let row = if let Some(km) = known_motifs.as_ref() {
                let (closest, relationship) =
                    self.get_closest_motif(&result.motif, km);
                row![
                    mod_code,
                    motif,
                    offset,
                    frac_mod,
                    high_count,
                    low_count,
                    mid_count,
                    relationship,
                    closest
                ]
            } else {
                row![
                    mod_code, motif, offset, frac_mod, high_count, low_count,
                    mid_count,
                ]
            };

            tab.add_row(row);
        }

        tab
    }

    fn format_human_readable_table(
        &self,
        results: &[EnrichedMotifData],
        known_motifs: Option<&HashMap<DnaBase, Vec<&EnrichedMotif>>>,
    ) -> Table {
        let mut tab = Table::new();
        tab.set_format(
            *prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE,
        );
        if known_motifs.is_some() {
            tab.set_titles(row![
                "motif",
                "frac_mod",
                "high_count",
                "low_count",
                "mid_count",
                "status",
                "closest_known_motif",
            ]);
        } else {
            tab.set_titles(row![
                "motif",
                "frac_mod",
                "high_count",
                "low_count",
                "mid_count",
            ]);
        }

        for result in results.iter() {
            let motif_repr = result.motif.to_string();
            let frac_mod = result.frac_modified();
            let high_count = result.total_high_count;
            let low_count = result.total_low_count;
            let mid_count = result.total_mid_count;
            let row = if let Some(km) = known_motifs.as_ref() {
                let (closest, relationship) =
                    self.get_closest_motif(&result.motif, km);
                row![
                    motif_repr,
                    frac_mod,
                    high_count,
                    low_count,
                    mid_count,
                    relationship,
                    closest
                ]
            } else {
                row![motif_repr, frac_mod, high_count, low_count, mid_count]
            };

            tab.add_row(row);
        }

        tab
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryRefineMotifs {
    #[clap(flatten)]
    input_args: InputArgs,
    #[clap(flatten)]
    known_motifs_args: KnownMotifsArgs,
    /// Machine-parsable table of refined motifs. Human-readable table always
    /// printed to stderr and log.
    #[arg(long = "out")]
    out_table: Option<PathBuf>,
    /// Minimum fraction of sites in the genome to be "high-modification"
    /// for a motif to be further refined, otherwise it will be discarded.
    #[arg(long = "min_refine_frac_mod", default_value_t = 0.6)]
    min_refine_frac_modified: f32,
    /// Minimum number of total sites in the genome required for a motif to be
    /// further refined, otherwise it will be discarded.
    #[arg(long, default_value_t = 300)]
    #[arg(long)]
    pub min_refine_sites: u64,
    #[clap(flatten)]
    refine_args: RefineArgs,
    /// Force override SAM specification of association of modification codes
    /// to primary sequence bases.
    #[arg(long = "force-override-spec", default_value_t = false)]
    override_spec: bool,
}

impl EntryRefineMotifs {
    fn refine_scored_motifs(
        &self,
        mod_db: &KmerModificationDb,
        scored_motifs: &[EnrichedMotifData],
        mpb: &MultiProgress,
        pool: ThreadPool,
    ) -> Vec<EnrichedMotifData> {
        let (motifs_to_refine, num_below_frac, num_below_min_sites, both) =
            scored_motifs.iter().fold(
                (Vec::new(), 0usize, 0usize, 0usize),
                |(mut agg, f, c, b), next| {
                    let enough_sites = next.enough_sites(self.min_refine_sites);
                    let enough_modified =
                        next.frac_modified() >= self.min_refine_frac_modified;
                    match (enough_sites, enough_modified) {
                        (true, true) => {
                            agg.push(next.motif.clone());
                            (agg, f, c, b)
                        }
                        (true, false) => (agg, f + 1, c, b),
                        (false, true) => (agg, f, c + 1, b),
                        (false, false) => (agg, f, c, b + 1),
                    }
                },
            );
        info!(
            "have {} motifs to refine, {} discarded",
            motifs_to_refine.len(),
            scored_motifs.len().saturating_sub(motifs_to_refine.len())
        );
        info!(
            "discard reasons:\n\tBelow fraction modified: \
             {num_below_frac}\n\tBelow min sites: \
             {num_below_min_sites}\n\tBelow both: {both}"
        );
        let refine_pb =
            mpb.add(get_subroutine_progress_bar(motifs_to_refine.len()));
        refine_pb.set_message("refining motifs");
        let cache = RwLock::new(FxHashMap::default());

        let refined_motifs = pool.install(|| {
            motifs_to_refine
                .par_iter()
                .progress_with(refine_pb)
                .map(|motif| {
                    let kmer_subset = mod_db.get_kmer_subset(
                        motif.canonical_base,
                        &KmerMask::default(),
                        motif.multi_sequence.mod_code,
                    );
                    motif.clone().refine(
                        &mod_db,
                        &cache,
                        &kmer_subset,
                        self.refine_args.min_sites,
                        self.refine_args.frac_sites_thresh,
                        self.refine_args.min_log_odds,
                    )
                })
                .collect::<Vec<EnrichedMotif>>()
        });

        info!("merging motifs");
        let merged_refined_motifs = merge_motifs(refined_motifs);
        info!(
            "have {} merged, refined motifs to score",
            merged_refined_motifs.len()
        );
        let re_score_pb =
            mpb.add(get_subroutine_progress_bar(merged_refined_motifs.len()));
        re_score_pb.set_message("scoring merged, refined, motifs");
        let refined_motifs = pool.install(|| {
            merged_refined_motifs
                .into_par_iter()
                .progress_with(re_score_pb)
                .map(|m| mod_db.get_enriched_motif_data(&m))
                .collect::<Vec<EnrichedMotifData>>()
        });
        refined_motifs
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.input_args.log_filepath.as_ref());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.input_args.threads)
            .build()
            .context("failed to make threadpool")?;
        let mpb = MultiProgress::new();
        if self.input_args.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        let context_bases = [
            self.refine_args.context_size[0],
            self.refine_args.context_size[1],
        ];
        let mod_db = load_bedmethyl_and_references(
            &self.input_args.reference_fasta,
            &self.input_args.in_bedmethyl,
            self.input_args.contig.clone(),
            self.refine_args.min_coverage,
            context_bases,
            self.refine_args.low_threshold,
            self.refine_args.high_threshold,
            &mpb,
            self.input_args.io_threads,
            &pool,
        )?;
        let inferred_mod_codes =
            mod_db.get_inferred_mod_code_associations(!self.override_spec)?;
        let motifs_to_evaluate = parse_known_motifs(
            &self.known_motifs_args,
            [context_bases[0] as usize, context_bases[1] as usize],
            &inferred_mod_codes,
        )?;
        if motifs_to_evaluate.is_empty() {
            bail!("failed to parse any motifs to evaluate")
        }

        info!("have {} motifs to evaluate", motifs_to_evaluate.len());

        let score_pb =
            mpb.add(get_subroutine_progress_bar(motifs_to_evaluate.len()));
        score_pb.set_message("scoring motifs");
        let scored_motifs = pool.install(|| {
            motifs_to_evaluate
                .par_iter()
                .progress_with(score_pb)
                .map(|motif| mod_db.get_enriched_motif_data(&motif))
                .collect::<Vec<EnrichedMotifData>>()
                .into_iter()
                .collect::<Vec<EnrichedMotifData>>()
        });

        let refined_motifs =
            self.refine_scored_motifs(&mod_db, &scored_motifs, &mpb, pool);

        let (refined_table, refined_mch_table) = make_tables(&refined_motifs);
        if let Some(p) = self.out_table.as_ref() {
            let writer = csv::WriterBuilder::new()
                .delimiter('\t' as u8)
                .from_path(p)
                .map_err(|e| {
                    anyhow!("failed to make TSV output writer at {p:?}, {e}")
                })?;
            refined_mch_table.to_csv_writer(writer)?;
        }
        info!("refined motifs:\n{refined_table}");
        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntrySearchMotifs {
    #[clap(flatten)]
    input_args: InputArgs,
    #[clap(flatten)]
    known_motifs_args: KnownMotifsArgs,
    /// Machine-parsable table of refined motifs. Human-readable table always
    /// printed to stderr and log.
    #[arg(long = "out")]
    out_table: Option<PathBuf>,
    /// Force override SAM specification of association of modification codes
    /// to primary sequence bases.
    #[arg(long = "force-override-spec", default_value_t = false)]
    override_spec: bool,
    /// Minimum coverage in the bedMethyl to consider a record valid.
    #[arg(long, default_value_t = 5)]
    min_coverage: u64,
    /// Upstream and downstream number of bases to search for a motif sequence
    /// around a modified base. Example: --context-size 12 12.
    #[arg(long, num_args=2, default_values_t=vec![12, 12])]
    context_size: Vec<u64>,
    /// Fraction modified threshold below which consider a genome location to
    /// be "low modification".
    #[arg(long = "low-thresh", default_value_t = 0.2)]
    low_threshold: f32,
    /// Fraction modified threshold above which consider a genome location to
    /// be "high modification" or enriched for modification.
    #[arg(long = "high-thresh", default_value_t = 0.6)]
    high_threshold: f32,
    /// Don't print final table to stderr (will still go to log file).
    #[arg(long, default_value_t = false)]
    suppress_table: bool,
}

impl EntrySearchMotifs {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.input_args.log_filepath.as_ref());
        if self.suppress_table
            && (self.out_table.is_none()
                && self.input_args.log_filepath.is_none())
        {
            bail!(
                "must provide an file to output table or a log file if \
                 suppressing human-readable table"
            )
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.input_args.threads)
            .build()
            .context("failed to make threadpool")?;
        let mpb = MultiProgress::new();
        if self.input_args.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        let context_bases = [self.context_size[0], self.context_size[1]];
        let mod_db = load_bedmethyl_and_references(
            &self.input_args.reference_fasta,
            &self.input_args.in_bedmethyl,
            self.input_args.contig.clone(),
            self.min_coverage,
            context_bases,
            self.low_threshold,
            self.high_threshold,
            &mpb,
            self.input_args.io_threads,
            &pool,
        )?;
        let inferred_mod_codes =
            mod_db.get_inferred_mod_code_associations(!self.override_spec)?;
        let motifs_to_evaluate = parse_known_motifs(
            &self.known_motifs_args,
            [context_bases[0] as usize, context_bases[1] as usize],
            &inferred_mod_codes,
        )?;
        if motifs_to_evaluate.is_empty() {
            bail!("failed to parse any motifs to evaluate")
        }

        info!("have {} motifs to evaluate", motifs_to_evaluate.len());

        let score_pb =
            mpb.add(get_subroutine_progress_bar(motifs_to_evaluate.len()));
        score_pb.set_message("evaluating motifs");
        let scored_motifs = pool.install(|| {
            motifs_to_evaluate
                .par_iter()
                .progress_with(score_pb)
                .map(|motif| mod_db.get_enriched_motif_data(&motif))
                .collect::<Vec<EnrichedMotifData>>()
                .into_iter()
                .collect::<Vec<EnrichedMotifData>>()
        });

        let (scored_table, scored_mch_table) = make_tables(&scored_motifs);

        if let Some(p) = self.out_table.as_ref() {
            let writer = csv::WriterBuilder::new()
                .delimiter('\t' as u8)
                .from_path(p)
                .map_err(|e| {
                    anyhow!("failed to make TSV output writer at {p:?}, {e}")
                })?;
            scored_mch_table.to_csv_writer(writer)?;
        }
        if self.suppress_table {
            debug!("evaluated motifs:\n{scored_table}");
        } else {
            info!("evaluated motifs:\n{scored_table}");
        }

        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryMotifBed {
    /// Input FASTA file
    fasta: PathBuf,
    /// Motif to search for within FASTA, e.g. CG
    motif: String,
    /// Offset within motif, e.g. 0
    offset: usize,
    /// Respect soft masking in the reference FASTA.
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
}

impl EntryMotifBed {
    fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(None);
        motif_bed(&self.fasta, &self.motif, self.offset, self.mask)
    }
}
