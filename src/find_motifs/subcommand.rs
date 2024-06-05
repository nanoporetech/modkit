use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use anyhow::bail;
use bio::io::fasta::Reader as FastaReader;
use clap::Args;
use indicatif::{MultiProgress, ParallelProgressIterator};
use itertools::Itertools;
use log::info;
use prettytable::{row, Table};
use rayon::prelude::*;

use crate::find_motifs::{
    find_motifs_for_mod, load_bedmethyl, EnrichedMotif, EnrichedMotifData,
    KmerModificationDb, MotifRelationship,
};
use crate::logging::init_logging;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::util::{get_subroutine_progress_bar, get_ticker};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryFindMotifs {
    /// Input bedmethyl table, can be used directly from modkit pileup.
    #[arg(short = 'i', long)]
    in_bedmethyl: PathBuf,
    /// Optionally output a machine-parsable TSV (human-readable table will
    /// always be output to the log).
    #[arg(short = 'o', long)]
    out_table: Option<PathBuf>,
    /// Optionally output machine parsable table with known motif
    /// modification frequencies that were not found during search.
    #[arg(long = "known-motifs-table", requires = "known_motifs")]
    out_known_table: Option<PathBuf>,
    /// Number of threads to use.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Reference sequence in FASTA format used for the pileup.
    /// Reference sequence in FASTA format used for the pileup.
    #[arg(short = 'r', long = "ref")]
    reference_fasta: PathBuf,
    /// Fraction modified threshold below which consider a genome location to
    /// be "low modification".
    #[arg(long = "low-thresh", default_value_t = 0.2)]
    low_threshold: f32,
    /// Fraction modified threshold above which consider a genome location to
    /// be "high modification" or enriched for modification.
    #[arg(long = "high-thresh", default_value_t = 0.6)]
    high_threshold: f32,
    /// Minimum log-odds to consider a motif sequence to be enriched.
    #[arg(long, default_value_t = 1.5)]
    min_log_odds: f32,
    /// Minimum log-odds to consider a motif sequence to be enriched when
    /// performing exhaustive search.
    #[arg(long, conflicts_with = "skip_search", default_value_t = 2.5)]
    exhaustive_seed_min_log_odds: f32,
    /// Exhaustive search seed length, increasing this value increases
    /// computational time.
    #[arg(long, conflicts_with = "skip_search", default_value_t = 3)]
    exhaustive_seed_len: usize,
    /// Skip the exhaustive search phase, saves time but the results may be
    /// less sensitive
    #[arg(long, default_value_t = false)]
    skip_search: bool,
    /// Minimum coverage in the bedMethyl to consider a record valid.
    #[arg(long, default_value_t = 5)]
    min_coverage: u64,
    /// Upstream and downstream number of bases to search for a motif sequence
    /// around a modified base. Example: --context-size 12 12.
    #[arg(long, num_args=2, default_values_t=vec![12, 12])]
    context_size: Vec<u64>,
    /// Initial "fixed" seed window size in base pairs around the modified
    /// base. Example: --init-context-size 2 2
    #[arg(long, num_args=2, default_values_t=vec![2, 2])]
    init_context_size: Vec<usize>,
    /// Minimum number of total sites in the genome required for a motif to be
    /// considered.
    #[arg(long, default_value_t = 300)]
    min_sites: u64,
    /// Minimum fraction of sites in the genome to be "high-modification" for a
    /// motif to be considered.
    #[arg(long = "min-frac-mod", default_value_t = 0.85)]
    frac_sites_thresh: f32,
    /// Gather enrichment information for a known motif as well as compare to
    /// discovered motifs. Format should be <sequence> <offset> <mod_code>.
    #[arg(long="known-motif", num_args = 3, action = clap::ArgAction::Append)]
    known_motifs: Option<Vec<String>>,
    /// Specify which modification codes to process, default will process all
    /// modification codes found in the input bedMethyl file
    #[arg(long = "mod-code")]
    mod_codes: Option<Vec<String>>,
    /// Force override SAM specification of association of modification codes
    /// to primary sequence bases.
    #[arg(long = "force-override-spec", default_value_t = false)]
    override_spec: bool,
    /// Output log to this file.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Disable the progress bars.
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
}

impl EntryFindMotifs {
    fn load_references(
        &self,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<HashMap<String, Vec<u8>>> {
        info!("loading references from {:?}", self.reference_fasta);
        let pb = multi_progress.add(get_ticker());
        pb.set_message("sequences read");
        let reader = FastaReader::from_file(&self.reference_fasta)?;

        let (contigs, n_fails) = reader.records().fold(
            (HashMap::new(), 0usize),
            |(mut agg, fails), record| match record {
                Ok(r) => {
                    let record_name = r.id().to_string();
                    let seq = r
                        .seq()
                        .iter()
                        .map(|&nt| nt.to_ascii_uppercase())
                        .collect::<Vec<u8>>();
                    agg.insert(record_name, seq);
                    pb.inc(1);
                    (agg, fails)
                }
                Err(_) => (agg, fails + 1),
            },
        );

        if n_fails > 0 {
            info!("failed to load {n_fails} record(s)");
        }

        if contigs.is_empty() {
            bail!("failed to read any reference sequences");
        } else {
            pb.finish_and_clear();
            info!("loaded {} sequence(s)", contigs.len());
            Ok(contigs)
        }
    }

    fn load_mod_db(
        &self,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<KmerModificationDb> {
        let reference_sequences = self.load_references(multi_progress)?;
        let context = [self.context_size[0], self.context_size[1]];
        for x in context {
            if x > 127u64 {
                bail!("context cannot be larger than 127x2 (255) bases")
            }
        }

        // as well
        load_bedmethyl(
            &self.in_bedmethyl,
            self.min_coverage,
            context,
            self.low_threshold,
            self.high_threshold,
            &reference_sequences,
            multi_progress,
        )
    }

    fn parse_known_motifs(
        &self,
        mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
    ) -> anyhow::Result<Option<Vec<EnrichedMotif>>> {
        let context_size =
            [self.context_size[0] as usize, self.context_size[1] as usize];
        self.known_motifs
            .as_ref()
            .map(|raw_motifs| {
                let all_motifs = raw_motifs
                    .chunks(3)
                    .map(|parts| {
                        EnrichedMotif::new_from_parts(
                            &parts[0],
                            &parts[2],
                            &parts[1],
                            context_size,
                            mod_code_lookup,
                        )
                    })
                    .collect::<anyhow::Result<Vec<EnrichedMotif>>>();

                all_motifs.map(|xs| {
                    xs.into_iter()
                        .collect::<HashSet<_>>()
                        .into_iter()
                        .collect::<Vec<_>>()
                })
            })
            .transpose()
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
        let _ = init_logging(self.log_filepath.as_ref());
        if self.context_size.len() != 2 {
            bail!("context-size must be 2 elements")
        }
        if self.init_context_size.len() != 2 {
            bail!("init-context-size must be 2 elements")
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let mod_db = self.load_mod_db(&mpb)?;
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
                            self.min_log_odds,
                            self.min_sites,
                            self.frac_sites_thresh,
                            self.skip_search,
                            self.exhaustive_seed_len,
                            self.exhaustive_seed_min_log_odds,
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
            let context_size =
                [self.context_size[0] as usize, self.context_size[1] as usize];
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
                    .map(|motif| {
                        let (total_high_count, total_low_count) =
                            mod_db.get_total_mod_counts(&motif);
                        let total_mid_count = mod_db.get_mid_counts(&motif);
                        EnrichedMotifData::new(
                            motif.clone(),
                            total_high_count,
                            total_low_count,
                            total_mid_count,
                        )
                    })
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
        let context_size =
            [self.context_size[0] as usize, self.context_size[1] as usize];
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
