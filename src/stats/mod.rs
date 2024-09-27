use derive_new::new;
use prettytable::{cell, row};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::bedmethyl::BedMethylLine;
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;
use crate::tabix::HtsTabixHandler;
use crate::util::{GenomeRegion, ModPositionInfo, StrandRule};

pub mod subcommand;

#[derive(new)]
struct MethylationStats {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub strand: StrandRule,
    pub per_mod_methylation: FxHashMap<ModCodeRepr, ModPositionInfo<u64>>,
}

impl MethylationStats {
    fn into_row(self, mod_codes: &[ModCodeRepr]) -> prettytable::Row {
        let name = self.name.unwrap_or_else(|| String::from("."));
        let mut row = row![self.chrom, self.start, self.end, name, self.strand];
        mod_codes.iter().for_each(|code| {
            let (count, valid, pct) =
                if let Some(count) = self.per_mod_methylation.get(code) {
                    (count.n_mod, count.n_valid, count.percent_modified())
                } else {
                    (0, 0, 0f32)
                };
            row.add_cell(cell!(count));
            row.add_cell(cell!(valid));
            row.add_cell(cell!(pct));
        });
        row
    }

    fn header(mod_codes: &[ModCodeRepr]) -> prettytable::Row {
        let mut row = row!["chrom", "start", "end", "name", "strand"];
        for code in mod_codes {
            row.add_cell(cell!(format!("count_{code}")));
            row.add_cell(cell!(format!("count_valid_{code}")));
            row.add_cell(cell!(format!("percent_{code}")));
        }
        row
    }
}

impl GenomeRegion {
    fn into_stats(
        self,
        index: &HtsTabixHandler<BedMethylLine>,
        min_coverage: u64,
        mod_codes: Option<&FxHashSet<ModCodeRepr>>,
        io_threads: usize,
        // todo allow combine
    ) -> anyhow::Result<MethylationStats> {
        let bm_lines = index.fetch_region(
            &self.chrom,
            &(self.start..self.end),
            self.strand,
            io_threads,
        )?;

        let mod_counts = bm_lines
            .par_iter()
            .filter(|rec| rec.valid_coverage >= min_coverage)
            .filter(|rec| rec.strand.overlaps(&self.strand))
            .filter(|rec| {
                mod_codes
                    .map(|codes| codes.contains(&rec.raw_mod_code))
                    .unwrap_or(true)
            })
            .fold(
                || FxHashMap::default(),
                |mut agg, next| {
                    let count = agg
                        .entry(next.raw_mod_code)
                        .or_insert(ModPositionInfo::zero());
                    count.n_mod += next.count_methylated;
                    count.n_valid += next.valid_coverage;

                    agg
                },
            )
            .reduce(|| FxHashMap::default(), |a, b| a.op(b));

        Ok(MethylationStats::new(
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.strand,
            mod_counts,
        ))
    }
}
