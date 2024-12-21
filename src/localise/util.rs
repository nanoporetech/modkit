use anyhow::{anyhow, bail};
use charming::component::{
    Axis, DataZoom, DataZoomType, Feature, Legend, Restore, SaveAsImage, Title,
    Toolbox, ToolboxDataZoom,
};
use charming::datatype::{CompositeValue, DataPoint, NumericValue};
use charming::element::{AxisType, LineStyle, Orient, Symbol};
use charming::series::Line;
use charming::{Chart, HtmlRenderer};
use clap::ValueEnum;
use indicatif::{MultiProgress, ProgressIterator};
use itertools::{Itertools, MinMaxResult};
use prettytable::row;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;
use crate::tabix::HtsTabixHandler;
use crate::util::{
    get_subroutine_progress_bar, GenomeRegion, ModPositionInfo, StrandRule,
};

#[derive(Default)]
pub(super) struct LocalizedModCounts {
    offsets: FxHashMap<ModCodeRepr, FxHashMap<i64, ModPositionInfo<u64>>>,
}

impl LocalizedModCounts {
    fn add_bedmethyl_record(
        &mut self,
        bed_methyl_line: &BedMethylLine,
        anchor_point: u64,
    ) {
        let pos = bed_methyl_line.start() as i64;
        let offset = (anchor_point as i64).saturating_sub(pos);
        let mod_pos_info = self
            .offsets
            .entry(bed_methyl_line.raw_mod_code)
            .or_insert(FxHashMap::default())
            .entry(offset)
            .or_insert(ModPositionInfo::new(0, 0));
        mod_pos_info.n_mod += bed_methyl_line.count_methylated;
        mod_pos_info.n_valid += bed_methyl_line.valid_coverage;
    }

    pub(super) fn get_table(&self, mpb: &MultiProgress) -> prettytable::Table {
        let pb = mpb.add(get_subroutine_progress_bar(self.offsets.len()));
        pb.set_message("compiling table");
        let mut tab = prettytable::Table::new();
        tab.set_titles(row![
            "mod_code",
            "offset",
            "n_valid",
            "n_mod",
            "percent_modified"
        ]);

        self.offsets
            .iter()
            .progress_with(pb)
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
            .fold(tab, |mut tab, (code, counts)| {
                counts
                    .iter()
                    .sorted_by(|(a, _), (b, _)| a.cmp(b))
                    .map(|(offset, info)| {
                        row![
                            *code,
                            *offset,
                            info.n_valid,
                            info.n_mod,
                            info.percent_modified()
                        ]
                    })
                    .for_each(|row| {
                        tab.add_row(row);
                    });
                tab
            })
    }

    pub(super) fn get_plot(
        &self,
        chart_name: Option<&String>,
    ) -> anyhow::Result<String> {
        let default_name = "modification_patterns".to_string();
        let xs = self
            .offsets
            .values()
            .flat_map(|counts| counts.keys().copied())
            .sorted()
            .collect::<Vec<i64>>();
        let (left, right) = match xs.iter().minmax() {
            MinMaxResult::MinMax(x, y) => (*x, *y),
            _ => bail!("should be at least one offset"),
        };
        let mut chart = Chart::new()
            .data_zoom(
                DataZoom::new()
                    .type_(DataZoomType::Slider)
                    .orient(Orient::Horizontal),
            )
            .data_zoom(
                DataZoom::new()
                    .type_(DataZoomType::Slider)
                    .orient(Orient::Vertical),
            )
            .legend(Legend::new())
            .title(Title::new().text(chart_name.unwrap_or(&default_name)))
            .toolbox(
                Toolbox::new().feature(
                    Feature::new()
                        .data_zoom(ToolboxDataZoom::new().y_axis_index("none"))
                        .restore(Restore::new())
                        .save_as_image(SaveAsImage::new()),
                ),
            )
            .x_axis(
                Axis::new()
                    .type_(AxisType::Value)
                    .min(left)
                    .max(right)
                    .name("offset"),
            )
            .y_axis(
                Axis::new().type_(AxisType::Value).name("percent modified"),
            );
        for (mod_code, counts) in self.offsets.iter() {
            let dat = counts
                .iter()
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
                .map(|(offset, info)| {
                    DataPoint::Value(CompositeValue::Array(vec![
                        CompositeValue::Number(NumericValue::Integer(*offset)),
                        CompositeValue::Number(NumericValue::Float(
                            info.frac_modified() as f64,
                        )),
                    ]))
                })
                .collect::<Vec<DataPoint>>();
            chart = chart.series(
                Line::new()
                    .name(format!("{mod_code}"))
                    .data(dat)
                    .symbol(Symbol::None)
                    .line_style(LineStyle::new().width(1.5)),
            );
        }

        HtmlRenderer::new(chart_name.unwrap_or(&default_name), 800, 800)
            .render(&chart)
            .map_err(|e| anyhow!("failed to render, {e:?}"))
    }
}

impl Moniod for LocalizedModCounts {
    fn zero() -> Self {
        Self::default()
    }

    fn op(self, other: Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        for (mod_code, counts) in other.offsets {
            let agg =
                self.offsets.entry(mod_code).or_insert(FxHashMap::default());
            for (offset, info) in counts {
                if let Some(x) = agg.get_mut(&offset) {
                    x.n_valid += info.n_valid;
                    x.n_mod += info.n_mod;
                } else {
                    agg.insert(offset, info);
                }
            }
        }
    }

    fn len(&self) -> usize {
        self.offsets.values().map(|x| x.keys().len()).sum::<usize>()
    }
}

impl GenomeRegion {
    pub(super) fn into_localized_mod_counts(
        self,
        index: &HtsTabixHandler<BedMethylLine>,
        strand_rule: Option<StrandRule>,
        stranded_features: Option<StrandedFeatures>,
        io_threads: usize,
    ) -> anyhow::Result<LocalizedModCounts> {
        let bedmethyl_records = index.fetch_region(
            &self.chrom,
            &(self.start..self.end),
            strand_rule.unwrap_or(self.strand),
            io_threads,
        )?;
        let anchor_point = self.midpoint();
        let loc_counts = bedmethyl_records
            .into_par_iter()
            .filter(|bm| {
                stranded_features
                    .map(|f| {
                        let overlaps = self.strand.overlaps(&bm.strand);
                        match f {
                            StrandedFeatures::Same => overlaps,
                            StrandedFeatures::Opposite => !overlaps,
                        }
                    })
                    .unwrap_or(true)
            })
            .fold(
                || LocalizedModCounts::zero(),
                |mut acc, next| {
                    acc.add_bedmethyl_record(&next, anchor_point);
                    acc
                },
            )
            .reduce(|| LocalizedModCounts::zero(), |a, b| a.op(b));
        Ok(loc_counts)
    }
}

#[derive(Debug, Copy, Clone, ValueEnum)]
pub(super) enum StrandedFeatures {
    #[clap(name = "same")]
    Same,
    #[clap(name = "opposite")]
    Opposite,
}
