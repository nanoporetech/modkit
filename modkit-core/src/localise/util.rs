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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FocusRegion {
    region: GenomeRegion,
    anchor_point: u64,
}

impl FocusRegion {
    pub(super) fn from_feature(
        mut region: GenomeRegion,
        window: u64,
        chrom_length: u64,
    ) -> Self {
        let anchor_point = region.midpoint();
        region.start = anchor_point.saturating_sub(window).min(chrom_length);
        region.end = anchor_point
            .saturating_add(window)
            .saturating_add(1)
            .min(chrom_length);
        Self { region, anchor_point }
    }

    pub(super) fn query_region(&self) -> &GenomeRegion {
        &self.region
    }
}

impl LocalizedModCounts {
    fn add_bedmethyl_record(
        &mut self,
        bed_methyl_line: &BedMethylLine,
        anchor_point: u64,
    ) {
        let pos = bed_methyl_line.start() as i64;
        let offset = pos.saturating_sub(anchor_point as i64);
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
        let (left, right, has_single_offset) = match xs.iter().minmax() {
            MinMaxResult::NoElements => {
                bail!("cannot create localize chart: no offsets available")
            }
            MinMaxResult::OneElement(x) => {
                let x = *x;
                (x.saturating_sub(1), x.saturating_add(1), true)
            }
            MinMaxResult::MinMax(x, y) if x == y => {
                let x = *x;
                (x.saturating_sub(1), x.saturating_add(1), true)
            }
            MinMaxResult::MinMax(x, y) => (*x, *y, false),
        };
        let x_axis = Axis::new()
            .type_(AxisType::Value)
            .min(left)
            .max(right)
            .name("offset");
        let x_axis =
            if has_single_offset { x_axis.min_interval(1) } else { x_axis };
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
            .x_axis(x_axis)
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
                            info.percent_modified() as f64,
                        )),
                    ]))
                })
                .collect::<Vec<DataPoint>>();
            let is_singleton = dat.len() == 1;
            let line = Line::new()
                .name(format!("{mod_code}"))
                .data(dat)
                .line_style(LineStyle::new().width(1.5));
            let line = if is_singleton {
                line.symbol(Symbol::Circle).show_symbol(true)
            } else {
                line.symbol(Symbol::None)
            };
            chart = chart.series(line);
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

impl FocusRegion {
    pub(super) fn into_localized_mod_counts(
        self,
        index: &HtsTabixHandler<BedMethylLine>,
        min_coverage: u64,
        strand_rule: Option<StrandRule>,
        stranded_features: Option<StrandedFeatures>,
        io_threads: usize,
    ) -> anyhow::Result<LocalizedModCounts> {
        let Self { region, anchor_point } = self;
        let bedmethyl_records = index.fetch_region(
            &region.chrom,
            &(region.start..region.end),
            strand_rule.unwrap_or(region.strand),
            io_threads,
        )?;
        let loc_counts = bedmethyl_records
            .into_par_iter()
            .filter(|bm| bm.valid_coverage >= min_coverage)
            .filter(|bm| {
                stranded_features
                    .map(|f| {
                        let overlaps = region.strand.overlaps(&bm.strand);
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

#[cfg(test)]
mod tests {
    use super::FocusRegion;
    use crate::util::{GenomeRegion, StrandRule};

    #[test]
    fn focus_region_uses_inclusive_window_and_preserves_anchor() {
        let cases = [
            (10, 11, 0, 100, 10, 11, 10),
            (10, 14, 2, 100, 10, 15, 12),
            (0, 1, 2, 100, 0, 3, 0),
            (99, 100, 2, 100, 97, 100, 99),
        ];

        for (
            feature_start,
            feature_end,
            window,
            chrom_length,
            expected_start,
            expected_end,
            expected_anchor,
        ) in cases
        {
            let feature = GenomeRegion::new(
                "chr1".to_string(),
                feature_start,
                feature_end,
                StrandRule::Both,
                None,
            );
            let focus =
                FocusRegion::from_feature(feature, window, chrom_length);

            assert_eq!(focus.region.start, expected_start);
            assert_eq!(focus.region.end, expected_end);
            assert_eq!(focus.anchor_point, expected_anchor);
        }
    }
    use serde_json::Value;

    use super::LocalizedModCounts;
    use crate::mod_base_code::{
        ModCodeRepr, HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE,
    };
    use crate::util::ModPositionInfo;

    fn counts(
        series: impl IntoIterator<Item = (ModCodeRepr, Vec<(i64, u64, u64)>)>,
    ) -> LocalizedModCounts {
        let offsets = series
            .into_iter()
            .map(|(code, values)| {
                let values = values
                    .into_iter()
                    .map(|(offset, n_valid, n_mod)| {
                        (offset, ModPositionInfo::new(n_valid, n_mod))
                    })
                    .collect();
                (code, values)
            })
            .collect();
        LocalizedModCounts { offsets }
    }

    fn chart_json(counts: &LocalizedModCounts) -> Value {
        let html = counts.get_plot(None).unwrap();
        let raw_chart = html
            .split_once("var option = ")
            .unwrap()
            .1
            .split_once(";\n          chart.setOption")
            .unwrap()
            .0;
        serde_json::from_str(raw_chart).unwrap()
    }

    fn series_named<'a>(chart: &'a Value, name: &str) -> &'a Value {
        chart["series"]
            .as_array()
            .unwrap()
            .iter()
            .find(|series| series["name"] == name)
            .unwrap()
    }

    #[test]
    fn plot_uses_percent_modified() {
        let counts = counts([(METHYL_CYTOSINE, vec![(-1, 4, 1), (1, 4, 3)])]);
        let chart = chart_json(&counts);
        let data = chart["series"][0]["data"].as_array().unwrap();

        assert_eq!(chart["yAxis"][0]["name"], "percent modified");
        assert_eq!(data[0][1].as_f64(), Some(25.0));
        assert_eq!(data[1][1].as_f64(), Some(75.0));
    }

    #[test]
    fn plot_one_point_has_padded_axis_and_visible_symbol() {
        let counts = counts([(METHYL_CYTOSINE, vec![(0, 4, 1)])]);
        let chart = chart_json(&counts);
        let x_axis = &chart["xAxis"][0];
        let series = series_named(&chart, "m");

        assert_eq!(x_axis["min"].as_i64(), Some(-1));
        assert_eq!(x_axis["max"].as_i64(), Some(1));
        assert_eq!(x_axis["minInterval"].as_f64(), Some(1.0));
        assert_eq!(series["symbol"], "circle");
        assert_eq!(series["showSymbol"], true);
        assert_eq!(series["data"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn plot_multiple_singleton_codes_share_padded_axis() {
        let counts = counts([
            (METHYL_CYTOSINE, vec![(7, 4, 1)]),
            (HYDROXY_METHYL_CYTOSINE, vec![(7, 5, 2)]),
        ]);
        let chart = chart_json(&counts);
        let x_axis = &chart["xAxis"][0];

        assert_eq!(x_axis["min"].as_i64(), Some(6));
        assert_eq!(x_axis["max"].as_i64(), Some(8));
        assert_eq!(x_axis["minInterval"].as_f64(), Some(1.0));
        for name in ["m", "h"] {
            let series = series_named(&chart, name);
            assert_eq!(series["symbol"], "circle");
            assert_eq!(series["showSymbol"], true);
        }
    }

    #[test]
    fn plot_singleton_symbol_does_not_change_multi_point_series() {
        let counts = counts([
            (METHYL_CYTOSINE, vec![(-2, 4, 1), (2, 4, 3)]),
            (HYDROXY_METHYL_CYTOSINE, vec![(0, 5, 2)]),
        ]);
        let chart = chart_json(&counts);
        let multi = series_named(&chart, "m");
        let singleton = series_named(&chart, "h");

        assert_eq!(multi["symbol"], "none");
        assert!(multi.get("showSymbol").is_none());
        assert_eq!(singleton["symbol"], "circle");
        assert_eq!(singleton["showSymbol"], true);
    }

    #[test]
    fn plot_preserves_multi_offset_bounds() {
        let counts = counts([(METHYL_CYTOSINE, vec![(-5, 4, 1), (9, 4, 3)])]);
        let chart = chart_json(&counts);
        let x_axis = &chart["xAxis"][0];

        assert_eq!(x_axis["min"].as_i64(), Some(-5));
        assert_eq!(x_axis["max"].as_i64(), Some(9));
        assert!(x_axis.get("minInterval").is_none());
    }

    #[test]
    fn plot_empty_counts_returns_clear_error() {
        let error = LocalizedModCounts::default().get_plot(None).unwrap_err();
        assert!(error.to_string().contains("no offsets"));
    }
}
