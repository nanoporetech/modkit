use crate::features::CountsFeatures;
use derive_new::new;
use log::debug;
use mod_kit::util::TAB;
use ndarray::{s, Array1};
use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, VecDeque};
use std::fmt::{Display, Formatter};
use std::ops::Div;

#[derive(new, Debug)]
pub(crate) struct OpenChromRecord {
    pub(crate) chrom: String,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) score: f32,
}

impl Display for OpenChromRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{TAB}{}{TAB}{}{TAB}{}",
            self.chrom, self.start, self.end, self.score
        )
    }
}

pub(crate) struct ProbMerger {
    preds: FxHashMap<String, BTreeMap<(u32, u32), f32>>,
    filter_threshold: f32,
}

impl ProbMerger {
    pub(crate) fn new(filter_threshold: f32) -> Self {
        Self { preds: FxHashMap::default(), filter_threshold }
    }

    pub(crate) fn add_preds(
        &mut self,
        counts: &[CountsFeatures<()>],
        probs: &[f32],
    ) {
        let iter = counts.iter().zip(probs);
        // .filter(|(_, p)| **p >= self.filter_threshold);
        for (feat, prob) in iter {
            let chrom_intervals = self
                .preds
                .entry(feat.chrom.to_owned())
                .or_insert_with(|| BTreeMap::new());
            let k = (feat.start, feat.end);
            chrom_intervals.insert(k, *prob);
        }
    }

    pub(crate) fn into_records(self) -> Vec<OpenChromRecord> {
        self.preds
            .into_iter()
            .flat_map(|(chrom, intervals)| {
                calc_base_level_probs(chrom, intervals, self.filter_threshold)
            })
            .collect()
    }

    pub(crate) fn size(&self) -> usize {
        self.preds.values().map(|x| x.len()).sum::<usize>()
    }
}

fn calc_base_level_probs(
    chrom: String,
    mut intervals: BTreeMap<(u32, u32), f32>,
    filter_threshold: f32,
) -> Vec<OpenChromRecord> {
    if intervals.is_empty() {
        return Vec::new();
    } else if intervals.len() == 1 {
        let ((st, en), p) = intervals.pop_first().unwrap();
        vec![OpenChromRecord::new(chrom.to_owned(), st, en, p)]
    } else {
        let ((leftmost, en), lp) = intervals.pop_first().unwrap();
        let ((st, rightmost), rp) = intervals.pop_last().unwrap();
        let width = rightmost.saturating_sub(leftmost) as usize;
        if width == 0 {
            debug!("zero/negative length width?");
            return Vec::new();
        }
        let mut mat = Array1::<f32>::zeros(width);
        let mut coverages = Array1::<f32>::zeros(width);
        // add the first one
        let first_interval_end = en.saturating_sub(leftmost) as usize;
        let slicer = s![0..first_interval_end];
        mat.slice_mut(&slicer).mapv_inplace(|x| x + lp);
        coverages.slice_mut(&slicer).mapv_inplace(|x| x + 1f32);
        // add the last one
        let last_interval_start = st.saturating_sub(leftmost) as usize;
        let slicer = s![last_interval_start..width];
        mat.slice_mut(&slicer).mapv_inplace(|x| x + rp);
        coverages.slice_mut(&slicer).mapv_inplace(|x| x + 1f32);
        // add the rest
        for ((st, en), p) in intervals {
            let st = st.saturating_sub(leftmost) as usize;
            let en = en.saturating_sub(leftmost) as usize;
            let slicer = s![st..en];
            mat.slice_mut(&slicer).mapv_inplace(|x| x + p);
            coverages.slice_mut(&slicer).mapv_inplace(|x| x + 1f32);
        }

        let vals =
            mat.div(coverages).to_vec().into_iter().collect::<VecDeque<f32>>();
        debug_assert!(!vals.is_empty());
        debug_assert_eq!(vals.len(), width);
        let trim = vals.iter().position(|p| !p.is_nan());
        if let Some(trim) = trim {
            let mut curr_prob = vals[trim];
            let mut curr_left = leftmost.saturating_add(trim as u32);
            let mut curr_right = curr_left.saturating_add(1);
            let mut records = Vec::<OpenChromRecord>::new();
            for prob in vals.iter().skip(trim.saturating_add(1)) {
                if *prob == curr_prob {
                    curr_right = curr_right.saturating_add(1);
                    continue;
                } else if prob.is_nan() {
                    if !curr_prob.is_nan() {
                        let record = OpenChromRecord::new(
                            chrom.to_owned(),
                            curr_left,
                            curr_right,
                            curr_prob,
                        );
                        records.push(record);
                    }
                    curr_left = curr_right;
                    curr_right = curr_right.saturating_add(1);
                    curr_prob = *prob;
                } else {
                    // changed into another run
                    if !curr_prob.is_nan() {
                        let record = OpenChromRecord::new(
                            chrom.to_owned(),
                            curr_left,
                            curr_right,
                            curr_prob,
                        );
                        records.push(record);
                    }
                    curr_left = curr_right;
                    curr_right = curr_right.saturating_add(1);
                    curr_prob = *prob;
                }
            }
            if !curr_prob.is_nan() {
                let record = OpenChromRecord::new(
                    chrom.to_owned(),
                    curr_left,
                    curr_right,
                    curr_prob,
                );
                records.push(record);
            }
            records.retain(|x| x.score >= filter_threshold);
            return records;
        } else {
            debug!("all NaN");
            return Vec::new();
        }
    }
}

#[cfg(test)]
mod prob_merger_tests {
    use crate::prob_merger::calc_base_level_probs;
    use std::collections::BTreeMap;

    #[test]
    fn test_prob_merger_merges() {
        let intervals = [
            ((210, 310), 0.9),
            ((215, 315), 0.9),
            ((220, 320), 0.8),
            ((230, 330), 0.7),
            ((400, 500), 0.9),
            ((405, 505), 0.8),
        ]
        .into_iter()
        .collect::<BTreeMap<(u32, u32), f32>>();
        let records = calc_base_level_probs("chr1".to_string(), intervals, 0.0);
        for record in records {
            println!("{record}");
        }
    }
}
