use std::fmt::{Display, Formatter};
use std::ops::Range;

use anyhow::bail;
use log_once::debug_once;

const STATE_NUM: usize = 2usize;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
#[repr(usize)]
pub(crate) enum States {
    Same = 0,
    Different = 1,
}

impl Into<usize> for States {
    fn into(self) -> usize {
        self as usize
    }
}

impl From<usize> for States {
    fn from(value: usize) -> Self {
        match value {
            0 => Self::Same,
            1 => Self::Different,
            _ => unreachable!("invalid {value}"),
        }
    }
}

impl Display for States {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let label = match self {
            States::Same => "same",
            States::Different => "different",
        };
        write!(f, "{label}")
    }
}

#[derive(Copy, Clone, Debug)]
struct DpCell {
    inner: [f64; STATE_NUM], // todo make the state num const generic
}

impl DpCell {
    #[inline]
    fn new_full(val: f64) -> Self {
        Self { inner: [val; STATE_NUM] }
    }

    fn new_empty() -> Self {
        Self::new_full(f64::NEG_INFINITY)
    }

    // fn total_probability(&self) -> f64 {
    //     rv::misc::logsumexp(&self.inner)
    // }

    fn get_value(&self, state: States) -> f64 {
        self.inner[state as usize]
    }

    fn get_value_mut(&mut self, state: States) -> &mut f64 {
        &mut self.inner[state as usize]
    }

    fn set_value(&mut self, state: States, value: f64) {
        assert!(
            value.is_finite() && !value.is_nan(),
            "cannot set {value} to state {state:?}"
        );
        self.inner[state as usize] = value;
    }

    fn argmax(&self) -> States {
        self.inner
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| States::from(i))
            .unwrap()
    }
}

#[derive(Debug)]
struct PointerCell {
    inner: [Option<States>; STATE_NUM],
}

impl PointerCell {
    fn empty() -> Self {
        Self { inner: [None; STATE_NUM] }
    }

    fn get_value(&self, state: States) -> Option<States> {
        self.inner[state as usize]
    }

    fn set_value(&mut self, state: States, value: States) {
        self.inner[state as usize] = Some(value);
    }
}

pub(crate) struct HmmModel {
    same_to_same: f64,
    // diff_to_diff: f64,
    same_to_diff: f64,
    // diff_to_same: f64,
    dmr_prior: f64,

    same_state_factor: f64,
    diff_state_factor: f64,
    significance_factor: f64,

    linear_proj: bool,
    projection: Projection,
}

impl HmmModel {
    fn prob_to_factor(fpr: f64) -> anyhow::Result<f64> {
        if fpr < 0f64 {
            bail!("fpr cannot be less than 0")
        } else if fpr >= 1.0 {
            bail!("fpr cannot be >= 1.0")
        } else {
            Ok((fpr / (1f64 - fpr)).ln())
        }
    }

    pub(crate) fn new(
        dmr_prior: f64,
        diff_stay: f64,
        same_state_factor: f64,
        diff_state_factor: f64,
        significance_factor: f64,
        decay_distance: u32,
        linear_proj: bool,
    ) -> anyhow::Result<Self> {
        let same_to_diff = dmr_prior.ln();
        let same_to_same = (1f64 - dmr_prior).ln();
        // let diff_to_diff = diff_stay.ln();
        // let diff_to_same = (1f64 - diff_stay).ln();

        let projection = Projection::new(decay_distance, diff_stay, dmr_prior)?;
        let significance_factor = Self::prob_to_factor(significance_factor)?;

        Ok(Self {
            same_to_same,
            same_to_diff,
            same_state_factor,
            dmr_prior,
            diff_state_factor,
            significance_factor,
            linear_proj,
            projection,
        })
    }

    pub(crate) fn viterbi_path(
        &self,
        scores: &[f64],
        positions: &[u64],
    ) -> Vec<States> {
        // P_s = e^(-score)
        let probs = scores
            .iter()
            .map(|&x| if x < 0f64 { 0f64 } else { x })
            .map(|x| (-1f64 * x).exp())
            .collect::<Vec<f64>>();

        let transitions =
            positions.windows(2).fold(vec![self.dmr_prior], |mut agg, wind| {
                assert_eq!(wind.len(), 2);
                assert!(wind[1] > wind[0]);
                let gap = (wind[1] - wind[0]) as f64;
                assert!(gap > 0f64, "gap should be greater than zero");
                let p_diff_to_diff = if self.linear_proj {
                    self.projection.linear_project_prob(gap)
                } else {
                    self.projection.ln_project_prob(gap)
                };
                agg.push(p_diff_to_diff);
                agg
            });
        assert_eq!(probs.len(), transitions.len());
        let (dp_matrix, pointers) = self.viterbi_forward(&probs, &transitions);
        let path = self.viterbi_decode(&dp_matrix, &pointers);
        assert_eq!(path.len(), scores.len() - 1);
        path
    }

    fn viterbi_decode(
        &self,
        dp_matrix: &[DpCell],
        pointers: &[PointerCell],
    ) -> Vec<States> {
        let final_state = dp_matrix.last().unwrap().argmax();
        // dbg!(final_state);
        let mut path = vec![final_state];
        let mut curr_pointer =
            pointers.last().unwrap().get_value(final_state).unwrap();
        for pointers in pointers.iter().rev().skip(1) {
            let pointer = pointers.get_value(curr_pointer);
            if let Some(pointer) = pointer {
                path.push(pointer);
                curr_pointer = pointer;
            } else {
                break;
            }
        }

        path.pop();
        path.reverse();
        path
    }

    fn viterbi_forward(
        &self,
        scores: &[f64],
        transitions: &[f64],
    ) -> (Vec<DpCell>, Vec<PointerCell>) {
        let first_cell = {
            let mut first_cell = DpCell::new_full(0f64);
            self.initialize_start_end_cell(&mut first_cell);
            first_cell
        };
        let first_pointers = PointerCell::empty();
        let (mut dp_matrix, pointers, last_cell) =
            scores.iter().zip(transitions).enumerate().fold(
                (Vec::new(), vec![first_pointers], first_cell),
                |(mut cells, mut pointers, prev_cell), (i, (x, t))| {
                    let mut next_cell = DpCell::new_empty();
                    let mut pointer_cell = PointerCell::empty();
                    self.forward(
                        &prev_cell,
                        &mut next_cell,
                        &mut pointer_cell,
                        *t,
                        *x,
                        i,
                    );
                    cells.push(prev_cell);
                    pointers.push(pointer_cell);
                    (cells, pointers, next_cell)
                },
            );
        dp_matrix.push(last_cell);
        assert_eq!(dp_matrix.len(), pointers.len());
        assert_eq!(dp_matrix.len(), scores.len() + 1);
        (dp_matrix, pointers)
    }

    #[inline]
    fn emission_probs(&self, p: f64, state: States) -> f64 {
        let p = if p == 0f64 {
            debug_once!("encountered 0 prob");
            1e-5
        } else {
            p
        };
        assert!(p <= 1f64, "p {p} cannot be greater than 1");
        let (factor, p) = match state {
            States::Same => (self.same_state_factor, p.ln()),
            States::Different => {
                (self.diff_state_factor, (1f64 - p + 1e-5).ln())
            }
        };
        let p = p - self.significance_factor;
        factor * p
    }

    fn forward(
        &self,
        prev_cell: &DpCell,
        current_cell: &mut DpCell,
        pointers: &mut PointerCell,
        p_diff2diff: f64,
        score: f64,
        _idx: usize,
    ) {
        // todo make the naming convention here less terrible!
        // emission probs
        let e_diff = self.emission_probs(score, States::Different);
        let e_same = self.emission_probs(score, States::Same);
        // "dynamic" transition probs
        assert!(p_diff2diff > 0f64, "p_diff2diff should not be zero");
        assert!(
            p_diff2diff < 1.0,
            "p_diff2diff should be less than zero {p_diff2diff}"
        );
        let lnp_diff2diff = p_diff2diff.ln();
        let lnp_diff_to_same = (1f64 - p_diff2diff).ln();
        // previous state
        let p_same = prev_cell.get_value(States::Same);
        let p_diff = prev_cell.get_value(States::Different);

        Self::check_emission_prob(e_diff, "e_d");
        Self::check_emission_prob(e_same, "e_s");
        Self::check_emission_prob(p_diff, "p_d");
        Self::check_emission_prob(p_same, "p_s");
        Self::check_emission_prob(p_diff2diff, "p_diff2diff");
        Self::check_emission_prob(lnp_diff2diff, "lnp_diff2diff");
        Self::check_emission_prob(lnp_diff_to_same, "lnp_diff_to_same");

        // Same-state
        let same2same = p_same + self.same_to_same;
        let diff2same = p_diff + lnp_diff_to_same;
        let (current_same, same_pointer) =
            [(same2same, States::Same), (diff2same, States::Different)]
                .into_iter()
                .max_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap())
                .unwrap();

        // Diff-state
        let diff2diff = p_diff + lnp_diff2diff;
        let same2diff = p_same + self.same_to_diff;

        let (current_diff, diff_pointer) =
            [(diff2diff, States::Different), (same2diff, States::Same)]
                .into_iter()
                .max_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap())
                .unwrap();

        Self::check_emission_prob(current_diff, "current_diff");
        Self::check_emission_prob(current_same, "current_same");

        current_cell.set_value(States::Same, current_same + e_same);
        current_cell.set_value(States::Different, current_diff + e_diff);
        pointers.set_value(States::Same, same_pointer);
        pointers.set_value(States::Different, diff_pointer);
    }

    // todo make this a compile time no-op
    #[inline(always)]
    fn check_emission_prob(x: f64, which: &str) {
        assert!(x.is_finite(), "{which} is not finite {x}");
        assert!(!x.is_nan(), "{which} is NaN {x}");
    }

    fn initialize_start_end_cell(&self, cell: &mut DpCell) {
        *cell.get_value_mut(States::Same) = self.same_to_same;
        *cell.get_value_mut(States::Different) = self.same_to_diff;
    }
}

struct Projection {
    prob_range: Range<f64>,
    distance_range: Range<f64>,
    prob_span: f64,
    ratio: f64,
}

impl Projection {
    fn new(
        max_distance: u32,
        max_diff_stay: f64,
        dmr_prob: f64,
    ) -> anyhow::Result<Self> {
        if max_diff_stay <= dmr_prob {
            bail!("max_diff_stay must be > switch_prob")
        }
        let low = 1f64 - max_diff_stay;
        let high = 1f64 - dmr_prob;

        let prob_range = low..high;
        let max_distance = max_distance as f64;
        let distance_range = 2f64..max_distance;
        let prob_span = prob_range.end - prob_range.start;
        let ratio = prob_span / (distance_range.end - distance_range.start);

        Ok(Self { prob_range, distance_range, prob_span, ratio })
    }

    #[inline]
    fn clamp_value(&self, x: f64) -> f64 {
        if x > self.distance_range.end {
            self.distance_range.end
        } else {
            x
        }
    }

    fn linear_project_prob(&self, x: f64) -> f64 {
        let x = self.clamp_value(x);
        let adjusted = ((x - self.distance_range.start) * self.ratio)
            + self.prob_range.start;

        1f64 - adjusted
    }

    fn ln_project_prob(&self, x: f64) -> f64 {
        if x == 1f64 {
            return 1f64 - self.prob_range.start;
        }
        let x = self.clamp_value(x);
        let ln_ratio =
            self.distance_range.end.ln() - self.distance_range.start.ln();
        let adjusted = ((x.ln() - self.distance_range.start.ln()) / ln_ratio)
            * (self.prob_span)
            + self.prob_range.start;
        let prob = 1f64 - adjusted;
        if prob > 1.0 {
            panic!(
                "prob should not be >1 x: {x}, prob: {prob}, adjusted \
                 {adjusted}"
            )
        }
        prob
    }
}

#[cfg(test)]
mod hmm_tests {
    use crate::hmm::HmmModel;

    #[test]
    fn test_prob_to_factor() {
        let sig_fact = 0.01;
        let fact = HmmModel::prob_to_factor(sig_fact).unwrap();
        dbg!(fact);
    }
}
