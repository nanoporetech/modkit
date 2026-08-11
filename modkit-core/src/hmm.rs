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
        assert_eq!(path.len(), scores.len());
        path
    }

    fn viterbi_decode(
        &self,
        dp_matrix: &[DpCell],
        pointers: &[PointerCell],
    ) -> Vec<States> {
        let final_state = dp_matrix.last().unwrap().argmax();
        let mut path = vec![final_state];
        let mut current_state = final_state;
        // The first pointer cell is the empty start cell, and the second
        // points back to the un-emitted start state. Decode only the emitted
        // states, from the final score back through the second score.
        for pointers in pointers.iter().skip(2).rev() {
            current_state = pointers.get_value(current_state).unwrap();
            path.push(current_state);
        }

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
    use crate::hmm::{HmmModel, States};

    fn test_model() -> HmmModel {
        HmmModel::new(0.1, 0.9, 0.3, -0.1, 0.01, 500, true).unwrap()
    }

    fn emission_score(model: &HmmModel, p: f64, state: States) -> f64 {
        let p = if p == 0f64 { 1e-5 } else { p };
        let (factor, log_probability) = match state {
            States::Same => (model.same_state_factor, p.ln()),
            States::Different => {
                (model.diff_state_factor, (1f64 - p + 1e-5).ln())
            }
        };
        factor * (log_probability - model.significance_factor)
    }

    fn transition_score(
        model: &HmmModel,
        previous: States,
        current: States,
        diff_stay: f64,
    ) -> f64 {
        match (previous, current) {
            (States::Same, States::Same) => model.same_to_same,
            (States::Same, States::Different) => model.same_to_diff,
            (States::Different, States::Same) => (1f64 - diff_stay).ln(),
            (States::Different, States::Different) => diff_stay.ln(),
        }
    }

    /// Exhaustively score the hidden start state and every emitted state.
    /// This deliberately does not use the dynamic-programming matrix or its
    /// back-pointers, so it independently checks both state order and length.
    fn brute_force_path(
        model: &HmmModel,
        scores: &[f64],
        positions: &[u64],
    ) -> Vec<States> {
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), positions.len());

        let probabilities = scores
            .iter()
            .map(|&score| (-score.max(0f64)).exp())
            .collect::<Vec<_>>();
        let diff_stays = positions.windows(2).fold(
            vec![model.dmr_prior],
            |mut transitions, window| {
                let gap = (window[1] - window[0]) as f64;
                transitions.push(if model.linear_proj {
                    model.projection.linear_project_prob(gap)
                } else {
                    model.projection.ln_project_prob(gap)
                });
                transitions
            },
        );

        let state_count = scores.len() + 1;
        let mut best: Option<(f64, Vec<States>)> = None;
        for encoded_path in 0..(1usize << state_count) {
            let states = (0..state_count)
                .map(|i| {
                    if encoded_path & (1usize << i) == 0 {
                        States::Same
                    } else {
                        States::Different
                    }
                })
                .collect::<Vec<_>>();
            let initial_score = match states[0] {
                States::Same => model.same_to_same,
                States::Different => model.same_to_diff,
            };
            let total_score = probabilities
                .iter()
                .zip(diff_stays.iter())
                .enumerate()
                .fold(initial_score, |total, (i, (&p, &diff_stay))| {
                    total
                        + transition_score(
                            model,
                            states[i],
                            states[i + 1],
                            diff_stay,
                        )
                        + emission_score(model, p, states[i + 1])
                });

            if best
                .as_ref()
                .map(|(best_score, _)| total_score > *best_score)
                .unwrap_or(true)
            {
                best = Some((total_score, states[1..].to_vec()));
            }
        }
        best.unwrap().1
    }

    #[test]
    fn test_prob_to_factor() {
        let sig_fact = 0.01;
        let fact = HmmModel::prob_to_factor(sig_fact).unwrap();
        dbg!(fact);
    }

    #[test]
    fn viterbi_path_matches_independent_oracle_for_tiny_sequences() {
        let model = test_model();
        let cases = [
            (vec![0.0], vec![10]),
            (vec![12.0], vec![10]),
            (vec![0.0, 12.0], vec![10, 20]),
            (vec![12.0, 0.0], vec![10, 20]),
            (vec![0.0, 0.0, 12.0, 12.0, 0.0], vec![10, 20, 30, 40, 50]),
        ];

        for (scores, positions) in cases {
            let expected = brute_force_path(&model, &scores, &positions);
            let actual = model.viterbi_path(&scores, &positions);
            assert_eq!(actual.len(), scores.len());
            assert_eq!(actual, expected, "scores: {scores:?}");
        }
    }

    #[test]
    fn viterbi_path_preserves_state_transitions_in_order() {
        let model = test_model();
        let scores = vec![0.0, 0.0, 12.0, 12.0, 0.0];
        let positions = vec![10, 20, 30, 40, 50];
        let expected = brute_force_path(&model, &scores, &positions);
        assert!(expected.windows(2).any(|states| states[0] != states[1]));

        assert_eq!(model.viterbi_path(&scores, &positions), expected);
    }
}
