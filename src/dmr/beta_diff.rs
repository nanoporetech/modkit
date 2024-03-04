use crate::dmr::llr_model::AggregatedCounts;
use anyhow::{anyhow, bail};
use log::info;
use rv::dist::Beta;
use rv::misc::gauss_legendre_quadrature;
use statrs::function::beta::ln_beta;
use std::fmt::{Display, Formatter};

const LOWER: f64 = 1e-5;
const UPPER: f64 = 1f64 - LOWER;

fn appell_f1_stable(x: f64, y: f64, a: f64, b1: f64, b2: f64, c: f64) -> f64 {
    let f = |u: f64| -> f64 {
        let numer = (a - 1f64) * u.ln() + (-a + c - 1f64) * (1f64 - u).ln();
        let denom = b1 * (1f64 - u * x).ln() + b2 * (1f64 - y * u).ln();
        (numer - denom).exp()
    };
    let correction = ln_beta(a, c - a);
    // let val = quad5(f, LOWER, UPPER);
    let val = gauss_legendre_quadrature(f, 16, (LOWER, UPPER));
    let val = val.ln();
    val - correction
}

#[derive(Copy, Clone)]
pub struct Counts {
    pub n_mod: usize,
    pub coverage: usize,
    pub frac_modified: f64,
}

impl Counts {
    pub fn new(n_mod: usize, coverage: usize) -> anyhow::Result<Self> {
        if n_mod > coverage {
            bail!("n_mod cannot be > coverage")
        } else {
            let frac_modified = n_mod as f64 / coverage as f64;
            Ok(Self { n_mod, coverage, frac_modified })
        }
    }

    fn empirical_effect_size(&self, other: &Self) -> f64 {
        self.frac_modified - other.frac_modified
    }

    fn resize(&self, max_coverage: usize) -> Self {
        if self.coverage > max_coverage {
            let n_mod = (self.frac_modified * max_coverage as f64).round();
            let frac_modified = n_mod / max_coverage as f64;
            Self {
                n_mod: n_mod as usize,
                coverage: max_coverage,
                frac_modified,
            }
        } else {
            Self {
                n_mod: self.n_mod,
                coverage: self.coverage,
                frac_modified: self.frac_modified,
            }
        }
    }

    fn n_canonical(&self) -> usize {
        self.coverage - self.n_mod
    }
}

impl TryFrom<&AggregatedCounts> for Counts {
    type Error = anyhow::Error;

    fn try_from(value: &AggregatedCounts) -> anyhow::Result<Self> {
        let n_mod = value.modified_counts();
        let coverage = value.total;
        Counts::new(n_mod, coverage)
    }
}

#[derive(Debug, Copy, Clone)]
pub struct BetaParams {
    alpha: f64,
    beta: f64,
}

impl BetaParams {
    pub(crate) fn new(alpha: f64, beta: f64) -> anyhow::Result<Self> {
        Beta::new(alpha, beta)
            .map(|dist| BetaParams { alpha: dist.alpha(), beta: dist.beta() })
            .map_err(|e| {
                anyhow!("invalid beta parameters {alpha}, {beta}, {e}")
            })
    }
}

impl Display for BetaParams {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let tmp = Beta::new_unchecked(self.alpha, self.beta);
        write!(f, "{tmp}")
    }
}

#[derive(Debug)]
pub struct EstimatedPMap {
    pub e_pmap: f64,
    pub effect_size: f64,
}

impl EstimatedPMap {
    fn new(e_pmap: f64, effect_size: f64) -> anyhow::Result<Self> {
        if e_pmap < 0f64 {
            bail!("e_pmap cannot be less than 0")
        } else {
            if e_pmap > 1f64 {
                Ok(Self::e_pmap_one(effect_size))
            } else {
                Ok(Self { e_pmap, effect_size })
            }
        }
    }

    fn e_pmap_one(effect_size: f64) -> Self {
        Self { e_pmap: 1.0f64, effect_size }
    }
}

const MAX_COV_ALLOWED: usize = 300usize;
pub struct PMapEstimator {
    max_coverages: [usize; 2],
    prior: BetaParams,
    rope: f64,
}

impl PMapEstimator {
    pub fn new(
        max_coverages: [usize; 2],
        a_num_reps: usize,
        b_num_reps: usize,
        prior: BetaParams,
        rope: f64,
        cap_coverages: bool,
    ) -> Self {
        let mut max_coverages = if cap_coverages {
            max_coverages
        } else {
            [max_coverages[0] * a_num_reps, max_coverages[1] * b_num_reps]
        };
        for x in max_coverages.iter_mut() {
            if *x > MAX_COV_ALLOWED {
                info!(
                    "calculated max coverage {x} is greater than maximum \
                     allowed ({MAX_COV_ALLOWED}), setting to {MAX_COV_ALLOWED}"
                );
                *x = MAX_COV_ALLOWED;
            }
        }

        Self { max_coverages, prior, rope }
    }

    fn calc_posterior_params(&self, counts: &Counts) -> BetaParams {
        let post_alpha = self.prior.alpha + counts.n_mod as f64;
        let post_beta = self.prior.beta + counts.n_canonical() as f64;
        // safe because we're adding
        BetaParams::new(post_alpha, post_beta).unwrap()
    }

    #[allow(non_snake_case)]
    pub fn calc_beta_diff(
        &self,
        d: f64,
        params1: &BetaParams,
        params2: &BetaParams,
    ) -> anyhow::Result<f64> {
        let ln_A = ln_beta(params1.alpha, params1.beta)
            + ln_beta(params2.alpha, params2.beta);
        if d.abs() < self.rope {
            if (params1.alpha + params2.alpha < 1f64)
                || (params1.beta + params2.beta < 1f64)
            {
                bail!(
                    "alpha1 + alpha2 <= 1 or beta1 + beta2 <= 1, params1 \
                     {params1:?}, params2 {params2:?}"
                )
            }
            let ln_p = ln_beta(
                params1.alpha + params2.alpha - 1f64,
                params1.beta + params2.beta - 1f64,
            ) - ln_A;
            Ok(ln_p)
        } else if d > 0f64 {
            let x = 1f64 - d;
            let y = 1f64 - d.powi(2);
            let a = params1.beta;
            let b1 =
                params1.alpha + params1.beta + params2.alpha + params2.beta
                    - 2f64;
            let b2 = 1f64 - params1.alpha;
            let c = params2.alpha + params1.beta;
            let f1 = appell_f1_stable(x, y, a, b1, b2, c);
            // dbg!(x, y, a, b1, b2, c, f1, ln_A);
            let ln_p = ln_beta(params2.alpha, params1.beta)
                + d.ln() * (params1.beta + params2.beta - 1f64)
                + (1f64 - d).ln() * (params2.alpha + params1.beta - 1f64)
                + f1
                - ln_A;
            Ok(ln_p)
        } else {
            let x = 1f64 - d.powi(2);
            let y = 1f64 + d;
            let a = params2.beta;
            let b1 = 1f64 - params2.alpha;
            let b2 =
                params1.alpha + params1.beta + params2.alpha + params2.beta
                    - 2f64;
            let c = params1.alpha + params2.beta;
            let f1 = appell_f1_stable(x, y, a, b1, b2, c);
            let ln_p = ln_beta(params1.alpha, params2.beta)
                + (-d).ln() * (params1.beta + params2.beta - 1f64)
                + (1f64 + d).ln() * (params1.alpha + params2.beta - 1f64)
                + f1
                - ln_A;
            Ok(ln_p)
        }
    }

    pub fn run(
        &self,
        a_counts: Counts,
        b_counts: Counts,
    ) -> anyhow::Result<EstimatedPMap> {
        // let a_counts: Counts = Counts::try_from(a_counts)?;
        // let b_counts: Counts = Counts::try_from(b_counts)?;
        let a_counts = a_counts.resize(self.max_coverages[0]);
        let b_counts = b_counts.resize(self.max_coverages[1]);
        let empirical_effect_size = a_counts.empirical_effect_size(&b_counts);
        if empirical_effect_size.abs() <= self.rope {
            // no difference
            return Ok(EstimatedPMap::e_pmap_one(empirical_effect_size));
        }
        // move away from 1 and -1
        // todo consider only doing this if the empirical effect size is E {-1,
        // 1}
        let adjusted_empirical_effect_size = if empirical_effect_size > 0f64 {
            empirical_effect_size - 0.005
        } else {
            empirical_effect_size + 0.005
        };

        let posterior_params_a = self.calc_posterior_params(&a_counts);
        let posterior_params_b = self.calc_posterior_params(&b_counts);

        let effect_prob = self.calc_beta_diff(
            adjusted_empirical_effect_size,
            &posterior_params_a,
            &posterior_params_b,
        )?;
        if effect_prob.exp() == 0.0f64 {
            Ok(EstimatedPMap::e_pmap_one(empirical_effect_size))
        } else {
            let null_prob = self.calc_beta_diff(
                0.0,
                &posterior_params_a,
                &posterior_params_b,
            )?;
            let ln_e_pmap = null_prob - effect_prob;
            let e_pmap = ln_e_pmap.exp();
            EstimatedPMap::new(e_pmap, empirical_effect_size)
        }
    }

    pub fn predict(
        &self,
        counts_a: &AggregatedCounts,
        counts_b: &AggregatedCounts,
    ) -> anyhow::Result<EstimatedPMap> {
        let mod_counts_a = Counts::try_from(counts_a)?;
        let mod_counts_b = Counts::try_from(counts_b)?;
        self.run(mod_counts_a, mod_counts_b)
    }
}

#[cfg(test)]
mod tests {
    use crate::dmr::beta_diff::{appell_f1_stable, LOWER, UPPER};
    use assert_approx_eq::assert_approx_eq;
    use rv::misc::gauss_legendre_quadrature;

    #[test]
    fn test_appell_f1_stable() {
        let answers = vec![
            3.4631730691211176,
            2.655223346206384,
            0.8708215438706287,
            0.4887961579016729,
            1.0,
        ];
        let xy: Vec<(f64, f64)> = vec![
            (0.9, 0f64),
            (0.7, 0.3),
            (-0.5, 0.2),
            (-0.9, -0.5),
            (0f64, 0f64),
        ];
        let a = 2f64;
        let b1 = 1f64;
        let b2 = 1f64;
        let c = 3f64;

        let correction = statrs::function::beta::beta(a, c - a);

        for (i, (x, y)) in xy.into_iter().enumerate() {
            let f = |u: f64| -> f64 {
                let numer = u.powf(a - 1f64) * (1f64 - u).powf(c - a - 1f64);
                let denom =
                    (1f64 - (u * x)).powf(b1) * (1f64 - (u * y)).powf(b2);
                numer / denom
            };
            let val = gauss_legendre_quadrature(f, 16, (LOWER, UPPER));
            let val = val / correction;
            let test_val = appell_f1_stable(x, y, a, b1, b2, c).exp();
            // println!("{x}, {y}, {val} {test_val}");
            assert_approx_eq!(val, test_val, 1e-4);
            assert_approx_eq!(test_val, answers[i], 1e-3);
        }
    }

    // #[test]
    // fn debug_test() {
    //     let estimator = PMapEstimator::new([10, 10], BetaParams::new(0.55,
    // 0.55), 0.05);     let p1 = BetaParams::new(13.55, 0.55);
    //     let p2 = BetaParams::new(9.55, 4.55);
    //     let ln_p = estimator.calc_beta_diff(0.3026923076923077, &p1, &p2);
    //     dbg!(ln_p);
    // }
}
