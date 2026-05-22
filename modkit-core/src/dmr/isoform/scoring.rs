use statrs::distribution::{ChiSquared, ContinuousCDF};
use statrs::function::gamma::ln_gamma;

// #[derive(Debug, Clone)]
// pub(super) struct IsoformMethylationCounts {
//     /// Counts for categories such as:
//     /// [methylated_a, methylated_b, unmethylated]
//     pub counts: Vec<u64>,
// }

#[derive(Debug, Clone)]
pub(super) struct DirichletMultinomialTestResult {
    pub lrt: f64,
    pub _df: f64,
    pub p_value: f64,
    pub pooled_counts: Vec<u32>,
    // pub isoform_proportions: Vec<Vec<f64>>,
    pub isoform_counts: Vec<Vec<u32>>,
}

/// Clamp probabilities away from exact 0 and 1.
fn clamp_prob(p: f64) -> f64 {
    p.clamp(1e-9, 1.0 - 1e-9)
}

/// Convert category probabilities + overdispersion rho
/// into Dirichlet concentration parameters.
///
/// rho is analogous to beta-binomial overdispersion:
///
///   rho ≈ 0       -> close to multinomial
///   larger rho    -> more overdispersion
///
/// alpha_sum = (1 / rho) - 1
///
/// alpha_k = p_k * alpha_sum
fn dirichlet_params_from_probs(probs: &[f64], rho: f64) -> Vec<f64> {
    let rho = rho.clamp(1e-9, 1.0 - 1e-9);
    let alpha_sum = (1.0 / rho) - 1.0;

    probs.iter().map(|&p| clamp_prob(p) * alpha_sum).collect()
}

/// Normalize counts into probabilities, with pseudocount smoothing.
pub(super) fn proportions_from_counts(
    counts: &[u32],
    pseudocount: f64,
) -> Vec<f64> {
    let k = counts.len() as f64;
    let total: f64 =
        counts.iter().map(|&x| x as f64).sum::<f64>() + pseudocount * k;

    counts.iter().map(|&x| ((x as f64) + pseudocount) / total).collect()
}

/// Log PMF of the Dirichlet-multinomial distribution.
///
/// x ~ DirichletMultinomial(n, alpha)
///
/// log P(x | alpha)
fn log_dirichlet_multinomial_pmf(counts: &[u32], alpha: &[f64]) -> f64 {
    assert_eq!(
        counts.len(),
        alpha.len(),
        "counts and alpha must have the same number of categories"
    );

    let n: u32 = counts.iter().sum();

    let n_f = n as f64;
    let alpha_sum: f64 = alpha.iter().sum();

    let log_multinomial_coeff = ln_gamma(n_f + 1.0)
        - counts.iter().map(|&x| ln_gamma(x as f64 + 1.0)).sum::<f64>();

    let log_dirichlet_ratio = ln_gamma(alpha_sum) - ln_gamma(n_f + alpha_sum)
        + counts
            .iter()
            .zip(alpha.iter())
            .map(|(&x, &a)| ln_gamma(x as f64 + a) - ln_gamma(a))
            .sum::<f64>();

    log_multinomial_coeff + log_dirichlet_ratio
}

/// Dirichlet-multinomial likelihood-ratio test across isoforms.
///
/// Each isoform has K mutually exclusive categories.
///
/// Example categories:
///   [m6A, m5C, unmodified]
///
/// Null:
///   all isoforms share the same category composition
///
/// Alternative:
///   each isoform has its own category composition
///
/// rho:
///   fixed overdispersion parameter
///
/// pseudocount:
///   smoothing for estimating proportions.
///   Use 0.5 or 1.0 for stability.
pub(super) fn dirichlet_multinomial_lrt(
    isoforms: &[Vec<u32>],
    rho: f64,
    pseudocount: f64,
) -> Option<DirichletMultinomialTestResult> {
    if isoforms.len() < 2 {
        return None;
    }

    let num_categories = isoforms.first()?.len();

    if num_categories < 2 {
        return None;
    }

    if isoforms.iter().any(|x| x.len() != num_categories) {
        return None;
    }

    let usable_isoforms: Vec<_> =
        isoforms.iter().filter(|iso| iso.iter().sum::<u32>() > 0).collect();

    if usable_isoforms.len() < 2 {
        return None;
    }

    let mut pooled_counts = vec![0_u32; num_categories];

    for iso in &usable_isoforms {
        for (j, &x) in iso.iter().enumerate() {
            pooled_counts[j] += x;
        }
    }

    let pooled_total: u32 = pooled_counts.iter().sum();

    if pooled_total == 0 {
        return None;
    }

    let pooled_proportions =
        proportions_from_counts(&pooled_counts, pseudocount);
    let alpha_null = dirichlet_params_from_probs(&pooled_proportions, rho);

    let mut log_l_null = 0.0;
    let mut log_l_alt = 0.0;
    let mut isoform_proportions = Vec::with_capacity(usable_isoforms.len());
    let mut isoform_counts = Vec::with_capacity(usable_isoforms.len());

    for iso in usable_isoforms {
        let props = proportions_from_counts(&iso, pseudocount);
        let alpha_alt = dirichlet_params_from_probs(&props, rho);

        log_l_null += log_dirichlet_multinomial_pmf(&iso, &alpha_null);
        log_l_alt += log_dirichlet_multinomial_pmf(&iso, &alpha_alt);

        isoform_proportions.push(props);
        isoform_counts.push((*iso).clone());
    }

    let lrt = 2.0 * (log_l_alt - log_l_null).max(0.0);

    let num_isoforms = isoform_proportions.len();

    let df = ((num_isoforms - 1) * (num_categories - 1)) as f64;

    let chi = ChiSquared::new(df).ok()?;
    let p_value = 1.0 - chi.cdf(lrt);

    Some(DirichletMultinomialTestResult {
        lrt,
        _df: df,
        p_value,
        pooled_counts,
        isoform_counts,
    })
}
