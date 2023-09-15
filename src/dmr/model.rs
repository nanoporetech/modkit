use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::Itertools;
use log::{debug, info};
use ndarray::{Array2, Axis};
use rv::prelude::*;

use crate::dmr::DmrInterval;

#[derive(Debug)]
pub(crate) struct AggregatedCounts {
    mod_code_counts: HashMap<char, usize>,
    total: usize,
}

impl AggregatedCounts {
    pub(crate) fn try_new(
        mod_code_counts: HashMap<char, usize>,
        total: usize,
    ) -> anyhow::Result<Self> {
        let total_modification_counts = mod_code_counts.values().sum::<usize>();
        if total_modification_counts > total {
            bail!(
                "total modification counts cannot be greater than total counts"
            )
        }
        Ok(Self {
            mod_code_counts,
            total,
        })
    }

    fn get_canonical_counts(&self) -> usize {
        // safe because we check at creation, could be more careful if there
        // was a chance that &mut self was available.
        self.total - self.mod_code_counts.values().sum::<usize>()
    }

    fn combine(&self, other: &Self) -> Self {
        let total = self.total + other.total;
        let mut counts = self.mod_code_counts.clone();
        other.mod_code_counts.iter().for_each(|(mod_code, count)| {
            *counts.entry(*mod_code).or_insert(0) += *count;
        });

        Self {
            mod_code_counts: counts,
            total,
        }
    }

    fn categorical_trials(
        &self,
        mod_codes_to_index: &HashMap<char, usize>,
    ) -> anyhow::Result<Vec<usize>> {
        let mut trials = self
            .mod_code_counts
            .iter()
            .try_fold(Vec::new(), |mut acc, (code, count)| {
                let index = *mod_codes_to_index.get(code)?;
                let mut trials = vec![index; *count];
                acc.append(&mut trials);
                Some(acc)
            })
            .ok_or(anyhow!("failed to make categotical trials"))?;
        let canonical_count = self.get_canonical_counts();
        trials.append(&mut vec![0usize; canonical_count]);
        Ok(trials)
    }

    fn string_counts(&self) -> String {
        let csv = self
            .mod_code_counts
            .iter()
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
            .fold(String::new(), |mut acc, (code, count)| {
                acc.push_str(&format!("{}:{},", code, count));
                acc
            });
        csv.chars().into_iter().take(csv.len() - 1).collect()
    }

    fn string_percentages(&self) -> String {
        let csv = self
            .mod_code_counts
            .iter()
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
            .fold(String::new(), |mut acc, (code, count)| {
                let frac = *count as f32 / self.total as f32;
                acc.push_str(&format!("{}:{:.2},", code, frac * 100f32));
                acc
            });
        csv.chars().into_iter().take(csv.len() - 1).collect()
    }
}

impl Display for AggregatedCounts {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.string_counts())
    }
}

#[derive(Debug)]
pub(crate) struct ModificationCounts {
    start: u64,
    stop: u64,
    control_counts: AggregatedCounts,
    // control_counts: HashMap<char, usize>,
    // control_total: usize,
    exp_counts: AggregatedCounts,
    // exp_counts: HashMap<char, usize>,
    // exp_total: usize,
    interval: DmrInterval,
    pub(crate) score: f64,
}

impl ModificationCounts {
    pub(crate) fn new(
        start: u64,
        stop: u64,
        control_counts: AggregatedCounts,
        exp_counts: AggregatedCounts,
        interval: DmrInterval,
    ) -> anyhow::Result<Self> {
        let score = llk_ratio(&control_counts, &exp_counts)?;
        Ok(Self {
            start,
            stop,
            control_counts,
            exp_counts,
            interval,
            score,
        })
    }

    pub(super) fn to_row(&self) -> anyhow::Result<String> {
        let sep = '\t';
        // let mle_llk_ratio = mle_llk_ratio(&self)?;
        // dbg!(&self);
        // let chi2_stat = chi_squared_test(&self)?;

        let line = format!(
            "\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}\n\
        ",
            self.interval.chrom,
            self.start,
            self.stop,
            self.interval.name,
            self.score,
            self.control_counts.string_counts(),
            self.control_counts.total,
            self.exp_counts.string_counts(),
            self.exp_counts.total,
        );
        Ok(line)
    }
}

// fn multi_counts_to_trials(
//     counts: &AggregatedCounts,
//     // counts: &HashMap<char, usize>,
//     code_to_index: &HashMap<char, usize>,
//     // total: usize,
// ) -> anyhow::Result<Vec<usize>> {
//     let mut trials = counts.mod_code_counts.iter().fold(
//         Vec::new(),
//         |mut acc, (code, count)| {
//             let index = *code_to_index.get(code).unwrap();
//             let mut trials = vec![index; *count];
//             acc.append(&mut trials);
//             acc
//         },
//     );
//     let canonical_count = counts.canonical_counts();
//     // let remaining = total
//     //     .checked_sub(trials.len())
//     //     .ok_or_else(|| anyhow!("total less that total methyl counts"))?;
//     trials.append(&mut vec![0; canonical_count]);
//     Ok(trials)
// }

fn dirichlet_llk(
    counts: &AggregatedCounts,
    prior: &Dirichlet,
    mod_codes_to_index: &HashMap<char, usize>,
) -> anyhow::Result<f64> {
    // categorical outputs, die rolls, etc.
    let xs = counts.categorical_trials(&mod_codes_to_index)?;
    let control_data = DataOrSuffStat::Data(&xs);
    let control_posterior = prior.posterior(&control_data);
    Ok(control_posterior.ln_m(&control_data))
}

fn llk_dirichlet(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let mods_to_index = control_counts
        .mod_code_counts
        .keys()
        .chain(exp_counts.mod_code_counts.keys())
        .copied()
        .collect::<HashSet<char>>()
        .into_iter()
        .sorted_by(|a, b| a.cmp(b))
        .enumerate()
        .map(|(i, c)| (c, i + 1))
        .collect::<HashMap<char, usize>>();

    let k = mods_to_index.len() + 1;
    let prior = Dirichlet::jeffreys(k)?;
    let llk_control = dirichlet_llk(&control_counts, &prior, &mods_to_index)?;
    // let control_trials =
    //     multi_counts_to_trials(&control_counts, &mods_to_index)?;
    // let control_data = DataOrSuffStat::Data(&control_trials);
    // let control_posterior = prior.posterior(&control_data);
    // let llk_control = control_posterior.ln_m(&control_data);

    // let exp_trials = multi_counts_to_trials(
    //     &modification_counts.exp_counts,
    //     &mods_to_index,
    //     modification_counts.exp_total,
    // )?;
    // let exp_data = DataOrSuffStat::Data(&exp_trials);
    // let exp_posterior = prior.posterior(&exp_data);
    let llk_exp = dirichlet_llk(&exp_counts, &prior, &mods_to_index)?;
    // let llk_exp = exp_posterior.ln_m(&exp_data);

    let combined_counts = control_counts.combine(exp_counts);
    let llk_combined = dirichlet_llk(&combined_counts, &prior, &mods_to_index)?;
    //
    // let all_trials = control_trials
    //     .into_iter()
    //     .chain(exp_trials.into_iter())
    //     .collect::<Vec<usize>>();
    // let all_data = DataOrSuffStat::Data(&all_trials);
    // let all_post = prior.posterior(&all_data);
    // let llk_same = all_post.ln_m(&all_data);

    // debug!("control dist:{} exp dist:{} combined dist:{}",
    //     control_posterior.to_string(),
    //     exp_posterior.to_string(),
    //     all_post.to_string());

    Ok(llk_control + llk_exp - llk_combined)
}

fn counts_to_trials(count_methyl: usize, count_canonical: usize) -> Vec<bool> {
    let mut x = vec![true; count_methyl];
    let mut y = vec![false; count_canonical];
    x.append(&mut y);
    x
}

fn beta_llk(count_methyl: usize, count_canonical: usize) -> f64 {
    let trials = counts_to_trials(count_methyl, count_canonical);
    let data = BernoulliData::Data(&trials);
    let prior = Beta::jeffreys();
    let posterior = prior.posterior(&data);
    posterior.ln_m(&data)
}

fn llk_beta(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let all_mods = control_counts
        .mod_code_counts
        .keys()
        .copied()
        .chain(exp_counts.mod_code_counts.keys().copied())
        .collect::<HashSet<char>>();
    if all_mods.len() != 1 {
        bail!("should have exactly one modification to use beta llk")
    }
    let raw_mod_code = all_mods.into_iter().take(1).collect::<Vec<char>>()[0];

    let control_methyls = *control_counts
        .mod_code_counts
        .get(&raw_mod_code)
        .unwrap_or(&0);
    let control_canonicals = control_counts.get_canonical_counts();
    // let control_canonicals =
    //     modification_counts.control_total - control_methyls;

    let llk_control = beta_llk(control_methyls, control_canonicals);
    // let control_trials = counts_to_trials(control_methyls, control_canonicals);
    // let control_data = BernoulliData::Data(&control_trials);
    // let prior = Beta::jeffreys();
    // let control_posterior = prior.posterior(&control_data);
    // let llk_control = control_posterior.ln_m(&control_data);

    let exp_methyls =
        *exp_counts.mod_code_counts.get(&raw_mod_code).unwrap_or(&0);
    // assert!(exp_methyls <= modification_counts.exp_total);
    let exp_canonicals = exp_counts.get_canonical_counts();
    // let exp_canonicals = modification_counts.exp_total - exp_methyls;
    let llk_exp = beta_llk(exp_methyls, exp_canonicals);
    // let exp_trials = counts_to_trials(exp_methyls, exp_canonicals);
    // let exp_data = DataOrSuffStat::Data(&exp_trials);
    // let exp_posterior = prior.posterior(&exp_data);
    // let llk_exp = exp_posterior.ln_m(&exp_data);

    // let all_trials = counts_to_trials(
    //     exp_methyls + control_methyls,
    //     exp_canonicals + control_canonicals,
    // );
    // let all_data = DataOrSuffStat::Data(&all_trials);
    // let combined_posterior = prior.posterior(&all_data);
    // let llk_same = combined_posterior.ln_m(&all_data);
    let llk_same = beta_llk(
        exp_methyls + control_methyls,
        exp_canonicals + control_canonicals,
    );

    Ok(llk_control + llk_exp - llk_same)
}

pub(crate) fn llk_ratio(
    control_counts: &AggregatedCounts,
    exp_counts: &AggregatedCounts,
) -> anyhow::Result<f64> {
    let n_categories = std::cmp::max(
        control_counts.mod_code_counts.keys().len(),
        exp_counts.mod_code_counts.keys().len(),
    ) + 1; // plus 1 for canonical
    if n_categories < 2 {
        return Ok(0f64);
    }
    if n_categories == 2 {
        llk_beta(control_counts, exp_counts)
    } else {
        llk_dirichlet(control_counts, exp_counts)
    }
}

fn get_counts_each(
    obs: &HashMap<char, usize>,
    mods_to_index: &HashMap<char, usize>,
    total: usize,
) -> Vec<u32> {
    let mut counts = vec![0u32; mods_to_index.len() + 1];
    for (code, count) in obs.iter() {
        let idx = *mods_to_index.get(code).unwrap();
        counts[idx] = *count as u32;
    }
    let count_modified = obs.values().sum::<usize>();
    assert!(count_modified <= total);
    let count_unmodified = (total - count_modified) as u32;
    counts[0] = count_unmodified;
    counts
}

// pub(crate) fn chi_squared_test(
//     modification_counts: &ModificationCounts,
// ) -> anyhow::Result<f64> {
//     if modification_counts.control_counts.total == 0
//         || modification_counts.exp_counts.total == 0
//     {
//         info!("site has zero total");
//         return Ok(1f64);
//     }
//     let n_conditions = 2;
//     let mods_to_index = modification_counts
//         .control_counts
//         .mod_code_counts
//         .keys()
//         .chain(modification_counts.exp_counts.mod_code_counts.keys())
//         .copied()
//         .collect::<HashSet<char>>()
//         .into_iter()
//         .sorted_by(|a, b| a.cmp(b))
//         .enumerate()
//         .map(|(i, c)| (c, i + 1))
//         .collect::<HashMap<char, usize>>();
//     let n_cats = mods_to_index.len() + 1;
//     let obs = vec![
//         get_counts_each(
//             &modification_counts.control_counts.mod_code_counts,
//             &mods_to_index,
//             modification_counts.control_counts.total,
//         ),
//         get_counts_each(
//             &modification_counts.exp_counts.mod_code_counts,
//             &mods_to_index,
//             modification_counts.exp_counts.total,
//         ),
//     ];
//
//     let data = obs.iter().flatten().copied().collect::<Vec<u32>>();
//     let contingency = Array2::from_shape_vec((n_conditions, n_cats), data)
//         .map(|a| a.t().to_owned().mapv(|x| f64::from(x + 1)))?;
//
//     // let condition_totals = contingency.sum_axis(Axis(0)).mapv(|x| f64::from(x));
//     let category_totals = contingency.sum_axis(Axis(1)).mapv(|x| f64::from(x));
//     let total: f64 = category_totals.sum();
//     let category_proportions = &category_totals / total;
//     let stat =
//         contingency
//             .columns()
//             .into_iter()
//             .fold(0f64, |acc, cond_counts| {
//                 let cond_total = cond_counts.sum();
//                 let expected = cond_total * &category_proportions;
//                 let stat =
//                     (&cond_counts - &expected).mapv(|x| x.powi(2)) / expected;
//                 let stat = stat.sum();
//                 acc + stat
//             });
//     let df = (n_conditions - 1) * (n_cats - 1);
//     let chi_dist = ChiSquared::new(df as f64)?;
//     let p_val = chi_dist.sf(&stat);
//
//     Ok(p_val)
// }

// assumed that mod_counts does not contain the canonical counts
fn modified_counts_to_frequencies(
    modified_counts: &HashMap<char, usize>,
    mod_code_to_index: &HashMap<char, usize>,
    total: usize,
) -> Vec<f64> {
    let mut frequencies = vec![0f64; mod_code_to_index.len() + 1];
    for (code, count) in modified_counts.iter() {
        let freq = *count as f64 / total as f64;
        let idx = *mod_code_to_index.get(code).unwrap();
        frequencies[idx] = freq;
    }

    let remaining = modified_counts.values().sum::<usize>();
    assert!(remaining <= total);
    let count_unmodified = total - remaining;
    let unmodified_freq = count_unmodified as f64 / total as f64;
    frequencies[0] = unmodified_freq;

    frequencies.into_iter().map(|x| x + 0.001).collect()
}

// pub(crate) fn mle_llk_ratio(
//     modification_counts: &ModificationCounts,
// ) -> anyhow::Result<f64> {
//     if modification_counts.control_total == 0
//         || modification_counts.exp_total == 0
//     {
//         info!("site has zero total");
//         return Ok(0f64);
//     }
//     let mods_to_index = modification_counts
//         .control_counts
//         .keys()
//         .chain(modification_counts.exp_counts.keys())
//         .copied()
//         .collect::<HashSet<char>>()
//         .into_iter()
//         .sorted_by(|a, b| a.cmp(b))
//         .enumerate()
//         .map(|(i, c)| (c, i + 1))
//         .collect::<HashMap<char, usize>>();
//
//     let alphas_control = modified_counts_to_frequencies(
//         &modification_counts.control_counts.mod_code_counts,
//         &mods_to_index,
//         modification_counts.control_counts.total,
//     );
//     let alphas_exp = modified_counts_to_frequencies(
//         &modification_counts.exp_counts.mod_code_counts,
//         &mods_to_index,
//         modification_counts.exp_counts.total,
//     );
//     let combined_counts = modification_counts
//         .control_counts
//         .iter()
//         .chain(modification_counts.exp_counts.iter())
//         .fold(HashMap::new(), |mut acc, (mod_code, count)| {
//             *acc.entry(*mod_code).or_insert(0) += *count;
//             acc
//         });
//     let alphas_combined = modified_counts_to_frequencies(
//         &combined_counts,
//         &mods_to_index,
//         modification_counts.exp_total + modification_counts.control_total,
//     );
//
//     // let control_trials = multi_counts_to_trials(
//     //     &modification_counts.control_counts,
//     //     &mods_to_index,
//     //     modification_counts.control_total,
//     // )?;
//     let control_trials = modification_counts
//         .control_counts
//         .categorical_trials(&mods_to_index);
//     let control_data = CategoricalData::Data(&control_trials);
//     let lk_control = Dirichlet::new(alphas_control.clone())
//         .with_context(|| {
//             format!(
//                 "invalid control Dir params {:?}, counts: {:?} {:?}",
//                 &alphas_control,
//                 &modification_counts.control_counts,
//                 &modification_counts.control_total
//             )
//         })?
//         .ln_m(&control_data);
//
//     let exp_trials = multi_counts_to_trials(
//         &modification_counts.exp_counts,
//         &mods_to_index,
//         modification_counts.exp_total,
//     )?;
//     let exp_data = CategoricalData::Data(&exp_trials);
//     let lk_exp = Dirichlet::new(alphas_exp.clone())
//         .with_context(|| format!("invalid exp Dir params {:?}", &alphas_exp))?
//         .ln_m(&exp_data);
//
//     let combined_trials = multi_counts_to_trials(
//         &combined_counts,
//         &mods_to_index,
//         modification_counts.exp_total + modification_counts.control_total,
//     )?;
//     let combined_data = CategoricalData::Data(&combined_trials);
//     let lk_combined = Dirichlet::new(alphas_combined.clone())
//         .with_context(|| {
//             format!("invalid combined Dir params {:?}", &alphas_combined)
//         })?
//         .ln_m(&combined_data);
//
//     Ok(lk_control + lk_exp - lk_combined)
// }

#[cfg(test)]
mod dmr_model_tests {
    use std::collections::HashMap;

    use crate::dmr::model::{chi_squared_test, ModificationCounts};
    use crate::dmr::DmrInterval;
    use crate::position_filter::Iv;

    #[test]
    fn test_chi2_stat() {
        let interval = DmrInterval::new(
            Iv {
                start: 1,
                stop: 10,
                val: (),
            },
            "xxx".to_string(),
            "name".to_string(),
        );

        let control_counts =
            HashMap::<char, usize>::from([('h', 18usize), ('m', 0)]);
        let control_total = 1507;
        let exp_counts =
            HashMap::<char, usize>::from([('h', 3usize), ('m', 0)]);
        let exp_total = 1079;
        let mod_counts = ModificationCounts::new(
            1,
            10,
            control_counts,
            control_total,
            exp_counts,
            exp_total,
            interval,
        );

        let stat = chi_squared_test(&mod_counts).unwrap();
    }
}
