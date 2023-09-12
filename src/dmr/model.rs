use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::Itertools;
use log::{debug, info};
use ndarray::{Array2, Axis};
use rv::prelude::*;

use crate::dmr::DmrInterval;

#[derive(new, Debug)]
pub(crate) struct ModificationCounts {
    start: u64,
    stop: u64,
    control_counts: HashMap<char, usize>,
    control_total: usize,
    exp_counts: HashMap<char, usize>,
    exp_total: usize,
    interval: DmrInterval,
}

impl ModificationCounts {
    fn string_counts(mapping: &HashMap<char, usize>) -> String {
        let csv = mapping.iter().sorted_by(|(a, _), (b, _)| a.cmp(b)).fold(
            String::new(),
            |mut acc, (code, count)| {
                acc.push_str(&format!("{}:{},", code, count));
                acc
            },
        );
        csv.chars().into_iter().take(csv.len() - 1).collect()
    }

    pub(super) fn to_row(&self, score: f64) -> anyhow::Result<String> {
        let sep = '\t';
        let mle_llk_ratio = mle_llk_ratio(&self)?;
        // dbg!(&self);
        let chi2_stat = chi_squared_test(&self)?;

        let line = format!(
            "\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {}{sep}\
        {score}{sep}\
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
            Self::string_counts(&self.control_counts),
            self.control_total,
            Self::string_counts(&self.exp_counts),
            self.exp_total,
            mle_llk_ratio,
            chi2_stat,
        );
        Ok(line)
    }
}

fn multi_counts_to_trials(
    counts: &HashMap<char, usize>,
    code_to_index: &HashMap<char, usize>,
    total: usize,
) -> anyhow::Result<Vec<usize>> {
    let mut trials =
        counts.iter().fold(Vec::new(), |mut acc, (code, count)| {
            let index = *code_to_index.get(code).unwrap();
            let mut trials = vec![index; *count];
            acc.append(&mut trials);
            acc
        });
    let remaining = total
        .checked_sub(trials.len())
        .ok_or_else(|| anyhow!("total less that total methyl counts"))?;
    trials.append(&mut vec![0; remaining]);
    Ok(trials)
}

fn llk_dirichlet(
    modification_counts: &ModificationCounts,
) -> anyhow::Result<f64> {
    let mods_to_index = modification_counts
        .control_counts
        .keys()
        .chain(modification_counts.exp_counts.keys())
        .copied()
        .collect::<HashSet<char>>()
        .into_iter()
        .sorted_by(|a, b| a.cmp(b))
        .enumerate()
        .map(|(i, c)| (c, i + 1))
        .collect::<HashMap<char, usize>>();

    let prior = Dirichlet::jeffreys(mods_to_index.len() + 1)?;
    let control_trials = multi_counts_to_trials(
        &modification_counts.control_counts,
        &mods_to_index,
        modification_counts.control_total,
    )?;
    let control_data = DataOrSuffStat::Data(&control_trials);
    let control_posterior = prior.posterior(&control_data);
    let llk_control = control_posterior.ln_m(&control_data);

    let exp_trials = multi_counts_to_trials(
        &modification_counts.exp_counts,
        &mods_to_index,
        modification_counts.exp_total,
    )?;
    let exp_data = DataOrSuffStat::Data(&exp_trials);
    let exp_posterior = prior.posterior(&exp_data);
    let llk_exp = exp_posterior.ln_m(&exp_data);

    let all_trials = control_trials
        .into_iter()
        .chain(exp_trials.into_iter())
        .collect::<Vec<usize>>();
    let all_data = DataOrSuffStat::Data(&all_trials);
    let all_post = prior.posterior(&all_data);
    let llk_same = all_post.ln_m(&all_data);

    // debug!("control dist:{} exp dist:{} combined dist:{}",
    //     control_posterior.to_string(),
    //     exp_posterior.to_string(),
    //     all_post.to_string());

    Ok(llk_control + llk_exp - llk_same)
}

fn counts_to_trials(count_methyl: usize, count_canonical: usize) -> Vec<bool> {
    let mut x = vec![true; count_methyl];
    let mut y = vec![false; count_canonical];
    x.append(&mut y);
    x
}

fn llk_beta(modification_counts: &ModificationCounts) -> anyhow::Result<f64> {
    let all_mods = modification_counts
        .control_counts
        .keys()
        .copied()
        .chain(modification_counts.exp_counts.keys().copied())
        .collect::<HashSet<char>>();
    if all_mods.len() != 1 {
        bail!("should have exactly one modification")
    }
    let raw_mod_code = all_mods.into_iter().take(1).collect::<Vec<char>>()[0];

    let control_methyls = *modification_counts
        .control_counts
        .get(&raw_mod_code)
        .unwrap_or(&0);
    assert!(control_methyls <= modification_counts.control_total);
    let control_canonicals =
        modification_counts.control_total - control_methyls;

    let control_trials = counts_to_trials(control_methyls, control_canonicals);
    let control_data = DataOrSuffStat::Data(&control_trials);
    let prior = Beta::jeffreys();
    let control_posterior = prior.posterior(&control_data);
    let llk_control = control_posterior.ln_m(&control_data);

    let exp_methyls = *modification_counts
        .exp_counts
        .get(&raw_mod_code)
        .unwrap_or(&0);
    assert!(exp_methyls <= modification_counts.exp_total);
    let exp_canonicals = modification_counts.exp_total - exp_methyls;
    let exp_trials = counts_to_trials(exp_methyls, exp_canonicals);
    let exp_data = DataOrSuffStat::Data(&exp_trials);
    let exp_posterior = prior.posterior(&exp_data);
    let llk_exp = exp_posterior.ln_m(&exp_data);

    let all_trials = counts_to_trials(
        exp_methyls + control_methyls,
        exp_canonicals + control_canonicals,
    );
    let all_data = DataOrSuffStat::Data(&all_trials);
    let combined_posterior = prior.posterior(&all_data);
    let llk_same = combined_posterior.ln_m(&all_data);

    Ok(llk_control + llk_exp - llk_same)
}

pub(crate) fn llk_ratio(
    modification_counts: &ModificationCounts,
) -> anyhow::Result<f64> {
    // plus one for canonical
    let n_categories = std::cmp::max(
        modification_counts.control_counts.keys().len(),
        modification_counts.exp_counts.keys().len(),
    ) + 1;
    if n_categories < 2 {
        // todo this should return an Option or something when counts is empty
        debug!("need more than one category, got {:?}", modification_counts);
        return Ok(0f64);
    }
    if n_categories == 2 {
        llk_beta(modification_counts)
    } else {
        llk_dirichlet(modification_counts)
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

pub(crate) fn chi_squared_test(
    modification_counts: &ModificationCounts,
) -> anyhow::Result<f64> {
    if modification_counts.control_total == 0
        || modification_counts.exp_total == 0
    {
        info!("site has zero total");
        return Ok(1f64);
    }
    let n_conditions = 2;
    let mods_to_index = modification_counts
        .control_counts
        .keys()
        .chain(modification_counts.exp_counts.keys())
        .copied()
        .collect::<HashSet<char>>()
        .into_iter()
        .sorted_by(|a, b| a.cmp(b))
        .enumerate()
        .map(|(i, c)| (c, i + 1))
        .collect::<HashMap<char, usize>>();
    let n_cats = mods_to_index.len() + 1;
    let obs = vec![
        get_counts_each(
            &modification_counts.control_counts,
            &mods_to_index,
            modification_counts.control_total,
        ),
        get_counts_each(
            &modification_counts.exp_counts,
            &mods_to_index,
            modification_counts.exp_total,
        ),
    ];

    let data = obs.iter().flatten().copied().collect::<Vec<u32>>();
    let contingency = Array2::from_shape_vec((n_conditions, n_cats), data)
        .map(|a| a.t().to_owned().mapv(|x| f64::from(x + 1)))?;

    // let condition_totals = contingency.sum_axis(Axis(0)).mapv(|x| f64::from(x));
    let category_totals = contingency.sum_axis(Axis(1)).mapv(|x| f64::from(x));
    let total: f64 = category_totals.sum();
    let category_proportions = &category_totals / total;
    let stat =
        contingency
            .columns()
            .into_iter()
            .fold(0f64, |acc, cond_counts| {
                let cond_total = cond_counts.sum();
                let expected = cond_total * &category_proportions;
                let stat =
                    (&cond_counts - &expected).mapv(|x| x.powi(2)) / expected;
                let stat = stat.sum();
                acc + stat
            });
    let df = (n_conditions - 1) * (n_cats - 1);
    let chi_dist = ChiSquared::new(df as f64)?;
    let p_val = chi_dist.sf(&stat);

    Ok(p_val)
}

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

pub(crate) fn mle_llk_ratio(
    modification_counts: &ModificationCounts,
) -> anyhow::Result<f64> {
    if modification_counts.control_total == 0
        || modification_counts.exp_total == 0
    {
        info!("site has zero total");
        return Ok(0f64);
    }
    let mods_to_index = modification_counts
        .control_counts
        .keys()
        .chain(modification_counts.exp_counts.keys())
        .copied()
        .collect::<HashSet<char>>()
        .into_iter()
        .sorted_by(|a, b| a.cmp(b))
        .enumerate()
        .map(|(i, c)| (c, i + 1))
        .collect::<HashMap<char, usize>>();

    let alphas_control = modified_counts_to_frequencies(
        &modification_counts.control_counts,
        &mods_to_index,
        modification_counts.control_total,
    );
    let alphas_exp = modified_counts_to_frequencies(
        &modification_counts.exp_counts,
        &mods_to_index,
        modification_counts.exp_total,
    );
    let combined_counts = modification_counts
        .control_counts
        .iter()
        .chain(modification_counts.exp_counts.iter())
        .fold(HashMap::new(), |mut acc, (mod_code, count)| {
            *acc.entry(*mod_code).or_insert(0) += *count;
            acc
        });
    let alphas_combined = modified_counts_to_frequencies(
        &combined_counts,
        &mods_to_index,
        modification_counts.exp_total + modification_counts.control_total,
    );

    let control_trials = multi_counts_to_trials(
        &modification_counts.control_counts,
        &mods_to_index,
        modification_counts.control_total,
    )?;
    let control_data = CategoricalData::Data(&control_trials);
    let lk_control = Dirichlet::new(alphas_control.clone())
        .with_context(|| {
            format!(
                "invalid control Dir params {:?}, counts: {:?} {:?}",
                &alphas_control,
                &modification_counts.control_counts,
                &modification_counts.control_total
            )
        })?
        .ln_m(&control_data);

    let exp_trials = multi_counts_to_trials(
        &modification_counts.exp_counts,
        &mods_to_index,
        modification_counts.exp_total,
    )?;
    let exp_data = CategoricalData::Data(&exp_trials);
    let lk_exp = Dirichlet::new(alphas_exp.clone())
        .with_context(|| format!("invalid exp Dir params {:?}", &alphas_exp))?
        .ln_m(&exp_data);

    let combined_trials = multi_counts_to_trials(
        &combined_counts,
        &mods_to_index,
        modification_counts.exp_total + modification_counts.control_total,
    )?;
    let combined_data = CategoricalData::Data(&combined_trials);
    let lk_combined = Dirichlet::new(alphas_combined.clone())
        .with_context(|| {
            format!("invalid combined Dir params {:?}", &alphas_combined)
        })?
        .ln_m(&combined_data);

    Ok(lk_control + lk_exp - lk_combined)
}

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
