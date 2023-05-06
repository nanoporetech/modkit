use crate::mod_bam::{BaseModCall, BaseModProbs};
use crate::mod_base_code::{DnaBase, ModCode};
use derive_new::new;
use log::debug;
use std::collections::HashMap;

#[derive(new)]
pub struct MultipleThresholdModCaller {
    per_base_thresholds: HashMap<DnaBase, f32>,
    per_mod_thresholds: HashMap<ModCode, f32>,
    default_threshold: f32,
}

impl MultipleThresholdModCaller {
    pub fn new_passthrough() -> Self {
        Self {
            per_base_thresholds: HashMap::new(),
            per_mod_thresholds: HashMap::new(),
            default_threshold: 0f32,
        }
    }

    /// Make a base modification call from the probabilities of each modification class.
    /// Result will be Err if the raw mod code cannot be parsed (this will change in the future,
    /// when BaseModProbs don't need to be converted to ModCodes.
    pub fn call(
        &self,
        canonical_base: &DnaBase,
        base_mod_probs: &BaseModProbs,
    ) -> anyhow::Result<BaseModCall> {
        let mod_code_to_probs = base_mod_probs
            .iter_probs()
            .map(|(raw_mod_code, p)| {
                ModCode::parse_raw_mod_code(*raw_mod_code).map(|mc| (mc, *p))
            })
            .collect::<Result<Vec<(ModCode, f32)>, _>>()?;

        let mut filtered_probs = mod_code_to_probs
            .into_iter()
            .filter_map(|(mod_code, p_mod)| {
                let threshold = self
                    .per_mod_thresholds
                    .get(&mod_code)
                    .or(self.per_base_thresholds.get(canonical_base))
                    .unwrap_or(&self.default_threshold);
                if p_mod >= *threshold {
                    Some(BaseModCall::Modified(p_mod, mod_code))
                } else {
                    None
                }
            })
            .collect::<Vec<BaseModCall>>();

        let canonical_threshold = self
            .per_mod_thresholds
            .get(&canonical_base.canonical_mod_code()?)
            .or(self.per_base_thresholds.get(canonical_base))
            .unwrap_or(&self.default_threshold);
        if base_mod_probs.canonical_prob() >= *canonical_threshold {
            filtered_probs
                .push(BaseModCall::Canonical(base_mod_probs.canonical_prob()))
        };

        Ok(filtered_probs
            .into_iter()
            .max()
            .unwrap_or(BaseModCall::Filtered))
    }

    /// Use thresholds to convert base modification probabilities into a "call", where
    /// the probabilities are 1.0 for the predicted class. None is returned when the
    /// probabilities all fail to meet the threshold requirements
    pub fn call_probs(
        &self,
        canonical_base: &DnaBase,
        mut base_mod_probs: BaseModProbs,
    ) -> anyhow::Result<Option<BaseModProbs>> {
        let base_mod_call = self.call(canonical_base, &base_mod_probs)?;
        match base_mod_call {
            BaseModCall::Modified(_, mod_code) => {
                Ok(Some(BaseModProbs::new(mod_code.char(), 1.0)))
            }
            BaseModCall::Canonical(_) => {
                base_mod_probs.iter_mut_probs().for_each(|p| *p = 0f32);
                Ok(Some(base_mod_probs))
            }
            BaseModCall::Filtered => Ok(None),
        }
    }
}

#[cfg(test)]
mod threshold_mod_caller_tests {
    use crate::mod_bam::{BaseModCall, BaseModProbs};
    use crate::mod_base_code::{DnaBase, ModCode};
    use crate::threshold_mod_caller::MultipleThresholdModCaller;
    use anyhow::anyhow;
    use indexmap::indexset;
    use std::collections::HashMap;

    fn assert_base_mod_call_canonical(
        base_mod_call: BaseModCall,
        p: f32,
    ) -> anyhow::Result<()> {
        match base_mod_call {
            BaseModCall::Canonical(q) => {
                if (p - q).abs() < f32::EPSILON {
                    Ok(())
                } else {
                    Err(anyhow!("expected {p} got {q}"))
                }
            }
            _ => Err(anyhow!("expected canonical got {:?}", base_mod_call)),
        }
    }

    fn assert_base_mod_call_modified(
        base_mod_call: BaseModCall,
        p: f32,
        mod_code: ModCode,
    ) -> anyhow::Result<()> {
        match base_mod_call {
            BaseModCall::Modified(q, code) => {
                if code != mod_code {
                    Err(anyhow!("expected {:?} got {:?}", mod_code, code))
                } else if (p - q).abs() < f32::EPSILON {
                    Ok(())
                } else {
                    Err(anyhow!("expected {p} got {q}"))
                }
            }
            _ => Err(anyhow!("expected modified got {:?}", base_mod_call)),
        }
    }

    #[test]
    fn test_multi_threshold_call_semantics() {
        // CASE A
        let per_mod_threshold = vec![(ModCode::A, 0.8), (ModCode::a, 0.9)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.8);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        // neither pass, Filtered call
        assert_eq!(call, BaseModCall::Filtered);
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Canonical(0.8));
        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.9, ModCode::a));

        // CASE B
        let per_mod_threshold = vec![(ModCode::A, 0.2), (ModCode::a, 0.9)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        // have to make this 0.79 because of some FP wobble
        let base_modprobs = BaseModProbs::new('a', 0.79);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        // call canonical, 'a' fails, but p_A is >= threshold
        assert_base_mod_call_canonical(call, 0.21).unwrap();

        let base_modprobs = BaseModProbs::new('a', 0.6);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        // same logic as above
        assert_base_mod_call_canonical(call, 0.4).unwrap();

        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_base_mod_call_canonical(call, 0.8).unwrap();

        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_base_mod_call_modified(call, 0.9, ModCode::a).unwrap();

        // CASE C
        let per_mod_threshold = vec![(ModCode::A, 0.2), (ModCode::a, 0.8)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.8);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.8, ModCode::a));
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Canonical(0.8));
        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.9, ModCode::a));
    }

    #[test]
    fn test_multi_threshold_passthrough() {
        let caller = MultipleThresholdModCaller::new_passthrough();
        let base_modprobs = BaseModProbs::new('a', 0.8);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.8, ModCode::a));
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Canonical(0.8));
    }

    #[test]
    fn test_multi_threshold_base_threshold() {
        let per_mod_threshold = vec![(ModCode::a, 0.8)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let per_base_threshold = vec![(DnaBase::A, 0.7)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            per_base_threshold,
            per_mod_threshold,
            0.75f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.75);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Filtered);
        let base_modprobs = BaseModProbs::new('a', 0.6);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Filtered);
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller.call(&DnaBase::A, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Canonical(0.8));
        let base_modprobs = BaseModProbs::new('m', 0.8);
        let call = caller.call(&DnaBase::C, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.8, ModCode::m));
        let base_modprobs = BaseModProbs::new('m', 0.72);
        let call = caller.call(&DnaBase::C, &base_modprobs).unwrap();
        assert_eq!(call, BaseModCall::Filtered);
    }

    #[test]
    fn test_multi_threshold_call_probs() {
        // CASE A
        let per_mod_threshold = vec![(ModCode::A, 0.8), (ModCode::a, 0.9)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.8);
        let call = caller.call_probs(&DnaBase::A, base_modprobs).unwrap();
        // neither pass, Filtered call
        assert!(call.is_none());

        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        // canonical call
        assert_eq!(call, BaseModProbs::new('a', 0f32));
        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        // modified call
        assert_eq!(call, BaseModProbs::new('a', 1.0));

        // CASE B
        let per_mod_threshold = vec![(ModCode::A, 0.2), (ModCode::a, 0.9)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.79);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        // canonical call
        assert_eq!(call, BaseModProbs::new('a', 0f32));
        let base_modprobs = BaseModProbs::new('a', 0.6);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        // same, canonical call
        assert_eq!(call, BaseModProbs::new('a', 0f32));
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('a', 0f32));
        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('a', 1.0));

        // CASE C
        let per_mod_threshold = vec![(ModCode::A, 0.2), (ModCode::a, 0.8)]
            .into_iter()
            .collect::<HashMap<_, _>>();
        let caller = MultipleThresholdModCaller::new(
            HashMap::new(),
            per_mod_threshold,
            0f32,
        );
        let base_modprobs = BaseModProbs::new('a', 0.8);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('a', 1.0));
        let base_modprobs = BaseModProbs::new('a', 0.2);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('a', 0f32));
        let base_modprobs = BaseModProbs::new('a', 0.9);
        let call = caller
            .call_probs(&DnaBase::A, base_modprobs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('a', 1.0));
    }

    #[test]
    fn test_multi_threshold_call_multiple_mods_semantics() {
        let per_mod_thresholds = vec![(ModCode::m, 0.7), (ModCode::h, 0.8)]
            .into_iter()
            .collect();
        let per_base_thresholds =
            vec![(DnaBase::C, 0.75)].into_iter().collect();
        let caller = MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds,
            0f32,
        );
        let mut base_mod_probs = BaseModProbs::new('m', 0.1);
        base_mod_probs.insert_base_mod_prob('h', 0.8);
        let call = caller.call(&DnaBase::C, &base_mod_probs).unwrap();
        assert_eq!(call, BaseModCall::Modified(0.8, ModCode::h));

        let mut base_mod_probs = BaseModProbs::new('m', 0.2);
        base_mod_probs.insert_base_mod_prob('h', 0.7);
        let call = caller.call(&DnaBase::C, &base_mod_probs).unwrap();
        assert_eq!(call, BaseModCall::Filtered);

        let per_mod_thresholds = vec![(ModCode::m, 0.7), (ModCode::h, 0.8)]
            .into_iter()
            .collect();
        let per_base_thresholds = vec![(DnaBase::C, 0.1)].into_iter().collect();
        let caller = MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds,
            0f32,
        );
        let mut base_mod_probs = BaseModProbs::new('m', 0.2);
        base_mod_probs.insert_base_mod_prob('h', 0.7);
        let call = caller.call(&DnaBase::C, &base_mod_probs).unwrap();
        assert_base_mod_call_canonical(call, 0.1).unwrap();
    }

    #[test]
    fn test_multi_threshold_call_probs_multiple_mods_semantics() {
        let per_mod_thresholds = vec![(ModCode::m, 0.7), (ModCode::h, 0.8)]
            .into_iter()
            .collect();
        let per_base_thresholds =
            vec![(DnaBase::C, 0.75)].into_iter().collect();
        let caller = MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds,
            0f32,
        );
        let mut base_mod_probs = BaseModProbs::new('m', 0.1);
        base_mod_probs.insert_base_mod_prob('h', 0.8);
        let call = caller
            .call_probs(&DnaBase::C, base_mod_probs)
            .unwrap()
            .unwrap();
        assert_eq!(call, BaseModProbs::new('h', 1.0));

        let mut base_mod_probs = BaseModProbs::new('m', 0.2);
        base_mod_probs.insert_base_mod_prob('h', 0.7);
        let call = caller.call_probs(&DnaBase::C, base_mod_probs).unwrap();
        assert!(call.is_none());

        let per_mod_thresholds = vec![(ModCode::m, 0.7), (ModCode::h, 0.8)]
            .into_iter()
            .collect();
        let per_base_thresholds = vec![(DnaBase::C, 0.1)].into_iter().collect();
        let caller = MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds,
            0f32,
        );
        let mut base_mod_probs = BaseModProbs::new('m', 0.2);
        base_mod_probs.insert_base_mod_prob('h', 0.7);
        let call = caller
            .call_probs(&DnaBase::C, base_mod_probs)
            .unwrap()
            .unwrap();
        let mut expected_base_mod_probs = BaseModProbs::new('m', 0f32);
        expected_base_mod_probs.insert_base_mod_prob('h', 0f32);
        assert_eq!(call, expected_base_mod_probs);
    }
}
