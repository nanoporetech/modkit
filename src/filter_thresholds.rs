use crate::mod_base_code::DnaBase;
use derive_new::new;
use std::collections::HashMap;

#[derive(Debug, new)]
pub struct FilterThresholds {
    default: f32,
    per_base_thresholds: HashMap<DnaBase, f32>,
}

impl FilterThresholds {
    pub fn new_passthrough() -> Self {
        Self::new(0f32, HashMap::new())
    }

    pub fn get(&self, base: &DnaBase) -> f32 {
        *self.per_base_thresholds.get(&base).unwrap_or(&self.default)
    }

    pub fn iter_thresholds(&self) -> impl Iterator<Item = (&DnaBase, &f32)> {
        self.per_base_thresholds.iter()
    }
}
