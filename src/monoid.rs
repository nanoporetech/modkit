use std::collections::HashMap;
use std::hash::Hash;

pub trait Moniod {
    fn zero() -> Self;
    fn op(self, other: Self) -> Self;
    fn op_mut(&mut self, other: Self);
    fn len(&self) -> usize;
}

impl<A, B> Moniod for HashMap<A, Vec<B>>
where
    A: Eq + Hash,
{
    fn zero() -> Self {
        HashMap::new()
    }

    fn op(self, other: Self) -> Self {
        let mut out = Self::zero();
        for (k, mut vs) in self.into_iter().chain(other.into_iter()) {
            let agg = out.entry(k).or_insert(Vec::new());
            agg.append(&mut vs);
        }
        out
    }

    fn op_mut(&mut self, other: Self) {
        for (k, mut vs) in other {
            let agg = self.entry(k).or_insert(Vec::new());
            agg.append(&mut vs);
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}
