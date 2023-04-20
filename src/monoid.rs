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

impl<A> Moniod for HashMap<A, u64>
where
    A: Eq + PartialEq + Hash,
{
    fn zero() -> Self {
        HashMap::new()
    }

    fn op(self, other: Self) -> Self {
        let mut agg = self;
        for (k, count) in other {
            *agg.entry(k).or_insert(0) += count;
        }
        agg
    }

    fn op_mut(&mut self, other: Self) {
        for (k, count) in other {
            *self.entry(k).or_insert(0) += count;
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl<A, B> Moniod for HashMap<A, B>
where
    A: Eq + Hash,
    B: Moniod,
{
    fn zero() -> Self {
        HashMap::new()
    }

    fn op(self, other: Self) -> Self {
        let mut agg = self;
        for (k, v) in other {
            agg.entry(k).or_insert(B::zero()).op_mut(v);
        }
        agg
    }

    fn op_mut(&mut self, other: Self) {
        for (k, v) in other {
            self.entry(k).or_insert(B::zero()).op_mut(v)
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}

#[cfg(test)]
mod moniod_tests {
    use crate::monoid::Moniod;
    use std::collections::HashMap;

    #[test]
    fn test_moniod_hashmap_counts() {
        let a = vec![('a', 1), ('b', 2), ('c', 10)]
            .into_iter()
            .collect::<HashMap<char, u64>>();
        let b = vec![('a', 9), ('b', 8), ('c', 10), ('d', 1)]
            .into_iter()
            .collect::<HashMap<char, u64>>();

        let expected = vec![('a', 10), ('b', 10), ('c', 20), ('d', 1)]
            .into_iter()
            .collect::<HashMap<char, u64>>();

        let c = a.op(b);
        assert_eq!(&c, &expected);

        let a = vec![('a', 1), ('b', 2), ('c', 10)]
            .into_iter()
            .collect::<HashMap<char, u64>>();
        let b = vec![('a', 9), ('b', 8), ('c', 10), ('d', 1)]
            .into_iter()
            .collect::<HashMap<char, u64>>();
        let mut c = HashMap::new();
        let mut d = HashMap::new();
        c.insert("foo".to_string(), a);
        d.insert("foo".to_string(), b.clone());
        d.insert("bar".to_string(), b.clone());
        let e = c.op(d);

        assert_eq!(e.get("foo").unwrap(), &expected);
        assert_eq!(e.get("bar").unwrap(), &b);
    }
}
