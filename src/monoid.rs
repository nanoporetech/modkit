use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;

pub trait Moniod {
    fn zero() -> Self;
    fn op(self, other: Self) -> Self;
    fn op_mut(&mut self, other: Self);
    fn len(&self) -> usize;
}

pub trait BorrowingMoniod {
    fn zero() -> Self;
    fn op(self, other: &Self) -> Self;
    fn op_mut(&mut self, other: &Self);
    fn len(&self) -> usize;
}

impl<A, B> Moniod for BTreeMap<A, Vec<B>>
where
    A: Eq + Hash + Ord + PartialOrd,
{
    fn zero() -> Self {
        BTreeMap::new()
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

impl<A, B> Moniod for HashMap<A, HashSet<B>>
where
    A: Eq + Hash,
    B: Eq + Hash,
{
    fn zero() -> Self {
        HashMap::new()
    }

    fn op(self, other: Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        other.into_iter().for_each(|(a, bs)| {
            self.entry(a).or_insert(HashSet::new()).extend(bs)
        })
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

impl<A> Moniod for FxHashMap<A, u32>
where
    A: Eq + Hash,
{
    fn zero() -> Self {
        FxHashMap::default()
    }

    fn op(self, other: Self) -> Self {
        let mut agg = self;
        for (k, v) in other {
            *agg.entry(k).or_insert(0) += v;
        }
        agg
    }

    fn op_mut(&mut self, other: Self) {
        for (k, v) in other {
            *self.entry(k).or_insert(0u32) += v;
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl<A> Moniod for FxHashMap<A, usize>
where
    A: Eq + Hash,
{
    fn zero() -> Self {
        FxHashMap::default()
    }

    fn op(self, other: Self) -> Self {
        let mut agg = self;
        for (k, v) in other {
            *agg.entry(k).or_insert(0usize) += v;
        }
        agg
    }

    fn op_mut(&mut self, other: Self) {
        for (k, v) in other {
            *self.entry(k).or_insert(0usize) += v;
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl<A, B> Moniod for FxHashMap<A, B>
where
    A: Eq + Hash,
    B: Moniod,
{
    fn zero() -> Self {
        FxHashMap::default()
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

impl<A> Moniod for FxHashMap<A, u64>
where
    A: Eq + PartialEq + Hash,
{
    fn zero() -> Self {
        FxHashMap::default()
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

impl<A> Moniod for Vec<A> {
    fn zero() -> Self {
        Vec::new()
    }

    fn op(self, mut other: Self) -> Self {
        let mut this = self;
        this.append(&mut other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        self.extend(other.into_iter());
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl<A> Moniod for FxHashSet<A>
where
    A: Eq + PartialEq + Hash,
{
    fn zero() -> Self {
        FxHashSet::default()
    }

    fn op(self, other: Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        other.into_iter().for_each(|a| {
            self.insert(a);
        });
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl<A, B> BorrowingMoniod for FxHashMap<A, HashSet<B>>
where
    A: Copy + Eq + Hash,
    B: Copy + Eq + PartialEq + Hash,
{
    fn zero() -> Self {
        FxHashMap::default()
    }

    fn op(self, other: &Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: &Self) {
        other.iter().for_each(|(a, bs)| {
            self.entry(*a).or_insert(HashSet::new()).extend(bs);
        })
    }

    fn len(&self) -> usize {
        todo!()
    }
}

#[cfg(test)]
mod moniod_tests {
    use crate::monoid::Moniod;
    use rustc_hash::FxHashMap;
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

    #[test]
    fn test_monoid_op_counts() {
        let xs: Vec<Vec<u8>> =
            vec![vec![1, 1, 1, 2], vec![1, 0, 2, 2], vec![0, 2, 2, 2]];
        let counts = xs
            .iter()
            .map(|x| {
                x.iter().enumerate().fold(
                    FxHashMap::<u8, FxHashMap<usize, u32>>::default(),
                    |mut counter, (pos, elem)| {
                        *counter
                            .entry(*elem)
                            .or_insert(FxHashMap::default())
                            .entry(pos)
                            .or_insert(0u32) += 1u32;
                        counter
                    },
                )
            })
            .reduce(|a, b| a.op(b))
            .unwrap();
        assert_eq!(*counts.get(&0).unwrap().get(&0usize).unwrap(), 1);
        assert_eq!(*counts.get(&0).unwrap().get(&1usize).unwrap(), 1);
        assert_eq!(*counts.get(&1).unwrap().get(&0usize).unwrap(), 2);
        assert!(counts.get(&2).unwrap().get(&0usize).is_none());
    }
}
