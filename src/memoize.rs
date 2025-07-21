use bit_set::BitSet;
use clap::ValueEnum;
use dashmap::DashMap;

use crate::{
    canonize::{CanonizeMode, Labeling, canonize},
    molecule::Molecule,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum CacheMode {
    None,
    Index,
    Savings,
    SavingsCanon,
}

#[derive(Eq, Hash, PartialEq)]
enum CacheType {
    Set(Vec<BitSet>),
    Canon(Labeling),
}

pub struct Cache {
    mode: CacheMode,
    cache: DashMap<CacheType, usize>,
    mol: Option<Molecule>,
}

impl Cache {
    pub fn new(mode: CacheMode, mol: &Molecule) -> Self {
        let mol = {
            match mode {
                CacheMode::SavingsCanon => Some(mol.clone()),
                _ => None,
            }
        };

        Self {
            mode,
            cache: DashMap::<CacheType, usize>::new(),
            mol,
        }
    }

    pub fn get(&self, fragments: &[BitSet], state_index: usize) -> Option<usize> {
        match self.mode {
            CacheMode::None => None,
            CacheMode::Index => self.index_get(fragments, state_index),
            CacheMode::Savings => self.savings_get(fragments, state_index),
            CacheMode::SavingsCanon => self.savings_canon_get(fragments, state_index),
        }
    }

    fn index_get(&self, fragments: &[BitSet], state_index: usize) -> Option<usize> {
        let mut frag_vec = fragments.to_vec();
        frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));

        if let Some(res) = self.cache.get(&CacheType::Set(frag_vec)) {
            if *res <= state_index {
                Some(state_index)
            }
            else {
                None
            }
        }
        else {
            None
        }
    }

    fn savings_get(&self, fragments: &[BitSet], state_index: usize) -> Option<usize> {
        let mut frag_vec = fragments.to_vec();
        frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));

        if let Some(res) = self.cache.get(&CacheType::Set(frag_vec)) {
            Some(state_index - *res)
        }
        else {
            None
        }
    }

    fn savings_canon_get(&self, fragments: &[BitSet], state_index: usize) -> Option<usize> {
        let mol = self.mol.as_ref().unwrap();
        let mut subgraph = BitSet::with_capacity(mol.graph().edge_count());
        for frag in fragments {
            subgraph.union_with(frag);
        }
        let canon = canonize(mol, &subgraph, CanonizeMode::Nauty);

        if let Some(res) = self.cache.get(&CacheType::Canon(canon)) {
            Some(state_index - *res)
        }
        else {
            None
        }
    }

    pub fn insert(&mut self, fragments: &[BitSet], state_index: usize, savings: usize) {
        match self.mode {
            CacheMode::None => (),
            CacheMode::Index => {
                let mut frag_vec = fragments.to_vec();
                frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));
                self.cache.insert(CacheType::Set(frag_vec), state_index);
            },
            CacheMode::Savings => {
                let mut frag_vec = fragments.to_vec();
                frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));
                self.cache.insert(CacheType::Set(frag_vec), savings);
            },
            CacheMode::SavingsCanon => {
                let mol = self.mol.as_ref().unwrap();
                let mut subgraph = BitSet::with_capacity(mol.graph().edge_count());
                for frag in fragments {
                    subgraph.union_with(frag);
                }
                let canon = canonize(mol, &subgraph, CanonizeMode::Nauty);

                self.cache.insert(CacheType::Canon(canon), savings);
            },
        };
    }
}