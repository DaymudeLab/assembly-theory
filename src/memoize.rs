use bit_set::BitSet;
use clap::ValueEnum;
use dashmap::DashMap;
use std::sync::Arc;
use crate::{
    canonize::{CanonizeMode, Labeling, canonize},
    molecule::Molecule,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum CacheMode {
    None,
    Index,
    IndexCanon,
}

#[derive(Eq, Hash, PartialEq)]
enum CacheType {
    Set(Vec<BitSet>),
    Canon(Vec<Labeling>),
}

#[derive(Clone)]
pub struct Cache {
    mode: CacheMode,
    cache: Arc<DashMap<CacheType, (usize, Vec<usize>)>>,
    frags_to_labels: Arc<DashMap<BitSet, Labeling>>,
}

impl Cache {
    pub fn new(mode: CacheMode) -> Self {
        Self {
            mode,
            cache: Arc::new(DashMap::<CacheType, (usize, Vec<usize>)>::new()),
            frags_to_labels: Arc::new(DashMap::<BitSet, Labeling>::new()),
        }
    }

    pub fn get(&self, mol: &Molecule, fragments: &[BitSet], state_index: usize, order: &Vec<usize>) -> Option<usize> {
        match self.mode {
            CacheMode::None => None,
            CacheMode::Index => self.index_get(fragments, state_index, order),
            CacheMode::IndexCanon => self.index_canon_get(mol, fragments, state_index, order),
        }
    }

    fn index_get(&self, fragments: &[BitSet], state_index: usize, order: &Vec<usize>) -> Option<usize> {
        let mut frag_vec = fragments.to_vec();
        frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));

        if let Some(res) = self.cache.get(&CacheType::Set(frag_vec)) {
            let (other_index, other_order) = (*res).clone();
            if other_index <= state_index && other_order <= *order {
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

    fn index_canon_get(&self, mol: &Molecule, fragments: &[BitSet], state_index: usize, order: &Vec<usize>) -> Option<usize>{
        let mut labels: Vec<Labeling> = fragments
            .iter()
            .map(|f| {
                if let Some(res) = self.frags_to_labels.get(f) {
                    (*res).clone()
                }
                else {
                    let label = canonize(mol, f, CanonizeMode::Nauty);
                    self.frags_to_labels.insert(f.clone(), label.clone());
                    label
                }
            })
            .collect();
        labels.sort_by(|a, b| a.cmp(b));

        if let Some(res) = self.cache.get(&CacheType::Canon(labels)) {
            let (other_index, other_order) = (*res).clone();
            if other_index <= state_index && other_order <= *order {
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

    pub fn insert(&mut self, fragments: &[BitSet], state_index: usize, order: &Vec<usize>) {
        match self.mode {
            CacheMode::None => (),
            CacheMode::Index => {
                let mut frag_vec = fragments.to_vec();
                frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));
                self.cache.insert(CacheType::Set(frag_vec), (state_index, order.clone()));
            },
            CacheMode::IndexCanon => {
                let mut labels: Vec<Labeling> = fragments
                    .iter()
                    .filter_map(|f| self.frags_to_labels.get(f))
                    .map(|res| (*res).clone())
                    .collect();
                labels.sort_by(|a, b| a.cmp(b));

                self.cache.insert(CacheType::Canon(labels), (state_index, order.clone()));
            },
        };
    }
}
