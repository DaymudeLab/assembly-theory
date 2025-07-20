use std::collections::HashMap;
use bit_set::BitSet;
use clap::ValueEnum;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum CacheMode {
    None,
    Index,
    Savings,
    IndexCanon,
    SavingsCanon,
}

pub struct Cache {
    mode: CacheMode,
    cache: HashMap<Vec<BitSet>, usize>,
}

impl Cache {
    pub fn new(mode: CacheMode) -> Self {
        let cache: HashMap<Vec<BitSet>, usize> = HashMap::new();
        Self {
            mode,
            cache,
        }
    }

    pub fn get(&self, fragments: &[BitSet], state_index: usize) -> Option<usize> {
        let mut frag_vec = fragments.to_vec();
        frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));

        match self.mode {
            CacheMode::None => None,
            CacheMode::Index => self.get_index_mode(frag_vec, state_index),
            CacheMode::Savings => self.get_savings_mode(frag_vec, state_index),
            _ => None,
        }
    }

    fn get_index_mode(&self, fragments: Vec<BitSet>, state_index: usize) -> Option<usize> {
        if let Some(res) = self.cache.get(&fragments) {
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

    fn get_savings_mode(&self, fragments: Vec<BitSet>, state_index: usize) -> Option<usize> {
        if let Some(res) = self.cache.get(&fragments) {
            Some(state_index - *res)
        }
        else {
            None
        }
    }

    pub fn insert(&mut self, fragments: &[BitSet], state_index: usize, savings: usize) {
        let mut frag_vec = fragments.to_vec();
        frag_vec.sort_by(|a, b| a.iter().next().cmp(&b.iter().next()));

        match self.mode {
            CacheMode::Index => {self.cache.insert(frag_vec, state_index);},
            CacheMode::Savings => {self.cache.insert(frag_vec, savings);},
            _ => (),
        };
    }
}