use bit_set::BitSet;

use crate::molecule::Molecule;

pub struct State {
    fragments: Vec<BitSet>,
    index: usize,
    removal_order: Vec<usize>,
    last_removed: usize,
}

impl State {
    pub fn new(mol: &Molecule) -> Self {
        Self {
            fragments: {
                let mut init = BitSet::new();
                init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
                vec![init]
            },
            index: mol.graph().edge_count(),
            removal_order: Vec::new(),
            last_removed: 0,
        }
    }

    pub fn update(&self, fragments: Vec<BitSet>, removed_len: usize, removed_idx: usize, last_removed: usize) -> Self {
        Self {
            fragments,
            index: self.index - removed_len + 1,
            removal_order: {
                let mut clone = self.removal_order.clone();
                clone.push(removed_idx);
                clone
            },
            last_removed,
        }
    }

    pub fn fragments(&self) -> &[BitSet] {
        &self.fragments
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn removal_order(&self) -> &Vec<usize> {
        &self.removal_order
    }

    pub fn last_removed(&self) -> usize {
        self.last_removed
    }
}