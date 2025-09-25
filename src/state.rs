use bit_set::BitSet;

use crate::molecule::Molecule;

pub struct State {
    fragments: Vec<BitSet>,
    index: usize,
    removal_order: Vec<usize>,
    largest_removed: usize,
    last_removed: isize,
}

impl State {
    pub fn new(mol: &Molecule) -> Self {
        Self {
            fragments: {
                let mut init = BitSet::new();
                init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
                vec![init]
            },
            index: mol.graph().edge_count() - 1,
            removal_order: Vec::new(),
            largest_removed: mol.graph().edge_count(),
            last_removed: -1,
        }
    }

    pub fn update(&self, fragments: Vec<BitSet>, remove_idx: usize, remove_len: usize) -> Self {
        Self {
            fragments,
            index: self.index - remove_len + 1,
            removal_order: {
                let mut clone = self.removal_order.clone();
                clone.push(remove_idx);
                clone
            },
            largest_removed: remove_len,
            last_removed: self.last_removed + (remove_idx as isize) + 1,
        }
    }

    pub fn fragments(&self) -> &Vec<BitSet> {
        &self.fragments
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn removal_order(&self) -> &Vec<usize> {
        &self.removal_order
    }

    pub fn largest_removed(&self) -> usize {
        self.largest_removed
    }

    pub fn last_removed(&self) -> isize {
        self.last_removed
    }
}