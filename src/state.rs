use bit_set::BitSet;

use crate::molecule::Molecule;

/// Struct containing information local to a state in the search space.
pub struct State {
    /// List of connected components in this state.
    fragments: Vec<BitSet>,
    /// This assembly state's upper bound on the assembly index,
    /// i.e., edges(mol) - 1 - [edges(subgraphs removed) - #(subgraphs removed)].
    index: usize,
    /// List of indices, one for each removed duplicate subgraph.
    /// Allows two states to be compared to determine which comes first in the serial removal order.
    /// Used to prevent memoization bugs.
    removal_order: Vec<usize>,
    /// Size of largest removed duplicatable subgraph up to this point.
    /// Since matches are removed in decreasing size, provides an upper bound
    /// on the largest subgraph to be removed from this state.
    largest_removed: usize,
    /// Index of the last removed duplicate subgraph match.
    /// Initialized to -1 before any match removals.
    last_removed: isize,
}

impl State {
    /// Initializes new state struct
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

    /// Creates a struct for a new child state after a duplicate subgraph match removal.
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