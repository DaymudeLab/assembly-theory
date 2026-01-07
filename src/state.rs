//! Assembly state of the top-down recursive search algorithm.

use bit_set::BitSet;

//use crate::molecule::Molecule;
use crate::object::AObject;
/// Assembly state of the top-down recursive search algorithm.
pub struct State {
    /// List of connected components in this state.
    fragments: Vec<BitSet>,
    /// The current upper bound on the assembly index, i.e., edges(mol) - 1 -
    /// [edges(subgraphs removed) - #(subgraphs removed)].
    index: usize,
    /// The indices of previously removed duplicate subgraphs. Used for
    /// disambiguating the serial order of two states during memoization.
    removal_order: Vec<usize>,
    /// Size of the largest duplicatable subgraph removed up to this point.
    /// Since matches are removed in decreasing size, this provides an upper
    /// bound on the largest subgraph that can be removed from this state.
    largest_removed: usize,
    /// Index of the last removed match; initially `-1`.
    last_removed: isize,
}

impl State {
    /// Construct the initial [`State`] for the given molecule.
    pub fn new<T: AObject>(mol: &T) -> Self {
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

    /// Construct the child [`State`] resulting from removing the specified
    /// match from this [`State`].
    pub fn update(
        &self,
        fragments: Vec<BitSet>,
        remove_ix: usize,
        match_ix: usize,
        remove_len: usize,
    ) -> Self {
        Self {
            fragments,
            index: self.index - remove_len + 1,
            removal_order: {
                let mut clone = self.removal_order.clone();
                clone.push(remove_ix);
                clone
            },
            largest_removed: remove_len,
            last_removed: match_ix as isize,
        }
    }

    /// Return a reference to this [`State`]'s list of connected components.
    pub fn fragments(&self) -> &Vec<BitSet> {
        &self.fragments
    }

    /// Return this [`State`]'s upper bound on the assembly index.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Return a reference to this [`State`]'s removal order (used for
    /// disambiguating state ordering during memoization).
    pub fn removal_order(&self) -> &Vec<usize> {
        &self.removal_order
    }

    /// Return the size of this [`State`]'s largest duplicatable subgraph
    /// removed so far.
    pub fn largest_removed(&self) -> usize {
        self.largest_removed
    }

    /// Return the index of this [`State`]'s last removed match.
    pub fn last_removed(&self) -> isize {
        self.last_removed
    }
}
