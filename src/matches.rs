//! Pairs of non-overlapping, isomorphic subgraphs in a molecule.

use bit_set::BitSet;
use dashmap::DashMap;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    assembly::ParallelMode,
    canonize::{canonize, CanonizeMode, Labeling},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    molecule::Molecule,
    state::State,
};

struct DagNode {
    fragment: BitSet,
    canonical_id: usize,
    children: Vec<usize>,
}

/// Pairs of non-overlapping, isomorphic subgraphs in a molecule, sorted to
/// guarantee deterministic iteration.
pub struct Matches {
    matches: Vec<(BitSet, BitSet)>,
}

impl DagNode {
    pub fn new(fragment: BitSet, canonical_id: usize) -> Self {
        Self {
            fragment,
            canonical_id,
            children: Vec::new(),
        }
    }

    pub fn fragment(&self) -> &BitSet {
        &self.fragment
    }

    pub fn len(&self) -> usize {
        self.fragment.len()
    }
}

impl Matches {
    /// Generate [`Matches`] from the given molecule with the specified modes.
    pub fn new(
        mol: &Molecule,
        enumerate_mode: EnumerateMode,
        canonize_mode: CanonizeMode,
        parallel_mode: ParallelMode,
    ) -> Self {
        // Enumerate all connected, non-induced subgraphs with at most |E|/2
        // edges and bin them into isomorphism classes using canonization.
        let isomorphism_classes = DashMap::<Labeling, Vec<BitSet>>::new();
        let bin_subgraph = |subgraph: &BitSet| {
            isomorphism_classes
                .entry(canonize(mol, subgraph, canonize_mode))
                .and_modify(|bucket| bucket.push(subgraph.clone()))
                .or_insert(vec![subgraph.clone()]);
        };
        if parallel_mode == ParallelMode::None {
            enumerate_subgraphs(mol, enumerate_mode)
                .iter()
                .for_each(bin_subgraph);
        } else {
            enumerate_subgraphs(mol, enumerate_mode)
                .par_iter()
                .for_each(bin_subgraph);
        }

        // In each isomorphism class, get non-overlapping pairs of subgraphs.
        let mut matches = Vec::new();
        for bucket in isomorphism_classes.iter() {
            for (i, first) in bucket.iter().enumerate() {
                for second in &bucket[i..] {
                    if first.is_disjoint(second) {
                        if first > second {
                            matches.push((first.clone(), second.clone()));
                        } else {
                            matches.push((second.clone(), first.clone()));
                        }
                    }
                }
            }
        }

        // Sort pairs in a deterministic order and return.
        matches.sort_by(|e1, e2| {
            let ord = [
                e2.0.len().cmp(&e1.0.len()), // Decreasing subgraph size.
                e1.0.cmp(&e2.0),             // First subgraph lexicographical.
                e1.1.cmp(&e2.1),             // Second subgraph lexicographical.
            ];
            let mut i = 0;
            while ord[i] == std::cmp::Ordering::Equal {
                i += 1;
            }
            ord[i]
        });

        Self { matches }
    }

    /// Return the number of matches.
    pub fn len(&self) -> usize {
        self.matches.len()
    }

    /// Return `true` if there are no matches.
    pub fn is_empty(&self) -> bool {
        self.matches.is_empty()
    }

    /// Return all matches that are later than the last-removed match in the
    /// given assembly state.
    pub fn later_matches(&self, state: &State) -> &[(BitSet, BitSet)] {
        let idx = (state.last_removed() + 1) as usize;
        &self.matches[idx..]
    }
}
