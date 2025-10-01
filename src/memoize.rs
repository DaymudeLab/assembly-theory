//! Memoize assembly states to avoid redundant recursive search.

use std::sync::{
    atomic::{AtomicUsize, Ordering::Relaxed},
    Arc,
};

use bit_set::BitSet;
use clap::ValueEnum;
use dashmap::DashMap;

use crate::{
    canonize::{canonize, CanonizeMode, Labeling},
    molecule::Molecule,
    state::State,
};

/// Strategy for memoizing assembly states in the search phase.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum MemoizeMode {
    /// Do not use memoization.
    None,
    /// Cache states by fragments' canonical labelings, allowing isomorphic
    /// assembly states to hash to the same value.
    CanonIndex,
}

/// Struct for the memoization cache.
#[derive(Clone)]
pub struct Cache {
    /// Memoization mode.
    memoize_mode: MemoizeMode,
    /// Canonization mode; only used with [`MemoizeMode::CanonIndex`].
    canonize_mode: CanonizeMode,
    /// A parallel-aware cache mapping keys (lists of canon ids) to their
    /// assembly index upper bounds and match removal order.
    cache: Arc<DashMap<Vec<usize>, (usize, Vec<usize>)>>,
    /// A parallel-aware map from a graph canonization labeling to a canon id.
    /// These ids are used as keys into the cache since hashing a usize is much
    /// faster than hashing the labelings directly.
    label_to_canon_id: Arc<DashMap<Labeling, usize>>,
    /// A parallel-aware map from fragments (represented as bitsets) to canon
    /// ids. Two fragments have the same id iff they are isomorphic.
    frag_to_canon_id: Arc<DashMap<BitSet, usize>>,
    /// A counter for what id should be given to the next unique labeling seen.
    next_id: Arc<AtomicUsize>,
}

impl Cache {
    /// Construct a new [`Cache`] with the specified modes.
    pub fn new(memoize_mode: MemoizeMode, canonize_mode: CanonizeMode) -> Self {
        Self {
            memoize_mode,
            canonize_mode,
            cache: Arc::new(DashMap::<Vec<usize>, (usize, Vec<usize>)>::new()),
            label_to_canon_id: Arc::new(DashMap::<Labeling, usize>::new()),
            frag_to_canon_id: Arc::new(DashMap::<BitSet, usize>::new()),
            next_id: Arc::new(AtomicUsize::from(0)),
        }
    }

    /// Create a key into the cache for the given assembly state.
    ///
    /// If using [`MemoizeMode::CanonIndex`], keys are sorted lists of canon ids.
    fn key(&mut self, mol: &Molecule, state: &State) -> Option<Vec<usize>> {
        match self.memoize_mode {
            MemoizeMode::None => None,
            MemoizeMode::CanonIndex => {
                let mut frag_ids = Vec::new();
                for frag in state.fragments() {
                    let id = self.get_canon_id(mol, frag);
                    frag_ids.push(id);
                }
                frag_ids.sort();

                Some(frag_ids)
            }
        }
    }

    /// Takes a fragment and canonizes it using the specified [`CanonizeMode`],
    /// then returns the corresponding canon id.
    fn get_canon_id(&mut self, mol: &Molecule, frag: &BitSet) -> usize {
        // If frag has id, use it
        if let Some(x) = self.frag_to_canon_id.get(frag) {
            *x
        }
        // Otherwise canonize to get labeling
        else {
            let canon = canonize(mol, frag, self.canonize_mode);
            // If label has id, use it
            if let Some(x) = self.label_to_canon_id.get(&canon) {
                let id = *x;
                self.frag_to_canon_id.insert(frag.clone(), id);
                id
            }
            // Otherwise asign new id
            else {
                let id = self.next_id.fetch_add(1, Relaxed);
                self.label_to_canon_id.insert(canon, id);
                self.frag_to_canon_id.insert(frag.clone(), id);
                id
            }
        }
    }

    /// Return `true` iff memoization is enabled and this assembly state is
    /// preempted by the cached assembly state.
    /// See https://github.com/DaymudeLab/assembly-theory/pull/95 for details.
    pub fn memoize_state(&mut self, mol: &Molecule, state: &State) -> bool {
        let state_index = state.index();
        let removal_order = state.removal_order();

        // true if this state has been memoized by the cache and thus we can
        // stop computation along this path.
        let mut prune_state = false;

        // If memoization is enabled, get this assembly state's cache key.
        if let Some(cache_key) = self.key(mol, state) {
            // Do all of the following atomically: Access the cache entry. If
            // the cached entry has a worse index upper bound or later removal
            // order than this state, or if it does not exist, then cache this
            // state's values and return `false`. Otherwise, the cached entry
            // preempts this assembly state, so return `true`.
            self.cache
                .entry(cache_key)
                .and_modify(|val| {
                    if val.0 > state_index || val.1 > *removal_order {
                        val.0 = state_index;
                        val.1 = removal_order.clone();
                    } else {
                        prune_state = true;
                    }
                })
                .or_insert((state_index, removal_order.clone()));
        }

        prune_state
    }
}
