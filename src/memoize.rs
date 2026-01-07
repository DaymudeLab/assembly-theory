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
    object::AObject,
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
    /// A parallel-aware cache mapping keys (lists of usize canonical IDs) to
    /// their assembly index upper bounds and match removal order.
    #[allow(clippy::type_complexity)]
    cache: Arc<DashMap<Vec<usize>, (usize, Vec<usize>)>>,
    /// A parallel-aware map from canonical labelings to canonical IDs. Lists
    /// of these IDs are used as memoization cache keys since usizes are much
    /// faster to hash than canonical labelings.
    labeling_to_id: Arc<DashMap<Labeling, usize>>,
    /// A parallel-aware map from fragments to canonical IDs. Two fragments
    /// have the same canonical ID iff they are isomorphic.
    fragment_to_id: Arc<DashMap<BitSet, usize>>,
    /// A parallel-aware counter for assigning a unique ID to the next unique
    /// canonical labeling seen.
    next_id: Arc<AtomicUsize>,
}

impl Cache {
    /// Construct a new [`Cache`] with the specified modes.
    pub fn new(memoize_mode: MemoizeMode, canonize_mode: CanonizeMode) -> Self {
        Self {
            memoize_mode,
            canonize_mode,
            cache: Arc::new(DashMap::<Vec<usize>, (usize, Vec<usize>)>::new()),
            labeling_to_id: Arc::new(DashMap::<Labeling, usize>::new()),
            fragment_to_id: Arc::new(DashMap::<BitSet, usize>::new()),
            next_id: Arc::new(AtomicUsize::from(0)),
        }
    }

    /// Create a memoization cache key for the given assembly state.
    ///
    /// If using [`MemoizeMode::CanonIndex`], keys are sorted lists of
    /// canonical IDs.
    fn key<T: AObject>(&mut self, mol: &T, state: &State) -> Option<Vec<usize>> {
        match self.memoize_mode {
            MemoizeMode::None => None,
            MemoizeMode::CanonIndex => {
                let mut fragment_ids: Vec<usize> = state
                    .fragments()
                    .iter()
                    .map(|fragment| self.canonical_id(mol, fragment))
                    .collect();
                fragment_ids.sort();

                Some(fragment_ids)
            }
        }
    }

    /// Obtain the canonical ID of the given fragment, canonizing it using the
    /// specified [`CanonizeMode`] if this has not already been done.
    fn canonical_id<T: AObject>(&mut self, mol: &T, fragment: &BitSet) -> usize {
        if let Some(id) = self.fragment_to_id.get(fragment) {
            *id
        } else {
            let labeling = canonize(mol, fragment, self.canonize_mode);
            if let Some(id) = self.labeling_to_id.get(&labeling) {
                self.fragment_to_id.insert(fragment.clone(), *id);
                *id
            } else {
                let id = self.next_id.fetch_add(1, Relaxed);
                self.fragment_to_id.insert(fragment.clone(), id);
                self.labeling_to_id.insert(labeling, id);
                id
            }
        }
    }

    /// Return `true` iff memoization is enabled and this assembly state is
    /// preempted by a cached assembly state. See
    /// <https://github.com/DaymudeLab/assembly-theory/pull/95> for details.
    pub fn memoize_state<T: AObject>(&mut self, mol: &T, state: &State) -> bool {
        let state_index = state.index();
        let removal_order = state.removal_order();
        let mut result = false;

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
                        result = true;
                    }
                })
                .or_insert((state_index, removal_order.clone()));
        }

        result
    }
}
