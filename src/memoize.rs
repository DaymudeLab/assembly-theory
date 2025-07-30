//! Memoize assembly states to avoid redundant recursive search.

use std::sync::Arc;

use bit_set::BitSet;
use clap::ValueEnum;
use dashmap::DashMap;

use crate::{
    canonize::{canonize, CanonizeMode, Labeling},
    molecule::Molecule,
};

/// Strategy for memoizing assembly states in the search phase.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum MemoizeMode {
    /// Do not use memoization.
    None,
    /// Cache states by fragments and store their assembly index upper bounds.
    FragsIndex,
    /// Like `FragsIndex`, but cache states by fragments' canonical labelings,
    /// allowing isomorphic assembly states to hash to the same value.
    CanonIndex,
}

/// Key type for the memoization cache.
#[derive(Clone, PartialEq, Eq, Hash)]
enum CacheKey {
    /// Use fragments as keys, as in [`MemoizeMode::FragsIndex`].
    Frags(Vec<BitSet>),
    /// Use fragments' canonical labelings as keys, as in
    /// [`MemoizeMode::CanonIndex`].
    Canon(Vec<Labeling>),
}

/// Struct for the memoization cache.
#[derive(Clone)]
pub struct Cache {
    /// Memoization mode.
    memoize_mode: MemoizeMode,
    /// Canonization mode; only used with [`MemoizeMode::CanonIndex`].
    canonize_mode: CanonizeMode,
    /// A parallel-aware cache mapping keys (either fragments or canonical
    /// labelings, depending on the memoization mode) to their assembly index
    /// upper bounds and match removal order.
    cache: Arc<DashMap<CacheKey, (usize, Vec<usize>)>>,
    /// A parallel-aware map from fragments to their canonical labelings; only
    /// used with [`MemoizeMode::CanonIndex`].
    fragment_labels: Arc<DashMap<BitSet, Labeling>>,
}

impl Cache {
    /// Construct a new [`Cache`] with the specified modes.
    pub fn new(memoize_mode: MemoizeMode, canonize_mode: CanonizeMode) -> Self {
        Self {
            memoize_mode,
            canonize_mode,
            cache: Arc::new(DashMap::<CacheKey, (usize, Vec<usize>)>::new()),
            fragment_labels: Arc::new(DashMap::<BitSet, Labeling>::new()),
        }
    }

    /// Create a [`CacheKey`] for the given assembly state.
    ///
    /// If using [`MemoizeMode::FragsIndex`], keys are the lexicographically
    /// sorted fragment [`BitSet`]s. If using [`MemoizeMode::CanonIndex`], keys
    /// are lexicographically sorted fragment canonical labelings created using
    /// the specified [`CanonizeMode`]. These labelings are stored for reuse.
    fn key(&self, mol: &Molecule, state: &[BitSet]) -> Option<CacheKey> {
        match self.memoize_mode {
            MemoizeMode::None => None,
            MemoizeMode::FragsIndex => {
                let mut fragments = state.to_vec();
                fragments.sort_by_key(|a| a.iter().next());
                Some(CacheKey::Frags(fragments))
            }
            MemoizeMode::CanonIndex => {
                let mut labelings: Vec<Labeling> = state
                    .iter()
                    .map(|fragment| {
                        self.fragment_labels
                            .entry(fragment.clone())
                            .or_insert(canonize(mol, fragment, self.canonize_mode))
                            .value()
                            .clone()
                    })
                    .collect();
                labelings.sort();
                Some(CacheKey::Canon(labelings))
            }
        }
    }

    /// Return `true` iff memoization is enabled and this assembly state is
    /// preempted by the cached assembly state.
    /// See https://github.com/DaymudeLab/assembly-theory/pull/95 for details.
    pub fn memoize_state(
        &self,
        mol: &Molecule,
        state: &[BitSet],
        state_index: usize,
        removal_order: &Vec<usize>,
    ) -> bool {
        // If memoization is enabled, get this assembly state's cache key.
        if let Some(cache_key) = self.key(mol, state) {
            // Do all of the following atomically: Access the cache entry. If
            // the cached entry has a worse index upper bound or later removal
            // order than this state, or if it does not exist, then cache this
            // state's values and return `false`. Otherwise, the cached entry
            // preempts this assembly state, so return `true`.
            let (cached_index, cached_order) = self
                .cache
                .entry(cache_key)
                .and_modify(|val| {
                    if val.0 > state_index || val.1 > *removal_order {
                        val.0 = state_index;
                        val.1 = removal_order.clone();
                    }
                })
                .or_insert((state_index, removal_order.clone()))
                .value()
                .clone();
            if cached_index <= state_index && cached_order < *removal_order {
                return true;
            }
        }
        false
    }
}
