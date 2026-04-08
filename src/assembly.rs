//! Compute assembly indices of molecules.
//!
//! # Example
//! ```
//! # use std::{fs, path::PathBuf};
//! use assembly_theory::assembly::index;
//! use assembly_theory::loader::parse_molfile_str;
//!
//! # fn main() -> Result<(), std::io::Error> {
//! // Load a molecule from a .mol file.
//! let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
//! let molfile = fs::read_to_string(path)?;
//! let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
//!
//! // Compute the molecule's assembly index.
//! assert_eq!(index(&anthracene), 6);
//! # Ok(())
//! # }
//! ```

use std::{
    collections::HashSet,
    sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc, Mutex,
    },
    time::Duration,
};

use bit_set::BitSet;
use clap::ValueEnum;
use rayon::iter::{ParallelBridge, ParallelIterator};
use tokio::{runtime::Runtime, sync::oneshot, time::timeout as tktimeout};

use crate::{
    bounds::{state_bounds, Bound},
    canonize::CanonizeMode,
    kernels::KernelMode,
    matches::Matches,
    memoize::{Cache, MemoizeMode},
    molecule::Molecule,
    pathway::{generate_pathway, PathwayStep},
    state::State,
    utils::connected_components_under_edges,
};

/// Parallelization strategy for the recursive search phase.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum ParallelMode {
    /// No parallelism.
    None,
    /// Create a task pool from the recursion's first level only.
    DepthOne,
    /// Spawn a new thread at every recursive call.
    Always,
}

/// Compute assembly depth; see
/// [Pagel et al. (2024)](https://arxiv.org/abs/2409.05993).
///
/// Note: This function is currently very (unusably) slow. It primarily exists
/// in this crate as a proof of concept.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::assembly::depth;
/// use assembly_theory::loader::parse_molfile_str;
///
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let benzene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly depth.
/// assert_eq!(depth(&benzene), 3);
/// # Ok(())
/// # }
/// ```
pub fn depth(mol: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in mol.partitions().unwrap() {
        let l = if left.is_basic_unit() {
            0
        } else {
            depth(&left)
        };

        let r = if right.is_basic_unit() {
            0
        } else {
            depth(&right)
        };

        ix = ix.min(l.max(r) + 1)
    }
    ix
}

/// Determine the fragments produced from the given assembly state by removing
/// the given pair of edge-disjoint, isomorphic subgraphs and then adding one
/// back; return `None` if not possible.
fn fragments(mol: &Molecule, state: &[BitSet], h1: &BitSet, h2: &BitSet) -> Option<Vec<BitSet>> {
    // Attempt to find fragments f1 and f2 containing h1 and h2, respectively;
    // if either do not exist, exit without further fragmentation.
    let f1 = state.iter().enumerate().find(|(_, c)| h1.is_subset(c));
    let f2 = state.iter().enumerate().find(|(_, c)| h2.is_subset(c));
    let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
        return None;
    };

    let mut fragments = state.to_owned();

    // If the same fragment f1 (== f2) contains both h1 and h2, replace this
    // one fragment f1 with all connected components comprising f1 - (h1 U h2).
    // Otherwise, replace fragments f1 and f2 with all connected components
    // comprising f1 - h1 and f2 - h2, respectively.
    if i1 == i2 {
        let mut union = h1.clone();
        union.union_with(h2);
        let mut difference = f1.clone();
        difference.difference_with(&union);
        let c = connected_components_under_edges(mol.graph(), &difference);
        fragments.extend(c);
        fragments.swap_remove(i1);
    } else {
        let mut diff1 = f1.clone();
        diff1.difference_with(h1);
        let c1 = connected_components_under_edges(mol.graph(), &diff1);
        fragments.extend(c1);

        let mut diff2 = f2.clone();
        diff2.difference_with(h2);
        let c2 = connected_components_under_edges(mol.graph(), &diff2);
        fragments.extend(c2);

        fragments.swap_remove(i1.max(i2));
        fragments.swap_remove(i1.min(i2));
    }

    // Drop any singleton fragments, add h1 as a fragment, and return.
    fragments.retain(|i| i.len() > 1);
    fragments.push(h1.clone());
    Some(fragments)
}

/// Shared state for collecting one or more optimal assembly pathways.
///
/// When `max_pathways` is `None`, behaves as a single-best-pathway tracker
/// (backward-compatible). When `Some(n)`, collects up to `n` distinct
/// pathways that achieve the optimal assembly index, deduplicated by the
/// *set* of match indices used (ignoring removal order).
pub struct PathwayCollector {
    /// Collected pathways: each is (match_sequence, leaf_fragments).
    pub pathways: Vec<(Vec<usize>, Vec<BitSet>)>,
    /// Set of canonical match-set keys already collected, for dedup.
    seen_keys: HashSet<Vec<usize>>,
    /// Maximum number of pathways to collect, or `None` for single-best mode.
    max_pathways: Option<usize>,
}

impl PathwayCollector {
    /// Create a new collector. `None` = single-best mode; `Some(n)` = collect
    /// up to `n` distinct alternative pathways.
    pub fn new(max_pathways: Option<usize>) -> Self {
        Self {
            pathways: Vec::new(),
            seen_keys: HashSet::new(),
            max_pathways,
        }
    }

    /// Compute a canonical key for deduplication: the sorted set of match indices.
    /// Two decompositions that remove the same set of matches (in any order)
    /// produce identical pathways, so we deduplicate by sorted match set.
    fn match_set_key(path: &[usize]) -> Vec<usize> {
        let mut key = path.to_vec();
        key.sort();
        key
    }
}

/// Recursive helper for [`index_search`], only public for benchmarking.
///
/// Inputs:
/// - `mol`: The molecule whose assembly index is being calculated.
/// - `matches`: Structural information about the molecule's matched fragments.
/// - `state`: The current assembly state.
/// - `best_index`: The smallest assembly index for all assembly states so far.
/// - `bounds`: The list of bounding strategies to apply.
/// - `cache`: Memoization cache storing previously searched assembly states.
/// - `parallel_mode`: The parallelism mode for this state's match iteration.
/// - `best_pathway`: When `Some`, tracks the `(match_ix_sequence, leaf_fragments)`
///    for the best decomposition path found so far. When `None`, pathway tracking
///    is skipped entirely.
/// - `path_so_far`: The sequence of match indices removed from root to this state.
///
/// Returns, from this assembly state and any of its descendents:
/// - `usize`: An updated upper bound on the assembly index. (Note: If this
///   state is pruned by bounds or deemed redundant by memoization, then the
///   upper bound returned is unchanged.)
/// - `usize`: The number of assembly states searched.
pub fn recurse_index_search(
    mol: &Molecule,
    matches: &Matches,
    state: &State,
    best_index: Arc<AtomicUsize>,
    bounds: &[Bound],
    cache: &mut Cache,
    parallel_mode: ParallelMode,
    best_pathway: Option<&Arc<Mutex<PathwayCollector>>>,
    path_so_far: &[usize],
) -> (usize, usize) {
    // If any bounds would prune this assembly state or if memoization is
    // enabled and this assembly state is preempted by the cached state, halt.
    if state_bounds(mol, state, best_index.load(Relaxed), bounds) || cache.memoize_state(mol, state)
    {
        return (state.index(), 1);
    }

    // Generate a list of matches (i.e., pairs of edge-disjoint, isomorphic
    // fragments) to remove from this state.
    let (intermediate_frags, matches_to_remove): (Vec<BitSet>, Vec<usize>) =
        matches.matches_to_remove(mol, state, best_index.load(Relaxed), bounds);

    // Keep track of the best assembly index found in any of this assembly
    // state's children and the number of states searched, including this one.
    let best_child_index = AtomicUsize::from(state.index());
    let states_searched = AtomicUsize::from(1);

    // Define a closure that handles recursing to a new assembly state based on
    // the given match.
    let recurse_on_match = |i: usize, match_ix: usize| {
        let (h1, h2) = matches.match_fragments(match_ix);

        if let Some(fragments) = fragments(mol, &intermediate_frags, h1, h2) {
            // If using depth-one parallelism, all descendant states should be
            // computed serially.
            let new_parallel = if parallel_mode == ParallelMode::DepthOne {
                ParallelMode::None
            } else {
                parallel_mode
            };

            // Build the extended path for this child.
            let child_path: Vec<usize> = if best_pathway.is_some() {
                let mut p = path_so_far.to_vec();
                p.push(match_ix);
                p
            } else {
                Vec::new()
            };

            let new_state = state.update(fragments, i, match_ix, h1.len());

            // If this is a leaf state (no more matches can be removed) and we
            // are tracking pathways, check if this is the new best.
            if best_pathway.is_some() {
                let child_index = new_state.index();
                if child_index < best_index.load(Relaxed) {
                    // This will be checked again after recursion, but for leaf
                    // states we capture here.
                }
            }

            // Recurse using the remaining matches and updated fragments.
            let (child_index, child_states_searched) = recurse_index_search(
                mol,
                matches,
                &new_state,
                best_index.clone(),
                bounds,
                &mut cache.clone(),
                new_parallel,
                best_pathway,
                &child_path,
            );

            // Update the best assembly indices (across children states and the
            // entire search) and the number of descendant states searched.
            let prev_best = best_child_index.fetch_min(child_index, Relaxed);
            let new_global = best_index.fetch_min(child_index, Relaxed);

            // If this child achieved a new global best (or tied for the best
            // when collecting alternatives) and we are tracking pathways,
            // update the shared pathway data.
            //
            // LEAF-ONLY GUARD: Only record when new_state.index() == child_index,
            // i.e., no deeper recursion improved the index. This ensures we
            // capture actual leaf decompositions, not intermediate ancestors
            // that merely propagated a descendant's result.
            if let Some(bp) = best_pathway {
                if new_state.index() == child_index {
                    let current_best = best_index.load(Relaxed);
                    let is_new_best = child_index < prev_best.min(new_global);
                    let guard_ref = bp.lock().unwrap();
                    let is_alternative = guard_ref.max_pathways.is_some()
                        && child_index == current_best
                        && !is_new_best;
                    drop(guard_ref);

                    if is_new_best || is_alternative {
                        let mut guard = bp.lock().unwrap();
                        let current_best = best_index.load(Relaxed);

                        // DUAL-MODE UPDATE:
                        //   New best → clear all previously collected pathways
                        //              and start fresh with this one.
                        //   Tied best → add as an alternative if distinct and
                        //               under the collection cap.
                        if is_new_best && child_index <= current_best {
                            guard.pathways.clear();
                            guard.seen_keys.clear();
                            let key = PathwayCollector::match_set_key(&child_path);
                            guard.seen_keys.insert(key);
                            guard
                                .pathways
                                .push((child_path, new_state.fragments().clone()));
                        } else if child_index == current_best {
                            // Tied with current best: add if distinct and under cap.
                            let key = PathwayCollector::match_set_key(&child_path);
                            let at_cap = guard
                                .max_pathways
                                .map_or(true, |n| guard.pathways.len() >= n);
                            if !at_cap && !guard.seen_keys.contains(&key) {
                                guard.seen_keys.insert(key);
                                guard
                                    .pathways
                                    .push((child_path, new_state.fragments().clone()));
                            }
                        }
                    }
                }
            }

            states_searched.fetch_add(child_states_searched, Relaxed);
        }
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    if parallel_mode == ParallelMode::None {
        matches_to_remove
            .iter()
            .enumerate()
            .for_each(|(i, match_ix)| recurse_on_match(i, *match_ix));
    } else {
        matches_to_remove
            .iter()
            .enumerate()
            .par_bridge()
            .for_each(|(i, match_ix)| recurse_on_match(i, *match_ix));
    }

    (
        best_child_index.load(Relaxed),
        states_searched.load(Relaxed),
    )
}

/// Compute a molecule's assembly index and related information using a
/// top-down recursive algorithm, parameterized by the specified options.
///
/// If `timeout` is `None`, run until the assembly index is found. Otherwise,
/// stop after `timeout` milliseconds and return the best upper bound on the
/// assembly index found so far.
///
/// See [`CanonizeMode`], [`ParallelMode`], [`KernelMode`], and [`Bound`] for
/// details on how to customize the algorithm. Notably, bounds are applied in
/// the order they appear in the `bounds` slice. It is generally better to
/// provide bounds that are quick to compute first.
///
/// The results returned are:
/// - The molecule's `u32` assembly index (or an upper bound if timed out).
/// - The molecule's `u32` number of edge-disjoint isomorphic subgraph pairs.
/// - The `usize` total number of assembly [`State`]s searched if search
///   completes, and `None` otherwise (i.e., if search timed out).
/// - The assembly pathway as a `Vec<PathwayStep>` if `extract_pathway` is
///   `true` and the search completes, and `None` otherwise.
/// - The raw decomposition `(match_sequence, remnants)` if `extract_pathway`
///   is `true` and the search completes, and `None` otherwise.
/// - Alternative assembly pathways as a `Vec<Vec<PathwayStep>>` if
///   `max_pathways` is `Some` and the search completes, and `None` otherwise.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::{
///     assembly::{index_search, ParallelMode},
///     bounds::Bound,
///     canonize::CanonizeMode,
///     kernels::KernelMode,
///     loader::parse_molfile_str,
///     memoize::MemoizeMode,
/// };
///
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly index without parallelism, memoization,
/// // kernelization, or bounds.
/// let (slow_index, _, _, _, _, _, _) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     MemoizeMode::None,
///     KernelMode::None,
///     &[],
///     false,
///     None,
/// );
///
/// // Compute the molecule's assembly index with parallelism, memoization, and
/// // some bounds.
/// let (fast_index, _, _, _, _, _, _) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::DepthOne,
///     MemoizeMode::CanonIndex,
///     KernelMode::None,
///     &[Bound::Log, Bound::Int],
///     false,
///     None,
/// );
///
/// // Limit search to 1 ms, which should time out.
/// let (index_bound, _, states_searched, _, _, _, _) = index_search(
///     &anthracene,
///     Some(1),
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     MemoizeMode::None,
///     KernelMode::None,
///     &[],
///     false,
///     None,
/// );
///
/// assert_eq!(slow_index, 6);
/// assert_eq!(fast_index, 6);
/// assert!(index_bound >= fast_index && states_searched == None);
/// # Ok(())
/// # }
/// ```
pub fn index_search(
    mol: &Molecule,
    timeout: Option<u64>,
    canonize_mode: CanonizeMode,
    parallel_mode: ParallelMode,
    memoize_mode: MemoizeMode,
    kernel_mode: KernelMode,
    bounds: &[Bound],
    extract_pathway: bool,
    max_pathways: Option<usize>,
) -> (
    u32,
    u32,
    Option<usize>,
    Option<Vec<PathwayStep>>,
    Option<(Vec<usize>, Vec<BitSet>)>,
    Option<Vec<Vec<PathwayStep>>>,
    Option<Vec<Vec<usize>>>,
) {
    // Catch not-yet-implemented modes.
    if kernel_mode != KernelMode::None {
        panic!("The chosen --kernel mode is not implemented yet!")
    }

    // Create the initial assembly state and memoization cache.
    let state = State::new(mol);
    let mut cache = Cache::new(memoize_mode, canonize_mode);

    // Enumerate matches (i.e., pairs of edge-disjoint isomorphic fragments).
    let matches = Matches::new(mol, canonize_mode);

    // Use an `Arc` to track the best assembly index across parallel threads.
    let best_index = Arc::new(AtomicUsize::from(mol.graph().edge_count() - 1));

    // Optionally create shared pathway tracking state.
    let best_pathway = if extract_pathway {
        Some(Arc::new(Mutex::new(PathwayCollector::new(max_pathways))))
    } else {
        None
    };

    // Search for the shortest assembly pathway recursively.
    if let Some(timeout) = timeout {
        // If a timeout is provided, we will search within an asynchronous task
        // that can be interrupted after the specified duration (see below). To
        // avoid subsequent scope issues, make copies of various variables.
        let best_index_copy = best_index.clone();
        let best_pathway_copy = best_pathway.clone();
        let parse_mol = mol.clone();
        let mol = mol.clone();
        let bounds = bounds.to_vec();
        let num_matches = matches.len();

        // Search within a dedicated asynchronous runtime.
        let rt = Runtime::new().unwrap();
        let result = rt.block_on(async {
            let (send, recv) = oneshot::channel();
            rayon::spawn(move || {
                let _ = send.send(recurse_index_search(
                    &mol,
                    &matches,
                    &state,
                    best_index_copy,
                    &bounds,
                    &mut cache,
                    parallel_mode,
                    best_pathway_copy.as_ref(),
                    &[],
                ));
            });
            tktimeout(Duration::from_millis(timeout), recv).await
        });

        // If the search completes before the timeout, return the true assembly
        // index. Otherwise, return the best upper bound on the assembly index
        // found before timing out.
        let (index, states_searched) = match result {
            Ok(Ok((index, states_searched))) => (index, Some(states_searched)),
            Err(_) => (best_index.load(Relaxed), None),
            _ => panic!("An unexpected error occurred in async index_search"),
        };

        // Reconstruct the pathway(s) if tracking was enabled and search completed.
        let (pathway, decomposition, alternative_pathways, alternative_decompositions) =
            match (&best_pathway, &states_searched) {
                (Some(bp), Some(_)) => {
                    let guard = bp.lock().unwrap();
                    if guard.pathways.is_empty() {
                        (None, None, None, None)
                    } else {
                        let mol_for_pathway = parse_mol.clone();
                        let matches_for_pathway = Matches::new(&mol_for_pathway, canonize_mode);
                        let (ref match_seq, ref remnants) = guard.pathways[0];
                        let primary = generate_pathway(
                            &mol_for_pathway,
                            &matches_for_pathway,
                            match_seq,
                            remnants,
                        );
                        let decomp = (match_seq.clone(), remnants.clone());

                        let (alt_steps, alt_seqs): (
                            Option<Vec<Vec<PathwayStep>>>,
                            Option<Vec<Vec<usize>>>,
                        ) = if guard.pathways.len() > 1 {
                            let (steps, seqs): (Vec<_>, Vec<_>) = guard.pathways[1..]
                                .iter()
                                .map(|(ms, rem)| {
                                    (
                                        generate_pathway(
                                            &mol_for_pathway,
                                            &matches_for_pathway,
                                            ms,
                                            rem,
                                        ),
                                        ms.clone(),
                                    )
                                })
                                .unzip();
                            (Some(steps), Some(seqs))
                        } else {
                            (None, None)
                        };

                        (Some(primary), Some(decomp), alt_steps, alt_seqs)
                    }
                }
                _ => (None, None, None, None),
            };

        (
            index as u32,
            num_matches as u32,
            states_searched,
            pathway,
            decomposition,
            alternative_pathways,
            alternative_decompositions,
        )
    } else {
        // Otherwise, if no timeout is provided, run the search normally.
        let (index, states_searched) = recurse_index_search(
            mol,
            &matches,
            &state,
            best_index,
            bounds,
            &mut cache,
            parallel_mode,
            best_pathway.as_ref(),
            &[],
        );

        // Reconstruct the pathway(s) if tracking was enabled.
        let (pathway, decomposition, alternative_pathways, alternative_decompositions) =
            match best_pathway {
                Some(bp) => {
                    let guard = bp.lock().unwrap();
                    if guard.pathways.is_empty() {
                        (None, None, None, None)
                    } else {
                        let (ref match_seq, ref remnants) = guard.pathways[0];
                        let primary = generate_pathway(mol, &matches, match_seq, remnants);
                        let decomp = (match_seq.clone(), remnants.clone());

                        let (alt_steps, alt_seqs): (
                            Option<Vec<Vec<PathwayStep>>>,
                            Option<Vec<Vec<usize>>>,
                        ) = if guard.pathways.len() > 1 {
                            let (steps, seqs): (Vec<_>, Vec<_>) = guard.pathways[1..]
                                .iter()
                                .map(|(ms, rem)| {
                                    (generate_pathway(mol, &matches, ms, rem), ms.clone())
                                })
                                .unzip();
                            (Some(steps), Some(seqs))
                        } else {
                            (None, None)
                        };

                        (Some(primary), Some(decomp), alt_steps, alt_seqs)
                    }
                }
                None => (None, None, None, None),
            };

        (
            index as u32,
            matches.len() as u32,
            Some(states_searched),
            pathway,
            decomposition,
            alternative_pathways,
            alternative_decompositions,
        )
    }
}

/// Compute a molecule's assembly index using an efficient default strategy.
///
/// To customize assembly index calculation beyond the default strategy, see
/// [`index_search`].
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::assembly::index;
/// use assembly_theory::loader::parse_molfile_str;
///
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly index.
/// assert_eq!(index(&anthracene), 6);
/// # Ok(())
/// # }
/// ```
pub fn index(mol: &Molecule) -> u32 {
    index_search(
        mol,
        None,
        CanonizeMode::TreeNauty,
        ParallelMode::DepthOne,
        MemoizeMode::CanonIndex,
        KernelMode::None,
        &[Bound::Int, Bound::MatchableEdges],
        false,
        None,
    )
    .0
}
