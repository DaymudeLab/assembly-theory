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
        Arc,
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
    pathway::Pathway,
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

/// Find all fragments in the given assembly state after removing the given
/// pair of edge-disjoint, isomorphic subgraphs and then adding one back.
fn fragments(mol: &Molecule, state: &[BitSet], h1: &BitSet, h2: &BitSet) -> Vec<BitSet> {
    // Find fragments f1 and f2 containing h1 and h2, respectively.
    let (i1, f1) = state
        .iter()
        .enumerate()
        .find(|(_, c)| h1.is_subset(c))
        .unwrap();
    let (i2, f2) = state
        .iter()
        .enumerate()
        .find(|(_, c)| h2.is_subset(c))
        .unwrap();

    let mut fragments = state.to_owned();

    // If the same fragment f1 (== f2) contains both h1 and h2, replace this
    // one fragment f1 with all connected components comprising f1 - (h1 U h2).
    // Otherwise, replace fragments f1 and f2 with all connected components
    // comprising f1 - h1 and f2 - h2, respectively.
    if i1 == i2 {
        let mut union = h1.clone();
        union.union_with(h2);
        let mut diff = f1.clone();
        diff.difference_with(&union);
        let c = connected_components_under_edges(mol.graph(), &diff);
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
    fragments.retain(|i| i.count() > 1);
    fragments.push(h1.clone());
    fragments
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
/// - `collect_removal_orders`: `true` iff match removal orders should be
///   collected for later assembly pathway reconstruction.
///
/// Returns, from this assembly state and any of its descendents:
/// - `usize`: An updated upper bound on the assembly index. (Note: If this
///   state is pruned by bounds or deemed redundant by memoization, then the
///   upper bound returned is unchanged.)
/// - `usize`: The number of assembly states searched.
/// - `Option<HashSet<Vec<usize>>>`: A set of assembly states' match removal
///   orders that attain the updated assembly index upper bound, or `None` if
///   `collect_removal_orders == false`.
pub fn recurse_index_search(
    mol: &Molecule,
    matches: &Matches,
    state: &State,
    best_index: Arc<AtomicUsize>,
    bounds: &[Bound],
    cache: &mut Cache,
    parallel_mode: ParallelMode,
    collect_removal_orders: bool,
) -> (usize, usize, Option<HashSet<Vec<usize>>>) {
    // If any bounds would prune this assembly state or if memoization is
    // enabled and this assembly state is preempted by the cached state, halt.
    if state_bounds(mol, state, best_index.load(Relaxed), bounds) || cache.memoize_state(mol, state)
    {
        return (
            state.index(),
            1,
            if collect_removal_orders {
                Some(HashSet::from([state.removal_order().clone()]))
            } else {
                None
            },
        );
    }

    // Generate a list of matches (i.e., pairs of edge-disjoint, isomorphic
    // fragments) to remove from this state.
    let (intermediate_frags, matches_to_remove): (Vec<BitSet>, Vec<usize>) =
        matches.matches_to_remove(mol, state, best_index.load(Relaxed), bounds);

    // Define a closure that handles recursing to a new assembly state based on
    // the given match.
    let recurse_on_match = |match_ix: usize| -> (usize, usize, Option<HashSet<Vec<usize>>>) {
        let (h1, h2) = matches.match_fragments(match_ix);
        let fragments = fragments(mol, &intermediate_frags, h1, h2);

        // If using depth-one parallelism, all descendant states should be
        // computed serially.
        let new_parallel = if parallel_mode == ParallelMode::DepthOne {
            ParallelMode::None
        } else {
            parallel_mode
        };

        // Recurse using the remaining matches and updated fragments.
        recurse_index_search(
            mol,
            matches,
            &state.update(fragments, match_ix, h1.count()),
            best_index.clone(),
            bounds,
            &mut cache.clone(),
            new_parallel,
            collect_removal_orders,
        )
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    let results: Vec<(usize, usize, Option<HashSet<Vec<usize>>>)> =
        if parallel_mode == ParallelMode::None {
            matches_to_remove
                .iter()
                .map(|match_ix| recurse_on_match(*match_ix))
                .collect()
        } else {
            matches_to_remove
                .iter()
                .par_bridge()
                .map(|match_ix| recurse_on_match(*match_ix))
                .collect()
        };

    // Compute the best assembly index bound and the total number of descendant
    // states searched across all children assembly states.
    let best_child_index = results.iter().map(|r| r.0).min().unwrap_or(state.index());
    let states_searched: usize = results.iter().map(|r| r.1).sum();

    // Update the globally best assembly index bound found so far.
    best_index.fetch_min(best_child_index, Relaxed);

    // Collect all descendant match removal orders attaining the best assembly
    // index bound across children assembly states.
    if collect_removal_orders {
        let mut best_removal_orders = if state.index() == best_child_index {
            HashSet::from([state.removal_order().clone()])
        } else {
            HashSet::<Vec<usize>>::new()
        };
        for r in results {
            if r.0 == best_child_index {
                best_removal_orders.extend(r.2.unwrap());
            }
        }
        (
            best_child_index,
            states_searched + 1,
            Some(best_removal_orders),
        )
    } else {
        (best_child_index, states_searched + 1, None)
    }
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
/// If `max_pathways` is `None`, skip minimum assembly pathway reconstruction.
/// Otherwise, reconstruct all such pathways from the completed search results
/// if `max_pathways == 0` or at most `max_pathways` such pathways otherwise.
///
/// The results returned are:
/// - The molecule's `u32` assembly index (or an upper bound if timed out).
/// - The molecule's `u32` number of edge-disjoint isomorphic subgraph pairs.
/// - The `usize` total number of assembly [`State`]s searched if search
///   completes, and `None` otherwise (i.e., if search timed out).
/// - A vector of of minimum assembly [`Pathway`]s.
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
/// let (slow_index, _, _, _) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     MemoizeMode::None,
///     KernelMode::None,
///     &[],
///     None,
/// );
/// assert_eq!(slow_index, 6);
///
/// // Compute the molecule's assembly index with parallelism, memoization, and
/// // some bounds.
/// let (fast_index, _, _, _) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::DepthOne,
///     MemoizeMode::CanonIndex,
///     KernelMode::None,
///     &[Bound::Log, Bound::Int],
///     None,
/// );
/// assert_eq!(fast_index, 6);
///
/// // Limit search to 1 ms, which should time out.
/// let (index_bound, _, states_searched, _) = index_search(
///     &anthracene,
///     Some(1),
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     MemoizeMode::None,
///     KernelMode::None,
///     &[],
///     None,
/// );
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
    max_pathways: Option<usize>,
) -> (u32, u32, Option<usize>, Vec<Pathway>) {
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

    // Search for the shortest assembly pathway recursively.
    let (index, states_searched, removal_orders) = if let Some(timeout) = timeout {
        // If a timeout is provided, we will search within an asynchronous task
        // that can be interrupted after the specified duration (see below). To
        // avoid subsequent scope issues, make copies of various variables.
        let mol = mol.clone();
        let matches = matches.clone();
        let best_index_copy = best_index.clone();
        let bounds = bounds.to_vec();

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
                    max_pathways.is_some(),
                ));
            });
            tktimeout(Duration::from_millis(timeout), recv).await
        });

        // If the search completes before the timeout, return the true assembly
        // index. Otherwise, return the best upper bound on the assembly index
        // found before timing out.
        match result {
            Ok(Ok((index, states_searched, removal_orders))) => {
                (index, Some(states_searched), removal_orders)
            }
            Err(_) => (best_index.load(Relaxed), None, None),
            _ => panic!("An unexpected error occurred in async index_search"),
        }
    } else {
        // Otherwise, if no timeout is provided, run the search normally.
        let (index, states_searched, removal_orders) = recurse_index_search(
            mol,
            &matches,
            &state,
            best_index,
            bounds,
            &mut cache,
            parallel_mode,
            max_pathways.is_some(),
        );
        (index, Some(states_searched), removal_orders)
    };

    // Generate at most the given number of minimum assembly pathways from the
    // returned match removal orders.
    let mut pathways = Vec::<Pathway>::new();
    if let Some(max_pathways) = max_pathways {
        if let Some(mut removal_orders) = removal_orders {
            for (ix, removal_order) in removal_orders.drain().enumerate() {
                if max_pathways > 0 && ix >= max_pathways {
                    break;
                } else {
                    pathways.push(Pathway::new(mol, &matches, &removal_order));
                }
            }
        }
    }

    (
        index as u32,
        matches.len() as u32,
        states_searched,
        pathways,
    )
}

/// Compute a molecule's assembly index using an efficient default strategy.
///
/// To customize assembly index calculation beyond the default strategy and/or
/// reconstruct minimum assembly pathways, see [`index_search`].
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
        None, // Run without a timeout.
        CanonizeMode::TreeNauty,
        ParallelMode::DepthOne,
        MemoizeMode::CanonIndex,
        KernelMode::None,
        &[Bound::Int, Bound::MatchableEdges],
        None, // Disable pathway reconstruction.
    )
    .0
}
