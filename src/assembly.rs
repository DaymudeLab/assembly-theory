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
//! To customize assembly index calculation beyond the default strategy, set a
//! timeout for long-running calculations, and/or reconstruct minimum assembly
//! pathways, see [`index_search`].

use std::{
    collections::BTreeSet,
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
/// - `max_pathways`: The number of distinct match removal orders attaining the best assembly index
///   bound to collect; `0` for all, `x > 0` for at most `x`, or `None` for none.
///
/// Returns, from this assembly state and any of its descendents:
/// - `usize`: An updated upper bound on the assembly index. (Note: If this
///   state is pruned by bounds or deemed redundant by memoization, then the
///   upper bound returned is unchanged.)
/// - `usize`: The number of assembly states searched.
/// - `BTreeSet<Vec<usize>>`: A set of assembly states' match removal orders that attain the
///   updated assembly index upper bound. Nonempty iff `max_pathways` is not `None`.
#[allow(clippy::too_many_arguments)]
pub fn recurse_index_search(
    mol: &Molecule,
    matches: &Matches,
    state: &State,
    best_index: Arc<AtomicUsize>,
    bounds: &[Bound],
    cache: &mut Cache,
    parallel_mode: ParallelMode,
    max_pathways: Option<usize>,
) -> (usize, usize, BTreeSet<Vec<usize>>) {
    // If any bounds would prune this assembly state or if memoization is
    // enabled and this assembly state is preempted by the cached state, halt.
    if state_bounds(mol, state, best_index.load(Relaxed), bounds) || cache.memoize_state(mol, state)
    {
        return (
            state.index(),
            1,
            if max_pathways.is_some() {
                BTreeSet::from([state.removal_order().clone()])
            } else {
                BTreeSet::<Vec<usize>>::new()
            },
        );
    }

    // Generate a list of matches (i.e., pairs of edge-disjoint, isomorphic
    // fragments) to remove from this state.
    let (intermediate_frags, matches_to_remove): (Vec<BitSet>, Vec<usize>) =
        matches.matches_to_remove(mol, state, best_index.load(Relaxed), bounds);

    // Define a closure that handles recursing to a new assembly state based on
    // the given match.
    let recurse_on_match = |match_ix: usize| -> (usize, usize, BTreeSet<Vec<usize>>) {
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
            max_pathways,
        )
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    let results: Vec<(usize, usize, BTreeSet<Vec<usize>>)> = if parallel_mode == ParallelMode::None
    {
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
    let mut best_removal_orders = BTreeSet::<Vec<usize>>::new();
    if let Some(max_pathways) = max_pathways {
        if state.index() == best_child_index {
            best_removal_orders.insert(state.removal_order().clone());
        }
        for mut r in results {
            if r.0 == best_child_index {
                if max_pathways == 0 || best_removal_orders.len() + r.2.len() <= max_pathways {
                    best_removal_orders.extend(r.2)
                } else {
                    while best_removal_orders.len() < max_pathways {
                        best_removal_orders.insert(r.2.pop_first().unwrap());
                    }
                    break; // Collected desired number of removal orders; no need to check more.
                }
            }
        }
    }
    (best_child_index, states_searched + 1, best_removal_orders)
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
/// - A vector of the molecule's minimum assembly [`Pathway`]s that is nonempty
///   if and only if search did not time out and `max_pathways` is not `None`.
///
/// # Examples
///
/// ## Customizing Assembly Index Search Options
///
/// The two examples below show different customizations of the assembly index
/// search options, enabling/disabling parallelism, memoization, and bounding.
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
/// // Calculate the molecule's assembly index without any search optimizations
/// // like parallelism, memoization, kernelization, or bounding.
/// let (slow_index, num_matches, states_searched, _) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,      // Disable parallelism.
///     MemoizeMode::None,       // Disable memoization.
///     KernelMode::None,        // Disable kernelization.
///     &[],                     // Do not use any bounds.
///     None,
/// );
/// assert_eq!(slow_index, 6);
/// assert_eq!(num_matches, 466);
/// assert_eq!(states_searched, Some(133814));
///
/// // Repeat the calculation with parallelism, memoization, and some bounds.
/// let (fast_index, num_matches, states_searched, _) = index_search(
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
/// assert_eq!(num_matches, 466);
/// // Note: states_searched is now nondeterministic because of parallelism.
/// # Ok(())
/// # }
/// ```
///
/// ## Setting a Timeout for Long Calculations
///
/// Assembly index calculation can be slow on large molecules, even with search
/// optimizations enabled. If `timeout` is set, two outcomes are possible:
/// 1. Search completes before the timeout and everything behaves as usual.
/// 2. The timeout is reached before search completes, search is interrupted,
///    and the best upper bound on the assembly index found so far is returned.
///    The number of states searched and assembly pathways are not returned.
/// ```
/// # use std::{fs, path::PathBuf};
/// # use assembly_theory::{
/// #     assembly::{index_search, ParallelMode},
/// #     bounds::Bound,
/// #     canonize::CanonizeMode,
/// #     kernels::KernelMode,
/// #     loader::parse_molfile_str,
/// #     memoize::MemoizeMode,
/// # };
/// #
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// # let molfile = fs::read_to_string(path)?;
/// # let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
/// #
/// // Limit search to 10 s, which is plenty for search to complete as usual,
/// // even with some search optimizations disabled.
/// let (true_index, num_matches, states_searched, pathways) = index_search(
///     &anthracene,
///     Some(10000),               // Set a 10 s timeout.
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,        // Disable parallelism.
///     MemoizeMode::None,         // Disable memoization.
///     KernelMode::None,          // Disable kernelization.
///     &[Bound::Log, Bound::Int], // Enable some bounds.
///     Some(1),                   // Reconstuct one pathway.
/// );
/// assert_eq!(true_index, 6);
/// assert_eq!(num_matches, 466);
/// assert_eq!(states_searched, Some(2454));
/// assert_eq!(pathways.len(), 1);
///
/// // Limit search to 1 ms, which will time out without more optimizations.
/// let (index_bound, num_matches, states_searched, pathways) = index_search(
///     &anthracene,
///     Some(1),                   // Set a 1 ms timeout.
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,        // Disable parallelism.
///     MemoizeMode::None,         // Disable memoization.
///     KernelMode::None,          // Disable kernelization.
///     &[Bound::Log, Bound::Int], // Enable some bounds.
///     Some(1),                   // Reconstuct one pathway.
/// );
/// assert!(index_bound >= true_index); // Index is an upper bound on timeout.
/// assert_eq!(num_matches, 466);       // The number of matches is unaffected.
/// assert_eq!(states_searched, None);  // States searched is not computed on timeout.
/// assert_eq!(pathways.len(), 0);      // Pathways are not computed on timeout.
/// # Ok(())
/// # }
/// ```
///
/// ## Reconstructing Minimum Assembly Pathways
///
/// An assembly [`Pathway`] records the fragment join operations producing the
/// target molecule. Set `max_pathways` to `Some(0)` to reconstruct all minimum
/// assembly pathways discovered by the search process, `Some(n)` for at most
/// the first `n` such pathways, or `None` to disable pathway reconstruction.
/// See [`Pathway::dag`] for details on the pathway data structure.
///
/// Note that the definition of `Some(0)`'s behavior is very carefully worded:
/// it finds all minimum assembly pathways *discovered by the search process*.
/// This is not to be taken as *all possible minimum assembly pathways*. When
/// search optimizations like parallelism, memoization, and bounding are used,
/// search prunes pathways that cannot yield a better assembly index. This may
/// include pathways that yield an *equal* (i.e., also minimum) assembly index.
/// Once pruned, these pathways will not survive to pathway reconstruction.
/// ```
/// # use std::{fs, path::PathBuf};
/// # use assembly_theory::{
/// #     assembly::{index_search, ParallelMode},
/// #     bounds::Bound,
/// #     canonize::CanonizeMode,
/// #     kernels::KernelMode,
/// #     loader::parse_molfile_str,
/// #     memoize::MemoizeMode,
/// # };
/// #
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// # let molfile = fs::read_to_string(path)?;
/// # let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
/// #
/// // Reconstruct a minimum assembly pathway discovered during search.
/// let (_, _, _, pathways) = index_search(
///     &anthracene,
///     None,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,        // Disable parallelism for determinism.
///     MemoizeMode::CanonIndex,
///     KernelMode::None,
///     &[Bound::Log, Bound::Int],
///     Some(1),                   // Reconstuct one pathway.
/// );
/// let pathway_str = r#"digraph {
///     0 [ label = "{14}" ]
///     1 [ label = "{15}" ]
///     2 [ label = "{14, 15}" ]
///     3 [ label = "{7, 14, 15}" ]
///     4 [ label = "{6, 7, 14, 15}" ]
///     5 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
///     6 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
///     7 [ label = "{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}" ]
///     0 -> 2 [ label = "{14}" ]
///     1 -> 2 [ label = "{15}" ]
///     0 -> 3 [ label = "{7}" ]
///     2 -> 3 [ label = "{14, 15}" ]
///     1 -> 4 [ label = "{6}" ]
///     3 -> 4 [ label = "{7, 14, 15}" ]
///     4 -> 5 [ label = "{6, 7, 14, 15}" ]
///     3 -> 5 [ label = "{3, 4, 5}" ]
///     2 -> 6 [ label = "{2, 13}" ]
///     5 -> 6 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
///     6 -> 7 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
///     5 -> 7 [ label = "{0, 1, 8, 9, 10, 11, 12}" ]
/// }
/// "#;
/// assert_eq!(format!("{}", pathways[0]), pathway_str);
/// # Ok(())
/// # }
/// ```
#[allow(clippy::too_many_arguments)]
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
                    max_pathways,
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
            Err(_) => (
                best_index.load(Relaxed),
                None,
                BTreeSet::<Vec<usize>>::new(),
            ),
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
            max_pathways,
        );
        (index, Some(states_searched), removal_orders)
    };

    // Generate at most the given number of minimum assembly pathways from the
    // returned match removal orders.
    let mut pathways = Vec::<Pathway>::new();
    if let Some(max_pathways) = max_pathways {
        for (ix, removal_order) in removal_orders.iter().enumerate() {
            if max_pathways > 0 && ix >= max_pathways {
                break;
            } else {
                pathways.push(Pathway::new(mol, &matches, removal_order));
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
/// To customize assembly index calculation beyond the default strategy, set a
/// timeout for long-running calculations, and/or reconstruct minimum assembly
/// pathways, see [`index_search`].
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
