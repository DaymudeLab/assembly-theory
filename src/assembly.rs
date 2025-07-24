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
    collections::HashMap,
    sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    },
};

use bit_set::BitSet;
use clap::ValueEnum;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    bounds::{bound_exceeded, Bound},
    canonize::{canonize, CanonizeMode, Labeling},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    kernels::KernelMode,
    molecule::Molecule,
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

/// Helper function for [`index_search`]; only public for benchmarking.
///
/// Return (1) a map of the given molecule's connected, non-induced subgraphs
/// with at most |E|/2 edges to their canonical labels and (2) all pairs of
/// non-overlapping, isomorphic subgraphs in the molecule, sorted to guarantee
/// deterministic iteration.
pub fn labels_matches(
    mol: &Molecule,
    enumerate_mode: EnumerateMode,
    canonize_mode: CanonizeMode,
) -> (HashMap<BitSet, Labeling>, Vec<(BitSet, BitSet)>) {
    // Enumerate all connected, non-induced subgraphs with at most |E|/2 edges
    // and bin them into isomorphism classes using canonization. Store these
    // subgraphs' canonical labels for later use in memoization (this way, we
    // only have to compute them once).
    let mut isomorphism_classes = HashMap::<Labeling, Vec<BitSet>>::new();
    let mut subgraph_labels = HashMap::<BitSet, Labeling>::new();
    for subgraph in enumerate_subgraphs(mol, enumerate_mode) {
        let label = canonize(mol, &subgraph, canonize_mode);
        isomorphism_classes
            .entry(label.clone())
            .and_modify(|bucket| bucket.push(subgraph.clone()))
            .or_insert(vec![subgraph.clone()]);
        subgraph_labels.insert(subgraph, label);
    }

    // In each isomorphism class, collect non-overlapping pairs of subgraphs.
    let mut matches = Vec::new();
    for bucket in isomorphism_classes.values() {
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

    (subgraph_labels, matches)
}

/// Determine the fragments produced by removing the given pair of duplicatable
/// subgraphs and then adding one back; return None if not possible.
fn fractures(
    mol: &Molecule,
    fragments: &[BitSet],
    h1: &BitSet,
    h2: &BitSet,
) -> Option<Vec<BitSet>> {
    // Attempt to find fragments f1 and f2 containing h1 and h2, respectively;
    // if either do not exist, exit without further fragmentation.
    let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
    let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));
    let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
        return None;
    };

    let mut fractures = fragments.to_owned();

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
        fractures.extend(c);
        fractures.swap_remove(i1);
    } else {
        let mut diff1 = f1.clone();
        diff1.difference_with(h1);
        let c1 = connected_components_under_edges(mol.graph(), &diff1);
        fractures.extend(c1);

        let mut diff2 = f2.clone();
        diff2.difference_with(h2);
        let c2 = connected_components_under_edges(mol.graph(), &diff2);
        fractures.extend(c2);

        fractures.swap_remove(i1.max(i2));
        fractures.swap_remove(i1.min(i2));
    }

    // Drop any singleton fragments, add h1 as a fragment, and return.
    fractures.retain(|i| i.len() > 1);
    fractures.push(h1.clone());
    Some(fractures)
}

/// Recursive helper for [`index_search`], only public for benchmarking.
///
/// Inputs:
/// - `mol`: The molecule whose assembly index is being calculated.
/// - `matches`: The remaining non-overlapping isomorphic subgraph pairs.
/// - `fragments`: TODO
/// - `state_index`: The assembly index of this assembly state.
/// - `best_index`: The smallest assembly index for all assembly states so far.
/// - `largest_remove`: An upper bound on the size of fragments that can be
///   removed from this or any descendant state.
/// - `bounds`: The list of bounding strategies to apply.
/// - `parallel_mode`: The parallelism mode to use for this state's children.
///
/// Returns, from this assembly state and any of its descendents:
/// - `usize`: The best assembly index found.
/// - `usize`: The number of assembly states searched.
#[allow(clippy::too_many_arguments)]
pub fn recurse_index_search(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    state_index: usize,
    best_index: Arc<AtomicUsize>,
    largest_remove: usize,
    bounds: &[Bound],
    parallel_mode: ParallelMode,
) -> (usize, usize) {
    // If any bounds would prune this assembly state, halt.
    if bound_exceeded(
        mol,
        fragments,
        state_index,
        best_index.load(Relaxed),
        largest_remove,
        bounds,
    ) {
        return (state_index, 1);
    }

    // Keep track of the best assembly index found in any of this assembly
    // state's children and the number of states searched, including this one.
    let best_child_index = AtomicUsize::from(state_index);
    let states_searched = AtomicUsize::from(1);

    // Define a closure that handles recursing to a new assembly state based on
    // the given (enumerated) pair of non-overlapping isomorphic subgraphs.
    let recurse_on_match = |i: usize, h1: &BitSet, h2: &BitSet| {
        if let Some(fractures) = fractures(mol, fragments, h1, h2) {
            // If using depth-one parallelism, all descendant states should be
            // computed serially.
            let new_parallel = if parallel_mode == ParallelMode::DepthOne {
                ParallelMode::None
            } else {
                parallel_mode
            };

            // Recurse using the remaining matches and updated fragments.
            let (child_index, child_states_searched) = recurse_index_search(
                mol,
                &matches[i + 1..],
                &fractures,
                state_index - h1.len() + 1,
                best_index.clone(),
                h1.len(),
                bounds,
                new_parallel,
            );

            // Update the best assembly indices (across children states and
            // the entire search).
            best_child_index.fetch_min(child_index, Relaxed);
            best_index.fetch_min(best_child_index.load(Relaxed), Relaxed);
            states_searched.fetch_add(child_states_searched, Relaxed);
        }
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    if parallel_mode == ParallelMode::None {
        matches
            .iter()
            .enumerate()
            .for_each(|(i, (h1, h2))| recurse_on_match(i, h1, h2));
    } else {
        matches
            .par_iter()
            .enumerate()
            .for_each(|(i, (h1, h2))| recurse_on_match(i, h1, h2));
    }

    (
        best_child_index.load(Relaxed),
        states_searched.load(Relaxed),
    )
}

/// Compute a molecule's assembly index and related information using a
/// top-down recursive algorithm, parameterized by the specified options.
///
/// See [`EnumerateMode`], [`CanonizeMode`], [`ParallelMode`], [`KernelMode`],
/// and [`Bound`] for details on how to customize the algorithm. Notably,
/// bounds are applied in the order they appear in the `bounds` slice. It is
/// generally better to provide bounds that are quick to compute first.
///
/// The results returned are:
/// - The molecule's `u32` assembly index.
/// - The molecule's `u32` number of non-overlapping isomorphic subgraph pairs.
/// - The `usize` total number of assembly states searched, where an assembly
///   state is a collection of fragments. Note that, depending on the algorithm
///   parameters used, some states may be searched/counted multiple times.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::{
///     assembly::{index_search, ParallelMode},
///     bounds::Bound,
///     canonize::CanonizeMode,
///     enumerate::EnumerateMode,
///     kernels::KernelMode,
///     loader::parse_molfile_str,
/// };
///
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly index without parallelism,
/// // kernelization, bounds, or memoization.
/// let (slow_index, _, _) = index_search(
///     &anthracene,
///     EnumerateMode::GrowErode,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     KernelMode::None,
///     &[],
///     false);
///
/// // Compute the molecule's assembly index with parallelism and some bounds.
/// let (fast_index, _, _) = index_search(
///     &anthracene,
///     EnumerateMode::GrowErode,
///     CanonizeMode::TreeNauty,
///     ParallelMode::DepthOne,
///     KernelMode::None,
///     &[Bound::Log, Bound::Int],
///     false);
///
/// assert_eq!(slow_index, 6);
/// assert_eq!(fast_index, 6);
/// # Ok(())
/// # }
/// ```
pub fn index_search(
    mol: &Molecule,
    enumerate_mode: EnumerateMode,
    canonize_mode: CanonizeMode,
    parallel_mode: ParallelMode,
    kernel_mode: KernelMode,
    bounds: &[Bound],
    memoize: bool,
) -> (u32, u32, usize) {
    // Catch not-yet-implemented modes.
    if kernel_mode != KernelMode::None {
        panic!("The chosen --kernel mode is not implemented yet!")
    }
    if memoize {
        panic!("--memoize is not implemented yet!")
    }

    // Enumerate non-overlapping isomorphic subgraph pairs.
    let (_, matches) = labels_matches(mol, enumerate_mode, canonize_mode);

    // Initialize the first fragment as the entire graph.
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Search for the shortest assembly pathway recursively.
    let edge_count = mol.graph().edge_count();
    let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
    let (index, states_searched) = recurse_index_search(
        mol,
        &matches,
        &[init],
        edge_count - 1,
        best_index,
        edge_count,
        bounds,
        parallel_mode,
    );

    (index as u32, matches.len() as u32, states_searched)
}

/// Compute a molecule's assembly index using an efficient default strategy.
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
        EnumerateMode::GrowErode,
        CanonizeMode::TreeNauty,
        ParallelMode::DepthOne,
        KernelMode::None,
        &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
        false,
    )
    .0
}
