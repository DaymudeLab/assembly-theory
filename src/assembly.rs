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
    canonize::{canonize, CanonizeMode},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    kernels::KernelMode,
    molecule::Molecule,
    utils::connected_components_under_edges,
};

/// Parallelization strategy for the search phase.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum ParallelMode {
    /// No parallelism.
    None,
    /// Create a task pool form the recursion's first level only.
    DepthOne,
    /// Spawn a new thread at every recursive call.
    Always,
}

static PARALLEL_MATCH_SIZE_THRESHOLD: usize = 100;

/// Computes assembly depth; see
/// [Pagel et al. (2024)](https://arxiv.org/abs/2409.05993).
pub fn assembly_depth(mol: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in mol.partitions().unwrap() {
        let l = if left.is_basic_unit() {
            0
        } else {
            assembly_depth(&left)
        };

        let r = if right.is_basic_unit() {
            0
        } else {
            assembly_depth(&right)
        };

        ix = ix.min(l.max(r) + 1)
    }
    ix
}

/// Return an iterator over all pairs of non-overlapping, isomorphic subgraphs
/// in the given molecule.
fn matches(
    mol: &Molecule,
    enumerate_mode: EnumerateMode,
    canonize_mode: CanonizeMode,
) -> impl Iterator<Item = (BitSet, BitSet)> {
    // Enumerate all connected, non-induced subgraphs and bin them into
    // isomorphism classes using canonization.
    let mut isomorphic_map = HashMap::<_, Vec<BitSet>>::new();
    for subgraph in enumerate_subgraphs(mol, enumerate_mode) {
        isomorphic_map
            .entry(canonize(mol, &subgraph, canonize_mode))
            .and_modify(|bucket| bucket.push(subgraph.clone()))
            .or_insert(vec![subgraph.clone()]);
    }

    // In each isomorphism class, enumerate non-overlapping pairs of subgraphs.
    let mut matches = Vec::new();
    for bucket in isomorphic_map.values() {
        for (i, first) in bucket.iter().enumerate() {
            for second in &bucket[i..] {
                if first.is_disjoint(second) {
                    matches.push((first.clone(), second.clone()));
                }
            }
        }
    }
    matches.into_iter()
}

/// Recursive helper for the serial version of index_search.
///
/// Inputs:
/// - `mol`: The molecule whose assembly index is being calculated.
/// - `matches`: The remaining non-overlapping isomorphic subgraph pairs.
/// - `fragments`: TODO
/// - `state_index`: The assembly index of this assembly state.
/// - `best_index`: The smallest assembly index for all assembly states so far.
/// - `largest_remove`: The size of the largest match removed so far.
/// - `bounds`: The list of bounding strategies to apply.
/// - `states_searched`: The number of assembly states searched so far.
#[allow(clippy::too_many_arguments)]
fn recurse_index_search_serial(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    state_index: usize,
    mut best_index: usize,
    largest_remove: usize,
    bounds: &[Bound],
    states_searched: &mut usize,
) -> usize {
    // Keep track of the best assembly index found in any of this assembly
    // state's children.
    let mut best_child_index = state_index;

    // Count this assembly state as searched.
    *states_searched += 1;

    // Check if any bounds apply; if so, halt this search branch.
    if bound_exceeded(
        mol,
        fragments,
        state_index,
        best_index,
        largest_remove,
        bounds,
    ) {
        return state_index;
    }

    // Search for duplicatable fragment
    for (i, (h1, h2)) in matches.iter().enumerate() {
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let largest_remove = h1.len();

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            continue;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));
        }

        fractures.retain(|i| i.len() > 1);
        fractures.push(h1.clone());

        best_child_index = best_child_index.min(recurse_index_search_serial(
            mol,
            &matches[i + 1..],
            &fractures,
            state_index - h1.len() + 1,
            best_index,
            largest_remove,
            bounds,
            states_searched,
        ));
        best_index = best_index.min(best_child_index);
    }

    best_child_index
}

/// Recursive helper for the parallel version of index_search.
///
/// Inputs:
/// - `mol`: The molecule whose assembly index is being calculated.
/// - `matches`: The remaining non-overlapping isomorphic subgraph pairs.
/// - `fragments`: TODO
/// - `state_index`: The assembly index of this assembly state.
/// - `best_index`: The smallest assembly index for all assembly states so far.
/// - `largest_remove`: The size of the largest match removed so far.
/// - `bounds`: The list of bounding strategies to apply.
/// - `states_searched`: The number of assembly states searched so far.
#[allow(clippy::too_many_arguments)]
fn recurse_index_search_parallel(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    state_index: usize,
    best_index: AtomicUsize,
    largest_remove: usize,
    bounds: &[Bound],
    states_searched: Arc<AtomicUsize>,
) -> usize {
    // Keep track of the best assembly index found in any of this assembly
    // state's children.
    let best_child_index = AtomicUsize::from(state_index);

    // Count this assembly state as searched.
    states_searched.fetch_add(1, Relaxed);

    // Check if any bounds apply; if so, halt this search branch.
    if bound_exceeded(
        mol,
        fragments,
        state_index,
        best_index.load(Relaxed),
        largest_remove,
        bounds,
    ) {
        return state_index;
    }

    // Search for duplicatable fragment
    matches.par_iter().enumerate().for_each(|(i, (h1, h2))| {
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let largest_remove = h1.len();

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            return;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));
        }

        fractures.retain(|i| i.len() > 1);
        fractures.push(h1.clone());

        let output = recurse_index_search_parallel(
            mol,
            &matches[i + 1..],
            &fractures,
            state_index - h1.len() + 1,
            best_index.load(Relaxed).into(),
            largest_remove,
            bounds,
            states_searched.clone(),
        );
        best_child_index.fetch_min(output, Relaxed);

        best_index.fetch_min(best_child_index.load(Relaxed), Relaxed);
    });

    best_child_index.load(Relaxed)
}

/// Computes a molecule's assembly index and related information using a top-
/// down recursive search, parameterized by the specified options.
///
/// See [`EnumerateMode`], [`CanonizeMode`], [`ParallelMode`], [`KernelMode`],
/// and [`Bound`] for details on how to customize the top-down algorithm.
///
/// Notably, bounds are applied in the order they appear in the `bounds` slice.
/// It is generally better to provide bounds that are quick to compute first.
///
/// The results returned are:
/// - `u32`: The molecule's assembly index.
/// - `u32`: The molecule's count of non-overlapping isomorphic subgraph pairs.
/// - `usize`: The total number of assembly states searched, where an assembly
/// state is a collection of fragments; note that some states may be searched
/// and thus counted by this value multiple times.
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
///     CanonizeMode::Nauty,
///     ParallelMode::None,
///     KernelMode::None,
///     &[],
///     false);
///
/// // Compute the molecule's assembly index with parallelism and some bounds.
/// let (fast_index, _, _) = index_search(
///     &anthracene,
///     EnumerateMode::GrowErode,
///     CanonizeMode::Nauty,
///     ParallelMode::Always,
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
    if parallel_mode == ParallelMode::DepthOne {
        panic!("The chosen --parallel mode is not implemented yet!")
    }
    if kernel_mode != KernelMode::None {
        panic!("The chosen --kernel mode is not implemented yet!")
    }
    if memoize {
        panic!("--memoize is not implemented yet!")
    }

    // Enumerate and sort array of non-overlapping, isomorphic subgraph pairs.
    let mut matches = matches(mol, enumerate_mode, canonize_mode).collect::<Vec<_>>();
    matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

    // Initialize fragments as all individual edges.
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    let edge_count = mol.graph().edge_count();

    // Use parallelism if specified by the parallel_mode and if the molecule's
    // number of non-overlapping, isomorphic subgraph pairs is large enough;
    // otherwise, run serially.
    let (index, states_searched) =
        if matches.len() > PARALLEL_MATCH_SIZE_THRESHOLD && parallel_mode == ParallelMode::Always {
            let states_searched = Arc::new(AtomicUsize::from(0));
            let index = recurse_index_search_parallel(
                mol,
                &matches,
                &[init],
                edge_count - 1,
                (edge_count - 1).into(),
                edge_count,
                bounds,
                states_searched.clone(),
            );
            let states_searched = states_searched.load(Relaxed);
            (index as u32, states_searched)
        } else {
            let mut states_searched = 0;
            let index = recurse_index_search_serial(
                mol,
                &matches,
                &[init],
                edge_count - 1,
                edge_count - 1,
                edge_count,
                bounds,
                &mut states_searched,
            );
            (index as u32, states_searched)
        };

    (index, matches.len() as u32, states_searched)
}

/// Computes a molecule's assembly index using an efficient default strategy.
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
        CanonizeMode::Nauty,
        ParallelMode::Always,
        KernelMode::None,
        &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
        false,
    )
    .0
}
