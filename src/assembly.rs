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
    sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    }, time::{Duration, Instant}
};

use bit_set::BitSet;
use clap::ValueEnum;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator, IndexedParallelIterator};

use crate::{
    bounds::Bound, canonize::CanonizeMode, enumerate::EnumerateMode, kernels::KernelMode, matches::Matches, memoize::{MemoizeMode, NewCache}, molecule:: Molecule, state::State, utils::connected_components_under_edges
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
/// the given pair of non-overlapping, isomorphic subgraphs and then adding one
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

fn _is_usable(state: &[BitSet], h1: &BitSet, h2: &BitSet, masks: &mut Vec<Vec<BitSet>>) -> bool {
    let f1 = state.iter().enumerate().find(|(_, c)| h1.is_subset(c));
    let f2 = state.iter().enumerate().find(|(_, c)| h2.is_subset(c));
    
    if let (Some((i1, _)), Some((i2, _))) = (f1, f2) {
        masks[h1.len() - 2][i1].union_with(h1);
        masks[h1.len() - 2][i2].union_with(h2);
        true
    } else {
        false
    }
}

/// Recursive helper for [`index_search`], only public for benchmarking.
///
/// Inputs:
/// - `mol`: The molecule whose assembly index is being calculated.
/// - `matches`: The remaining non-overlapping isomorphic subgraph pairs.
/// - `graph`: If clique is enabled, the graph of compatible isomorphic subgraph paris.
/// - `subgraph`: A bitset with length of matches. Has a 1 for every match to be searched in this state.
/// - `removal_order`: TODO
/// - `state`: The current assembly state, i.e., a list of fragments.
/// - `state_index`: This assembly state's upper bound on the assembly index,
///   i.e., edges(mol) - 1 - [edges(subgraphs removed) - #(subgraphs removed)].
/// - `best_index`: The smallest assembly index for all assembly states so far.
/// - `bounds`: The list of bounding strategies to apply.
/// - `cache`: TODO
/// - `parallel_mode`: The parallelism mode for this state's match iteration.
/// - `kernel_mode`: The kernelization mode for this state.
///
/// Returns, from this assembly state and any of its descendents:
/// - `usize`: An updated upper bound on the assembly index. (Note: If this
///   state is pruned by bounds or deemed redundant by memoization, then the
///   upper bound returned is unchanged.)
/// - `usize`: The number of assembly states searched.
#[allow(clippy::too_many_arguments)]
pub fn recurse_index_search(
    mol: &Molecule,
    matches: &Matches,
    state: &State,
    best_index: Arc<AtomicUsize>,
    bounds: &[Bound],
    cache: &mut NewCache,
    parallel_mode: ParallelMode,
) -> (usize, usize) {
    // if removal_order.len() == 1 {println!("{:?}", removal_order);}
    // println!("{:?}", state);

    // Memoization
    if cache.memoize_state(mol, state) {
        return (state.index(), 1);
    }

    // Generate matches
    let (frags, valid_matches) = matches.generate_matches(mol, state, best_index.load(Relaxed), bounds);

    // Keep track of the best assembly index found in any of this assembly
    // state's children and the number of states searched, including this one.
    let best_child_index = AtomicUsize::from(state.index());
    let states_searched = AtomicUsize::from(1);
    
    // Define a closure that handles recursing to a new assembly state based on
    // the given (enumerated) pair of non-overlapping isomorphic subgraphs.
    let recurse_on_match = |i: usize, v: (usize, usize)| {
        //let (h1, h2) = &matches[v];
        let h1 = matches.get_frag(v.0);
        let h2 = matches.get_frag(v.1);
        
        // TODO: try keeping track of fragment to elimate some some search
        if let Some(fragments) = fragments(mol, &frags, h1, h2) {
            // If using depth-one parallelism, all descendant states should be
            // computed serially.
            let new_parallel = if parallel_mode == ParallelMode::DepthOne {
                ParallelMode::None
            } else {
                parallel_mode
            };

            let new_state = state.update(fragments, h1.len(), i, matches.match_id(&v).unwrap());

            // Recurse using the remaining matches and updated fragments.
            let (child_index, child_states_searched) = recurse_index_search(
                mol,
                matches,
                &new_state,
                best_index.clone(),
                bounds,
                &mut cache.clone(),
                new_parallel,
            );

            // Update the best assembly indices (across children states and
            // the entire search) and the number of descendant states searched.
            best_child_index.fetch_min(child_index, Relaxed);
            best_index.fetch_min(best_child_index.load(Relaxed), Relaxed);
            states_searched.fetch_add(child_states_searched, Relaxed);
        }
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    if parallel_mode == ParallelMode::None {
        valid_matches
            .iter()
            .enumerate()
            .for_each(|(i, v)| recurse_on_match(i, *v));
    } else {
        valid_matches
            .par_iter()
            .enumerate()
            .for_each(|(i, v)| recurse_on_match(i, *v));
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
/// let (slow_index, _, _) = index_search(
///     &anthracene,
///     EnumerateMode::GrowErode,
///     CanonizeMode::TreeNauty,
///     ParallelMode::None,
///     MemoizeMode::None,
///     KernelMode::None,
///     &[],
/// );
///
/// // Compute the molecule's assembly index with parallelism, memoization, and
/// // some bounds.
/// let (fast_index, _, _) = index_search(
///     &anthracene,
///     EnumerateMode::GrowErode,
///     CanonizeMode::TreeNauty,
///     ParallelMode::DepthOne,
///     MemoizeMode::CanonIndex,
///     KernelMode::None,
///     &[Bound::Log, Bound::Int],
/// );
///
/// assert_eq!(slow_index, 6);
/// assert_eq!(fast_index, 6);
/// # Ok(())
/// # }
/// ```
pub fn index_search(
    mol: &Molecule,
    _enumerate_mode: EnumerateMode,
    canonize_mode: CanonizeMode,
    parallel_mode: ParallelMode,
    memoize_mode: MemoizeMode,
    _kernel_mode: KernelMode,
    bounds: &[Bound],
    _clique: bool,
) -> (u32, u32, usize) {
    // Enumerate non-overlapping isomorphic subgraph pairs.
    let matches = Matches::new(mol, canonize_mode);

    // Create memoization cache.
    let mut cache = NewCache::new(memoize_mode);

    // Initialize state
    let state = State::new(mol);

    // Search for the shortest assembly pathway recursively.
    let best_index = Arc::new(AtomicUsize::from(mol.graph().edge_count() - 1));
    let (index, states_searched) = recurse_index_search(
        mol,
        &matches,
        &state,
        best_index,
        bounds,
        &mut cache,
        parallel_mode,
    );

    // println!("Bounded: {}", cache.count());

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
        MemoizeMode::CanonIndex,
        KernelMode::None,
        &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
        false,
    )
    .0
}
