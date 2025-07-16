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

use std::sync::{
    atomic::{AtomicUsize, Ordering::Relaxed},
    Arc,
};

use bit_set::BitSet;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    bounds::{
        Bound,
        log_bound,
        int_bound,
        vec_simple_bound,
        vec_small_frags_bound
    },
    molecule::Molecule,
    utils::connected_components_under_edges,
};

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

#[allow(clippy::too_many_arguments)]
fn recurse_index_search(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    ix: usize,
    largest_remove: usize,
    mut best: usize,
    bounds: &[Bound],
    states_searched: &mut usize,
) -> usize {
    let mut cx = ix;

    *states_searched += 1;

    // Branch and Bound
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::Int => ix - int_bound(fragments, largest_remove) >= best,
            Bound::VecSimple => ix - vec_simple_bound(fragments, largest_remove, mol) >= best,
            Bound::VecSmallFrags => {
                ix - vec_small_frags_bound(fragments, largest_remove, mol) >= best
            }
            // TODO: Remove after all bounds are implemented.
            _ => false,
        };
        if exceeds {
            return ix;
        }
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

        cx = cx.min(recurse_index_search(
            mol,
            &matches[i + 1..],
            &fractures,
            ix - h1.len() + 1,
            largest_remove,
            best,
            bounds,
            states_searched,
        ));
        best = best.min(cx);
    }

    cx
}

#[allow(clippy::too_many_arguments)]
fn parallel_recurse_index_search(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    ix: usize,
    largest_remove: usize,
    best: AtomicUsize,
    bounds: &[Bound],
    states_searched: Arc<AtomicUsize>,
) -> usize {
    let cx = AtomicUsize::from(ix);

    states_searched.fetch_add(1, Relaxed);

    // Branch and Bound
    for bound_type in bounds {
        let best = best.load(Relaxed);
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::Int => ix - int_bound(fragments, largest_remove) >= best,
            Bound::VecSimple => ix - vec_simple_bound(fragments, largest_remove, mol) >= best,
            Bound::VecSmallFrags => {
                ix - vec_small_frags_bound(fragments, largest_remove, mol) >= best
            }
            // TODO: Remove after all bounds are implemented.
            _ => false,
        };
        if exceeds {
            return ix;
        }
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

        let output = parallel_recurse_index_search(
            mol,
            &matches[i + 1..],
            &fractures,
            ix - h1.len() + 1,
            largest_remove,
            best.load(Relaxed).into(),
            bounds,
            states_searched.clone(),
        );
        cx.fetch_min(output, Relaxed);

        best.fetch_min(cx.load(Relaxed), Relaxed);
    });

    cx.load(Relaxed)
}

/// Computes a molecule's assembly index and related information using the
/// specified strategy.
///
/// The first result in the returned tuple is the assembly index of the molecule. The second result
/// gives the number of duplicatable subgraphs (pairs of disjoint and isomorphic subgraphs) in the
/// molecule. The third result is the number of states searched where a new state is considered to
/// be searched each time a duplicatable subgraph is removed.
///
/// If the search space of the molecule is large (>100) parallelization will be used.
///
/// Bounds will be used in the order provided in the `bounds` slice. Execution along a search path
/// will halt immediately after finding a bound that exceeds the current best assembly pathway. It
/// is generally better to provide bounds that are quick to compute first.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::assembly::index_search;
/// use assembly_theory::bounds::Bound;
/// use assembly_theory::loader::parse_molfile_str;
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly index with no bounds.
/// let (slow_index, _, _) = index_search(&anthracene, &[]);
///
/// // Compute the molecule's assembly index with the log and integer bounds.
/// let (fast_index, _, _) =
///     index_search(&anthracene, &[Bound::Log, Bound::Int]);
///
/// assert_eq!(slow_index, 6);
/// assert_eq!(fast_index, 6);
/// # Ok(())
/// # }
/// ```
pub fn index_search(mol: &Molecule, bounds: &[Bound]) -> (u32, u32, usize) {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Create and sort matches array
    let mut matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

    let edge_count = mol.graph().edge_count();

    let (index, total_search) = if matches.len() > PARALLEL_MATCH_SIZE_THRESHOLD {
        let total_search = Arc::new(AtomicUsize::from(0));
        let index = parallel_recurse_index_search(
            mol,
            &matches,
            &[init],
            edge_count - 1,
            edge_count,
            (edge_count - 1).into(),
            bounds,
            total_search.clone(),
        );
        let total_search = total_search.load(Relaxed);
        (index as u32, total_search)
    } else {
        let mut total_search = 0;
        let index = recurse_index_search(
            mol,
            &matches,
            &[init],
            edge_count - 1,
            edge_count,
            edge_count - 1,
            bounds,
            &mut total_search,
        );
        (index as u32, total_search)
    };

    (index, matches.len() as u32, total_search)
}

/// Like [`index_search`], but no parallelism is used.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::assembly::serial_index_search;
/// use assembly_theory::bounds::Bound;
/// use assembly_theory::loader::parse_molfile_str;
/// # fn main() -> Result<(), std::io::Error> {
/// // Load a molecule from a .mol file.
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
///
/// // Compute the molecule's assembly index with no bounds.
/// let (slow_index, _, _) = serial_index_search(&anthracene, &[]);
///
/// // Compute the molecule's assembly index with the log and integer bounds.
/// let (fast_index, _, _) = 
///     serial_index_search(&anthracene, &[Bound::Log, Bound::Int]);
///
/// assert_eq!(slow_index, 6);
/// assert_eq!(fast_index, 6);
/// # Ok(())
/// # }
/// ```
pub fn serial_index_search(mol: &Molecule, bounds: &[Bound]) -> (u32, u32, usize) {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Create and sort matches array
    let mut matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

    let edge_count = mol.graph().edge_count();
    let mut total_search = 0;
    let index = recurse_index_search(
        mol,
        &matches,
        &[init],
        edge_count - 1,
        edge_count,
        edge_count - 1,
        bounds,
        &mut total_search,
    );
    (index as u32, matches.len() as u32, total_search)
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
    index_search(mol, &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags]).0
}
