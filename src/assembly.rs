//! Compute assembly indices of molecules.
//! # Example
//! ```
//! # use std::fs;
//! # use std::path::PathBuf;
//! # use assembly_theory::*;
//! # fn main() -> Result<(), std::io::Error> {
//! # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
//! // Read a molecule data file
//! let molfile = fs::read_to_string(path)?;
//! let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
//!
//! // Compute assembly index of benzene
//! assert_eq!(assembly::index(&benzene), 3);
//! # Ok(())
//! # }
//! ```
use std::{
    collections::BTreeSet,
    sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    },
};

use bit_set::BitSet;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    molecule::Bond, molecule::Element, molecule::Molecule, utils::connected_components_under_edges,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeType {
    bond: Bond,
    ends: (Element, Element),
}

static PARALLEL_MATCH_SIZE_THRESHOLD: usize = 100;

/// Enum to represent the different bounds available during the computation of molecular assembly
/// indices.
/// Bounds are used by `index_search()` to speed up assembly index computations.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Bound {
    /// `Log` bounds by the logarithm base 2 of remaining edges
    Log,
    /// `IntChain` bounds by the length of the smallest addition chain to create the remaining
    /// fragments
    IntChain,
    /// 'VecChainSimple' bounds using addition chain length with the information of the edge types
    /// in a molecule
    VecChainSimple,
    /// 'VecChainSmallFrags' bounds using information on the number of fragments of size 2 in the
    /// molecule
    VecChainSmallFrags,
}

pub fn naive_assembly_depth(mol: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in mol.partitions().unwrap() {
        let l = if left.is_basic_unit() {
            0
        } else {
            naive_assembly_depth(&left)
        };

        let r = if right.is_basic_unit() {
            0
        } else {
            naive_assembly_depth(&right)
        };

        ix = ix.min(l.max(r) + 1)
    }
    ix
}

fn recurse_naive_index_search(
    mol: &Molecule,
    matches: &BTreeSet<(BitSet, BitSet)>,
    fragments: &[BitSet],
    ix: usize,
) -> usize {
    let mut cx = ix;
    for (h1, h2) in matches {
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

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
            fractures.push(h1.clone());
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

            fractures.push(h1.clone());
        }
        cx = cx.min(recurse_naive_index_search(
            mol,
            matches,
            &fractures,
            ix - h1.len() + 1,
        ));
    }
    cx
}

/// Calculates the assembly index of a molecule without using any bounding strategy or
/// parallelization. This function is very inefficient and should only be used as a performance
/// benchmark against other strategies.
pub fn naive_index_search(mol: &Molecule) -> u32 {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    recurse_naive_index_search(
        mol,
        &mol.matches().collect(),
        &[init],
        mol.graph().edge_count() - 1,
    ) as u32
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
            Bound::IntChain => ix - addition_bound(fragments, largest_remove) >= best,
            Bound::VecChainSimple => ix - vec_bound_simple(fragments, largest_remove, mol) >= best,
            Bound::VecChainSmallFrags => {
                ix - vec_bound_small_frags(fragments, largest_remove, mol) >= best
            }
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
            Bound::IntChain => ix - addition_bound(fragments, largest_remove) >= best,
            Bound::VecChainSimple => ix - vec_bound_simple(fragments, largest_remove, mol) >= best,
            Bound::VecChainSmallFrags => {
                ix - vec_bound_small_frags(fragments, largest_remove, mol) >= best
            }
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

/// Computes information related to the assembly index of a molecule using the provided bounds.
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
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// use assembly_theory::assembly::{Bound, index_search};
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene naively, with no bounds.
/// let (slow_index, _, _) = index_search(&benzene, &[]);
///
/// // Compute assembly index of benzene with the log and integer chain bounds
/// let (fast_index, _, _) = index_search(&benzene, &[Bound::Log, Bound::IntChain]);
///
/// assert_eq!(slow_index, 3);
/// assert_eq!(fast_index, 3);
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
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// use assembly_theory::assembly::{Bound, serial_index_search};
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene naively, with no bounds.
/// let (slow_index, _, _) = serial_index_search(&benzene, &[]);
///
/// // Compute assembly index of benzene with the log and integer chain bounds
/// let (fast_index, _, _) = serial_index_search(&benzene, &[Bound::Log, Bound::IntChain]);
///
/// assert_eq!(slow_index, 3);
/// assert_eq!(fast_index, 3);
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

fn log_bound(fragments: &[BitSet]) -> usize {
    let mut size = 0;
    for f in fragments {
        size += f.len();
    }

    size - (size as f32).log2().ceil() as usize
}

fn addition_bound(fragments: &[BitSet], m: usize) -> usize {
    let mut max_s: usize = 0;
    let mut frag_sizes: Vec<usize> = Vec::new();

    for f in fragments {
        frag_sizes.push(f.len());
    }

    let size_sum: usize = frag_sizes.iter().sum();

    // Test for all sizes m of largest removed duplicate
    for max in 2..m + 1 {
        let log = (max as f32).log2().ceil();
        let mut aux_sum: usize = 0;

        for len in &frag_sizes {
            aux_sum += (len / max) + (len % max != 0) as usize
        }

        max_s = max_s.max(size_sum - log as usize - aux_sum);
    }

    max_s
}

// Count number of unique edges in a fragment
// Helper function for vector bounds
fn unique_edges(fragment: &BitSet, mol: &Molecule) -> Vec<EdgeType> {
    let g = mol.graph();
    let mut nodes: Vec<Element> = Vec::new();
    for v in g.node_weights() {
        nodes.push(v.element());
    }
    let edges: Vec<petgraph::prelude::EdgeIndex> = g.edge_indices().collect();
    let weights: Vec<Bond> = g.edge_weights().copied().collect();

    // types will hold an element for every unique edge type in fragment
    let mut types: Vec<EdgeType> = Vec::new();
    for idx in fragment.iter() {
        let bond = weights[idx];
        let e = edges[idx];

        let (e1, e2) = g.edge_endpoints(e).expect("bad");
        let e1 = nodes[e1.index()];
        let e2 = nodes[e2.index()];
        let ends = if e1 < e2 { (e1, e2) } else { (e2, e1) };

        let edge_type = EdgeType { bond, ends };

        if types.iter().any(|&t| t == edge_type) {
            continue;
        } else {
            types.push(edge_type);
        }
    }

    types
}

fn vec_bound_simple(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
    // Calculate s (total number of edges)
    // Calculate z (number of unique edges)
    let mut s = 0;
    for f in fragments {
        s += f.len();
    }

    let mut union_set = BitSet::new();
    for f in fragments {
        union_set.union_with(f);
    }
    let z = unique_edges(&union_set, mol).len();

    (s - z) - ((s - z) as f32 / m as f32).ceil() as usize
}

fn vec_bound_small_frags(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
    let mut size_two_fragments: Vec<BitSet> = Vec::new();
    let mut large_fragments: Vec<BitSet> = fragments.to_owned();
    let mut indices_to_remove: Vec<usize> = Vec::new();

    // Find and remove fragments of size 2
    for (i, frag) in fragments.iter().enumerate() {
        if frag.len() == 2 {
            indices_to_remove.push(i);
        }
    }
    for &index in indices_to_remove.iter().rev() {
        let removed_bitset = large_fragments.remove(index);
        size_two_fragments.push(removed_bitset);
    }

    // Compute z = num unique edges of large_fragments NOT also in size_two_fragments
    let mut fragments_union = BitSet::new();
    let mut size_two_fragments_union = BitSet::new();
    for f in fragments {
        fragments_union.union_with(f);
    }
    for f in size_two_fragments.iter() {
        size_two_fragments_union.union_with(f);
    }
    let z = unique_edges(&fragments_union, mol).len()
        - unique_edges(&size_two_fragments_union, mol).len();

    // Compute s = total number of edges in fragments
    // Compute sl = total number of edges in large fragments
    let mut s = 0;
    let mut sl = 0;
    for f in fragments {
        s += f.len();
    }
    for f in large_fragments {
        sl += f.len();
    }

    // Find number of unique size two fragments
    let mut size_two_types: Vec<(EdgeType, EdgeType)> = Vec::new();
    for f in size_two_fragments.iter() {
        let mut types = unique_edges(f, mol);
        types.sort();
        if types.len() == 1 {
            size_two_types.push((types[0], types[0]));
        } else {
            size_two_types.push((types[0], types[1]));
        }
    }
    size_two_types.sort();
    size_two_types.dedup();

    s - (z + size_two_types.len() + size_two_fragments.len())
        - ((sl - z) as f32 / m as f32).ceil() as usize
}

/// Computes the assembly index of a molecule using an effecient bounding strategy
/// # Example
/// ```
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene
/// assert_eq!(assembly::index(&benzene), 3);
/// # Ok(())
/// # }
/// ```
pub fn index(m: &Molecule) -> u32 {
    index_search(
        m,
        &[
            Bound::IntChain,
            Bound::VecChainSimple,
            Bound::VecChainSmallFrags,
        ],
    )
    .0
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs, path::PathBuf};

    use csv::ReaderBuilder;

    use crate::loader;

    use super::*;

    // Read Master CSV
    fn read_dataset_index(dataset: &str) -> HashMap<String, u32> {
        let path = format!("./data/{dataset}/ma-index.csv");
        let mut reader = ReaderBuilder::new()
            .from_path(path)
            .expect("ma-index.csv does not exist.");
        let mut index_records = HashMap::new();
        for result in reader.records() {
            let record = result.expect("ma-index.csv is malformed.");
            let record = record.iter().collect::<Vec<_>>();
            index_records.insert(
                record[0].to_string(),
                record[1]
                    .to_string()
                    .parse::<u32>()
                    .expect("Assembly index is not an integer."),
            );
        }
        index_records
    }

    // Read Test CSV
    fn test_molecule<F>(function: F, dataset: &str, filename: &str)
    where
        F: Fn(&Molecule) -> u32,
    {
        let path = PathBuf::from(format!("./data/{dataset}/{filename}"));
        let molfile = fs::read_to_string(path).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        let dataset = read_dataset_index(dataset);
        let ground_truth = dataset
            .get(filename)
            .expect("Index dataset has no ground truth value");
        let index = function(&molecule);
        assert_eq!(index, *ground_truth);
    }

    #[test]
    fn all_bounds_benzene() {
        test_molecule(index, "checks", "benzene.mol");
    }

    #[test]
    fn all_bounds_aspirin() {
        test_molecule(index, "checks", "aspirin.mol");
    }

    #[test]
    #[ignore = "expensive test"]
    fn all_bounds_morphine() {
        test_molecule(index, "checks", "morphine.mol");
    }

    #[test]
    fn naive_method_benzene() {
        test_molecule(naive_index_search, "checks", "benzene.mol");
    }

    #[test]
    fn naive_method_aspirin() {
        test_molecule(naive_index_search, "checks", "aspirin.mol");
    }

    #[test]
    #[ignore = "expensive test"]
    fn naive_method_morphine() {
        test_molecule(naive_index_search, "checks", "morphine.mol");
    }
}
