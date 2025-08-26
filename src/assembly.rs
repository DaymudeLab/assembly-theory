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
    collections::{HashMap, HashSet}, sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    }
};

use bit_set::BitSet;
use clap::ValueEnum;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use petgraph::graph::EdgeIndex;

use crate::{
    bounds::{bound_exceeded, Bound},
    canonize::{canonize, CanonizeMode, Labeling},
    enumerate::{EnumerateMode},
    kernels::{KernelMode, deletion_kernel, inclusion_kernel},
    memoize::{Cache, MemoizeMode, NewCache},
    molecule::Molecule,
    utils::{connected_components_under_edges, edge_neighbors},
    reductions::CompatGraph,
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
/// Return all pairs of non-overlapping, isomorphic subgraphs in the molecule,
/// sorted to guarantee deterministic iteration.
pub fn matches(mol: &Molecule) -> Vec<(BitSet, BitSet)> {
    let num_edges = mol.graph().edge_count();
    let mut subgraphs: HashSet<BitSet> = HashSet::new();
    let mut size_one_subgraphs: HashSet<BitSet> = HashSet::new();
    let mut matches: Vec<(BitSet, BitSet)> = Vec::new();
    let mut buckets: HashMap<Labeling, Vec<BitSet>> = HashMap::new();

    // Generate subgraphs with one edge
    for i in 0..num_edges {
        let mut new_bitset = BitSet::with_capacity(num_edges);
        new_bitset.insert(i);
        size_one_subgraphs.insert(new_bitset);
    }

    // Generate subgraphs with two edges
    for subgraph in size_one_subgraphs.iter() {
        let mut neighborhood = BitSet::with_capacity(num_edges);
        for e in subgraph {
            neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
        }
        neighborhood.difference_with(subgraph);

        for e in neighborhood.iter() {
            let mut new_subgraph = subgraph.clone();
            new_subgraph.insert(e);
            subgraphs.insert(new_subgraph);
        }
    }

    while subgraphs.len() > 0 {
        buckets.clear();
        // Generate buckets
        for g in subgraphs.iter() {
            let canon = canonize(mol, g, CanonizeMode::TreeNauty);
            buckets.entry(canon)
                .and_modify(|v| v.push(g.clone()))
                .or_insert(vec![g.clone()]);
        }

        subgraphs.clear();
        for (_, bucket) in buckets.iter() {
            // Generate matches
            let mut has_match = BitSet::with_capacity(bucket.len());
            for (i, first) in bucket.iter().enumerate() {
                for (j, second) in bucket[i+1..].iter().enumerate() {
                    if first.intersection(second).count() == 0 {
                        if first > second {
                        matches.push((first.clone(), second.clone()));
                        } else {
                            matches.push((second.clone(), first.clone()));
                        }
                        has_match.insert(i);
                        has_match.insert(i + 1 + j);
                    }
                }
            }

            // Generate new subgraphs
            // Only extend if subgraph has a match
            for i in has_match.iter() {
                let subgraph = &bucket[i];
                let mut neighborhood = BitSet::with_capacity(num_edges);
                for e in subgraph {
                    neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                }
                neighborhood.difference_with(subgraph);

                for e in neighborhood.iter() {
                    let mut new_subgraph = subgraph.clone();
                    new_subgraph.insert(e);
                    subgraphs.insert(new_subgraph);
                }
            }
        }
    }

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

    matches
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

fn is_usable(state: &[BitSet], h1: &BitSet, h2: &BitSet, masks: &mut Vec<Vec<BitSet>>) -> bool {
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
    matches: &Vec<(BitSet, BitSet)>,
    graph: &Option<CompatGraph>,
    mut subgraph: BitSet,
    removal_order: Vec<usize>,
    state: &[BitSet],
    state_index: usize,
    best_index: Arc<AtomicUsize>,
    bounds: &[Bound],
    cache: &mut NewCache,
    parallel_mode: ParallelMode,
    kernel_mode: KernelMode,
) -> (usize, usize) {
    //println!("{:?}", removal_order);
    //println!("{:?}", state);

    // An upper bound on the size of fragments that can be
    // removed from this or any descendant assembly state.
    let largest_remove = {
        if let Some(v) = subgraph.iter().next() {
            matches[v].0.len()
        }
        else {
            return (state_index, 1);
        }
    };

    let mut masks: Vec<Vec<BitSet>> = vec![vec![BitSet::with_capacity(mol.graph().edge_count()); state.len()]; largest_remove - 1];
    for v in subgraph.iter() {
        let m = &matches[v];
        is_usable(state, &m.0, &m.1, &mut masks);
    }
    /*while masks.iter().last().unwrap().len() == 0 {
        masks.pop();
    }*/
    /*let mut state = Vec::new();
    for (i, m) in masks[0].iter().enumerate() {
        state[i] = m.clone();
    }*/
    /*let mut temp = Vec::new();
    let state = {
        if removal_order.len() == 0 {
            for m in masks[0].iter() {
                if m.len() != 0 {
                    temp.extend(connected_components_under_edges(mol.graph(), m));
                }
            }
            temp.retain(|i| i.len() > 1);
            &temp
        }
        else {
            state
        }
    };*/

    /*let mut masks: Vec<Vec<BitSet>> = vec![vec![BitSet::with_capacity(mol.graph().edge_count()); state.len()]; largest_remove - 1];
    for v in subgraph.iter() {
        let m = &matches[v];
        is_usable(&state, &m.0, &m.1, &mut masks);
    }*/

    // If any bounds would prune this assembly state or if memoization is
    // enabled and this assembly state is preempted by the cached state, halt.
    if bound_exceeded(
        mol,
        matches,
        graph,
        &state,
        &subgraph,
        state_index,
        best_index.load(Relaxed),
        largest_remove,
        bounds,
        &masks,
    ) || cache.memoize_state(mol, &state, state_index)
    {
        return (state_index, 1);
    }

    let mut state = Vec::new();
    for (i, m) in masks[0].iter().enumerate() {
        state.push(m.clone());
    }

    // Apply kernels
    let mut must_include = matches.len();
    if kernel_mode != KernelMode::None {
        let g = graph.as_ref().unwrap();

        // Deletion kernel
        // Returns a subgraph without nodes that never occur in an optimal solution.
        subgraph = deletion_kernel(matches, g, subgraph);
        
        // Inclusion kernel.
        // must_include is the first match in the matches list that will be included in an optimal solution.
        must_include = inclusion_kernel(matches, g, &subgraph);
    }

    // Keep track of the best assembly index found in any of this assembly
    // state's children and the number of states searched, including this one.
    let best_child_index = AtomicUsize::from(state_index);
    let states_searched = AtomicUsize::from(1);

    // Define a closure that handles recursing to a new assembly state based on
    // the given (enumerated) pair of non-overlapping isomorphic subgraphs.
    let recurse_on_match = |v: usize| {
        let (h1, h2) = &matches[v];
        if let Some(fragments) = fragments(mol, &state, h1, h2) {
            // If using depth-one parallelism, all descendant states should be
            // computed serially.
            let new_parallel = if parallel_mode == ParallelMode::DepthOne {
                ParallelMode::None
            } else {
                parallel_mode
            };

            // If kernelizing once, do not kernelize again. If kernelizing at depth-one,
            // kernelize once more at the beginning of each descendent.
            let new_kernel = if kernel_mode == KernelMode::Once {
                KernelMode::None
            } else if kernel_mode == KernelMode::DepthOne {
                KernelMode::Once
            } else {
                kernel_mode
            };

            // Update subgraph
            // If using CompatGraph, only consider the neighborhood of the node removed.
            // Otherwise, only consider the matches that occur after the removed match in the list.
            let mut sub_clone;
            if let Some(g) = graph {
                sub_clone = subgraph.clone();
                sub_clone.intersect_with(&g.forward_neighbors(v, &subgraph));
            }
            else {
                sub_clone = BitSet::with_capacity(matches.len());
                for j in (v+1)..matches.len() {
                    sub_clone.insert(j);
                }
            }

            // Recurse using the remaining matches and updated fragments.
            let (child_index, child_states_searched) = recurse_index_search(
                mol,
                matches,
                graph,
                sub_clone,
                {
                    let mut clone = removal_order.clone();
                    clone.push(v);
                    clone
                },
                &fragments,
                state_index - h1.len() + 1,
                best_index.clone(),
                bounds,
                &mut cache.clone(),
                new_parallel,
                new_kernel,
            );

            // Update the best assembly indices (across children states and
            // the entire search) and the number of descendant states searched.
            best_child_index.fetch_min(child_index, Relaxed);
            best_index.fetch_min(best_child_index.load(Relaxed), Relaxed);
            states_searched.fetch_add(child_states_searched, Relaxed);
        }
    };

    // Use the iterator type corresponding to the specified parallelism mode.
    // Only search on nodes with v <= must_include since any searches after that would
    // not use must_include.
    if parallel_mode == ParallelMode::None {
        subgraph
            .iter()
            .filter(|v| *v <= must_include)
            .for_each(|v| recurse_on_match(v));
    } else {
        let _ = subgraph
            .iter()
            .filter(|v| *v <= must_include)
            .collect::<Vec<usize>>()
            .par_iter()
            .for_each(|v| recurse_on_match(*v));
    }

    if removal_order.len() == 1 {
        println!("{:?}, {}", removal_order, states_searched.load(Relaxed));
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
    kernel_mode: KernelMode,
    bounds: &[Bound],
    clique: bool,
) -> (u32, u32, usize) {
    // Catch not-yet-implemented modes.

    // Enumerate non-overlapping isomorphic subgraph pairs.
    let matches = matches(mol);
    /*for (i, m) in matches.iter().enumerate() {
        println!("{i}: {:?}", m);
    }*/

    // Create memoization cache.
    //let mut cache = Cache::new(memoize_mode, canonize_mode);
    let mut cache = NewCache::new();

    let graph = {
        if clique {
            Some(CompatGraph::new(&matches))
        }
        else {
            None
        }
    };

    // Initialize the first fragment as the entire graph.
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    let mut subgraph = BitSet::with_capacity(matches.len());
    for i in 0..matches.len() {
        subgraph.insert(i);
    }

    // Search for the shortest assembly pathway recursively.
    let edge_count = mol.graph().edge_count();
    let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
    //let best_index = Arc::new(AtomicUsize::from(22));
    let (index, states_searched) = recurse_index_search(
        mol,
        &matches,
        &graph,
        subgraph,
        Vec::new(),
        &[init],
        edge_count - 1,
        best_index,
        bounds,
        &mut cache,
        parallel_mode,
        kernel_mode,
    );

    println!("Bounded: {}", cache.count());

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
