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
    collections::{BTreeMap, HashMap, HashSet}, sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    }
};

use bit_set::BitSet;
use clap::ValueEnum;
use graph_canon::canon;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use petgraph::graph::EdgeIndex;
use itertools::Itertools;
use std::cmp::{max, min};

use crate::{
    bounds::{bound_exceeded, Bound},
    canonize::{canonize, CanonizeMode, Labeling},
    enumerate::{EnumerateMode},
    kernels::{KernelMode, deletion_kernel, inclusion_kernel},
    memoize::{MemoizeMode, NewCache},
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

#[derive(Debug)]
pub struct DagNode {
    fragment: BitSet,
    canon_id: usize,
    // added_edge: usize,
    children: Vec<usize>,
}

impl DagNode {
    pub fn new(fragment: BitSet, canon_id: usize) -> Self {
        Self {
            fragment,
            canon_id,
            children: Vec::new(),
        }
    }
}


pub fn matches(mol: &Molecule, canonize_mode: CanonizeMode) -> (HashMap<(usize, usize), usize>, Vec<DagNode>) {
    let num_edges = mol.graph().edge_count();

    let mut subgraphs: Vec<usize> = Vec::new();
    let mut constructed_frags: HashSet<BitSet> = HashSet::new();
    let mut frag_to_id: HashMap<BitSet, usize> = HashMap::new();

    let mut dag: Vec<DagNode> = Vec::with_capacity(num_edges);
    let mut matches: Vec<(usize, usize)> = Vec::new();

    let mut next_canon_id = 0;

    // Generate subgraphs with one edge
    for i in 0..num_edges {
        let mut new_bitset = BitSet::with_capacity(num_edges);
        new_bitset.insert(i);
        
        let node = DagNode::new(new_bitset, 0);

        dag.push(node);
        subgraphs.push(i);
    }

    next_canon_id += 1;

    while subgraphs.len() > 0 {
        // Initialize new buckets and subgraph list
        let mut buckets: BTreeMap<Labeling, Vec<(BitSet, usize)>> = BTreeMap::new();
        let mut child_subgraphs: Vec<usize> = Vec::new();

        // Extend and bucket new subgraphs
        for parent_idx in subgraphs.iter() {
            // Get fragment to be extended
            let fragment = &dag[*parent_idx].fragment;

            // Construct edge neighborhood of fragment
            let mut neighborhood = BitSet::with_capacity(num_edges);
            for e in fragment {
                neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
            }
            neighborhood.difference_with(fragment);

            // Extend fragment with edges from neighborhood
            for e in neighborhood.iter() {
                // Create extended fragment
                let mut new_fragment = fragment.clone();
                new_fragment.insert(e);

                // Check if fragment has already been created
                // If yes, continue
                if constructed_frags.contains(&new_fragment) {
                    continue;
                }

                // Bucket subgraph
                // Gives ownership of new_fragment to buckets
                let canon = canonize(mol, &new_fragment, canonize_mode);
                buckets.entry(canon)
                    .and_modify(|bucket| bucket.push((new_fragment.clone(), *parent_idx)))
                    .or_insert(vec![(new_fragment.clone(), *parent_idx)]);
            
                // Add to constructed_frags
                constructed_frags.insert(new_fragment);
            }
        }

        // Go through buckets to create matches
        for fragments in buckets.values() {
            // Variables to track if a match has been found
            let mut frag_has_match = BitSet::with_capacity(fragments.len());

            // Loop through pairs of bitsets to find matches
            for pair in fragments.iter().enumerate().combinations(2) {
                // Extract pair of fragments
                let frag1 = &pair[0].1.0;
                let frag2 = &pair[1].1.0;

                // Check for disjointness
                // If the fragments are disjoint, they form a match
                if frag1.is_disjoint(&frag2) {
                    // Extract indices of fragments
                    let idx1 = pair[0].0;
                    let idx2 = pair[1].0;

                    // If this is the first match for a fragment, create a dag node
                    if !frag_has_match.contains(idx1) {
                        // Create new dag node
                        let new_dag_node = DagNode::new(frag1.clone(), next_canon_id);
                        let child_idx = dag.len();

                        // Add child to dag list
                        dag.push(new_dag_node);

                        // Add index to parent's child list
                        let parent_idx = pair[0].1.1;
                        let mut parent = &mut dag[parent_idx];
                        parent.children.push(child_idx);

                        // Add dag node to next subgraphs list
                        child_subgraphs.push(child_idx);

                        // Add to frag_to_id map
                        frag_to_id.insert(frag1.clone(), child_idx);
                    }
                    if !frag_has_match.contains(idx2) {
                        // Create new dag node
                        let new_dag_node = DagNode::new(frag2.clone(), next_canon_id);
                        let child_idx = dag.len();

                        // Add child to dag list
                        dag.push(new_dag_node);

                        // Add index to parent's child list
                        let parent_idx = pair[1].1.1;
                        let mut parent = &mut dag[parent_idx];
                        parent.children.push(child_idx);

                        // Add dag node to next subgraphs list
                        child_subgraphs.push(child_idx);

                        // Add to frag_to_id map
                        frag_to_id.insert(frag2.clone(), child_idx);
                    }

                    // Get fragment ids and add to matches
                    let frag1_id = frag_to_id.get(&frag1).unwrap();
                    let frag2_id = frag_to_id.get(&frag2).unwrap();

                    if frag1_id > frag2_id {
                        matches.push((*frag1_id, *frag2_id));
                    }
                    else {
                        matches.push((*frag2_id, *frag1_id));
                    }

                    // Mark that these fragments have matches
                    frag_has_match.insert(idx1);
                    frag_has_match.insert(idx2);
                }
            }

            // If there was a match, increment canon_id counter
            if !frag_has_match.is_empty() {
                next_canon_id += 1;
            }
        }

        // Update to the new subgraphs
        subgraphs = child_subgraphs;
    }

    // Sort matches
    // Larger frag_ids get a smaller match_id
    // Thus fragments with many edges have small ids
    matches.sort();
    matches.reverse();

    // give matches ids
    let mut next_match_id = 0;
    let mut match_to_id: HashMap<(usize, usize), usize> = HashMap::new();
    for m in matches {
        match_to_id.insert(m, next_match_id);
        next_match_id += 1;
    }

    (match_to_id, dag)
}


/// Helper function for [`index_search`]; only public for benchmarking.
///
/// Return all pairs of non-overlapping, isomorphic subgraphs in the molecule,
/// sorted to guarantee deterministic iteration.
/*pub fn _matches(mol: &Molecule) -> (Vec<(usize, usize)>, DagNode) {
    let num_edges = mol.graph().edge_count();

    let mut subgraphs: HashMap<BitSet, usize> = HashMap::new();
    let mut size_one_subgraphs: HashMap<BitSet, usize> = HashMap::new();
    let mut buckets: HashMap<Labeling, Vec<DagNode>> = HashMap::new();

    let mut dag: Vec<DagNode> = Vec::with_capacity(num_edges);
    let mut canon_ids: HashMap<Labeling, usize> = HashMap::new();
    let mut matches: HashMap<(usize, usize), usize> = Hashmap::new();

    let mut next_frag_id = 1;
    let mut next_match_id = 0;
    let mut next_canon_id = 1;

    // Generate subgraphs with one edge
    for i in 0..num_edges {
        let mut new_bitset = BitSet::with_capacity(num_edges);
        new_bitset.insert(i);
        size_one_subgraphs.insert(new_bitset, i);

        dag.push(DagNode::new(new_bitset, 0, 0));
    }

    // Generate subgraphs with two edges
    for (subgraph, idx) in size_one_subgraphs.iter() {
        let mut neighborhood = BitSet::with_capacity(num_edges);
        for e in subgraph {
            neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
        }
        neighborhood.difference_with(subgraph);

        for e in neighborhood.iter() {
            let mut new_subgraph = subgraph.clone();
            new_subgraph.insert(e);
            subgraphs.entry(new_subgraph)
                .or_insert(*idx);
        }
    }

    // Generate subgraphs with more than two edges
    while subgraphs.len() > 0 {
        buckets.clear();
        // Generate buckets
        for (g, parent) in subgraphs.iter() {
            let canon = canonize(mol, g, CanonizeMode::TreeNauty);            
            buckets.entry(canon)
                .and_modify(|v| v.push((g.clone(), *parent)))
                .or_insert(vec![(g.clone(), *parent)]);
        }

        subgraphs.clear();
        let mut frag_to_id: HashMap<BitSet, usize> = HashMap::new();
        for (_, bucket) in buckets.iter() {
            // Generate matches
            let mut has_match = BitSet::with_capacity(bucket.len());
            for (i, (first, first_parent)) in bucket.iter().enumerate() {
                for (j, (second, second_parent)) in bucket[i+1..].iter().enumerate() {
                    if first.intersection(second).count() == 0 {
                        if !has_match.contains(i) {
                            frag_to_id.insert(first.clone(), next_id);

                            dag.insert(next_id, (first.clone(), Vec::new()));
                            dag.entry(*first_parent)
                                .and_modify(|(_, children)| children.push(next_id));

                            let mut neighborhood = BitSet::with_capacity(num_edges);
                            for e in first {
                                neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                            }
                            neighborhood.difference_with(first);

                            for e in neighborhood.iter() {
                                let mut new_subgraph = first.clone();
                                new_subgraph.insert(e);
                                subgraphs.insert(new_subgraph, next_id);
                            }

                            next_id += 1;
                        }
                        if !has_match.contains(i + 1 + j) {
                            frag_to_id.insert(second.clone(), next_id);

                            dag.insert(next_id, (second.clone(), Vec::new()));
                            dag.entry(*second_parent)
                                .and_modify(|(_, children)| children.push(next_id));

                            let mut neighborhood = BitSet::with_capacity(num_edges);
                            for e in second {
                                neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                            }
                            neighborhood.difference_with(second);

                            for e in neighborhood.iter() {
                                let mut new_subgraph = second.clone();
                                new_subgraph.insert(e);
                                subgraphs.insert(new_subgraph, next_id);
                            }

                            next_id += 1;
                        }

                        let id1 = frag_to_id.get(first).unwrap();
                        let id2 = frag_to_id.get(second).unwrap();
                        if first > second {
                            matches.push((*id1, *id2));
                        } else {
                            matches.push((*id1, *id2));
                        }

                        has_match.insert(i);
                        has_match.insert(i + 1 + j);
                    }
                }
            }
        }
    }

    // Sort matches
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

    (matches, dag)
}*/

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
    matches: &HashMap<(usize, usize), usize>,
    dag: &Vec<DagNode>,
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
    
    let num_edges = mol.graph().edge_count();
    let mut enum_matches: Vec<(usize, usize)> = Vec::new();
    
    // MATCHES METHOD 1
    // EXTEND ALWAYS AND DO MATCHING AT END
    let mut subgraphs: Vec<usize> = Vec::new();
    let mut buckets: HashMap<usize, Vec<usize>> = HashMap::new();

    // create subgraphs of size 1
    for fragment in state.iter() {
        for edge_idx in fragment.iter() {
            subgraphs.push(edge_idx);
        }
    }

    while subgraphs.len() > 0 {
        // Subgraphs to extend in next loop
        let mut new_subgraphs = Vec::new();

        // Extend each subgraph and bucket
        for frag_id in subgraphs {
            // Get ids of extensions of this subgraph
            let children_ids = &dag[frag_id].children;

            for child_id in children_ids {
                // Check if this extension is valid
                // TODO: save fragment id
                // TODO: save edge id in DAG
                let child_frag = &dag[*child_id].fragment;
                let valid = state.iter().any(|frag| child_frag.is_subset(frag));

                if valid {
                    let canon_id = &dag[*child_id].canon_id;
                    
                    // Add fragment to bucket
                    buckets.entry(*canon_id)
                        .and_modify(|bucket| bucket.push(*child_id))
                        .or_insert(vec![*child_id]);
                    
                    // Add fragment to new subgraphs
                    new_subgraphs.push(*child_id);
                }
            }
        }

        // Update to new subgraphs
        subgraphs = new_subgraphs;
    }

    // Generate matches from buckets of ismorphic fragments
    for fragments in buckets.values() {
        // Loop over pairs of fragments
        for pair in fragments.iter().combinations(2) {
            let frag1 = *max(pair[0], pair[1]);
            let frag2 = *min(pair[0], pair[1]);

            // Check for valid match
            // If valid, add to matches to iterate over
            if let Some(x) = matches.get(&(frag1, frag2)) {
                enum_matches.push((frag1, frag2));
            }
        }
    }

    // If there are no matches, return
    if enum_matches.len() == 0 {
        return (state_index, 1);
    }

    enum_matches.sort();
    enum_matches.reverse();

    // Modify mask
    // TODO: try storing fragment so search can be reduced
    // TODO: try with a set for frag_ids to only do intersection once per frag
    let largest_remove = dag[enum_matches[0].0].fragment.len();
    let mut masks = vec![
        vec![BitSet::with_capacity(num_edges); state.len()];
        largest_remove - 1
    ];

    for m in enum_matches.iter() {
        // Extract fragments
        let frag1 = &dag[m.0].fragment;
        let frag2 = &dag[m.1].fragment;
        let len = frag1.len() - 2;

        // Find indices of fragments that contain frag1 and frag2
        let (frag1_idx, _) = state.iter()
            .enumerate()   
            .find(|(_, state_frag)| frag1.is_subset(state_frag))
            .unwrap();
        let (frag2_idx, _) = state.iter()
            .enumerate()
            .find(|(_, state_frag)| frag2.is_subset(state_frag))
            .unwrap();

        // Union with masks to indicate that these edges are used
        masks[len][frag1_idx].union_with(frag1);
        masks[len][frag2_idx].union_with(frag2);
    }

    // Let new state be only the edges which can be matched
    let mut state = Vec::new();

    for frag in masks[0].iter() {
        // Split into connected components
        let connected_frags = connected_components_under_edges(mol.graph(), &frag);
        
        // Add connected components to new state
        state.extend(connected_frags);
    }
                

    // MATCHES METHOD 2
    // ONLY EXTEND IF A SUBGRAPH HAS A MATCH


    // Memoization
    if cache.memoize_state(mol, &state, state_index, &removal_order) {
        return (state_index, 1);
    }

    // Keep track of the best assembly index found in any of this assembly
    // state's children and the number of states searched, including this one.
    let best_child_index = AtomicUsize::from(state_index);
    let states_searched = AtomicUsize::from(1);
    
    // Define a closure that handles recursing to a new assembly state based on
    // the given (enumerated) pair of non-overlapping isomorphic subgraphs.
    let  recurse_on_match = |v: (usize, usize)| {
        //let (h1, h2) = &matches[v];
        let h1 = &dag[v.0].fragment;
        let h2 = &dag[v.1].fragment;
        
        // TODO: try keeping track of fragment to elimate some some search
        if let Some(fragments) = fragments(mol, &state, h1, h2) {
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
                matches,
                dag,
                BitSet::new(),
                {
                    let mut clone = removal_order.clone();
                    clone.push(*matches.get(&v).unwrap());
                    clone
                },
                &fragments,
                state_index - h1.len() + 1,
                best_index.clone(),
                bounds,
                &mut cache.clone(),
                new_parallel,
                kernel_mode,
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
        enum_matches
            .iter()
            .for_each(|v| recurse_on_match(*v));
    } else {
        enum_matches
            .par_iter()
            .for_each(|v| recurse_on_match(*v));
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
    let (matches, dag) = matches(mol, canonize_mode);

    /*let mut match_list: Vec<((usize, usize), usize)> = Vec::new();
    for m in matches.iter() {
        match_list.push((*m.0, *m.1));
    }

    for (idx, d) in dag.iter().enumerate() {
        println!("{}: {:?}", idx, d);
    }

    match_list.sort_by_key(|m| m.1);
    for m in match_list {
        println!("{:?}", m);
    }
    println!("{}", matches.len());*/

    // Create memoization cache.
    //let mut cache = Cache::new(memoize_mode, canonize_mode);
    let mut cache = NewCache::new(memoize_mode);

    /*let graph = {
        if clique {
            Some(CompatGraph::new(&matches, &dag))
        }
        else {
            None
        }
    };*/

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
    let (index, states_searched) = recurse_index_search(
        mol,
        &matches,
        &dag,
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
