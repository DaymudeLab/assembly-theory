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
    }, time::{Duration, Instant}
};

use bit_set::BitSet;
use clap::ValueEnum;
use graph_canon::canon;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator, IndexedParallelIterator};
use petgraph::graph::EdgeIndex;
use itertools::Itertools;
use std::cmp::{max, min};
use ahash::RandomState;
use fxhash::FxHashMap;

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

/// Helper function for [`index_search`]; only public for benchmarking.
///
/// Return all pairs of non-overlapping, isomorphic subgraphs in the molecule,
/// sorted to guarantee deterministic iteration.
pub fn matches(mol: &Molecule, canonize_mode: CanonizeMode) -> (FxHashMap<(usize, usize), usize>, Vec<DagNode>) {
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
                        let parent = &mut dag[parent_idx];
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
                        let parent = &mut dag[parent_idx];
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
    let mut match_to_id: FxHashMap<(usize, usize), usize> = FxHashMap::default();
    for m in matches {
        match_to_id.insert(m, next_match_id);
        next_match_id += 1;
    }

    (match_to_id, dag)
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
    matches: &FxHashMap<(usize, usize), usize>,
    dag: &Vec<DagNode>,
    last_removed: usize,
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
    
    // Generate matches
    let num_edges = mol.graph().edge_count();

    let mut valid_matches: Vec<(usize, usize)> = Vec::new();
    let mut subgraphs: Vec<(usize, usize)> = Vec::new();
    let mut masks: Vec<Vec<BitSet>> = Vec::new();

    // Create subgraphs of size 1
    for (state_id, fragment) in state.iter().enumerate() {
        for edge_idx in fragment.iter() {
            subgraphs.push((edge_idx, state_id));
        }
    }

    while subgraphs.len() > 0 {
        // Subgraphs to extend in next loop
        let mut buckets: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let mut new_subgraphs = Vec::new();
        let mut used_edge_mask = vec![BitSet::with_capacity(num_edges); state.len()];

        // Extend each subgraph and bucket
        for (frag_id, state_id) in subgraphs {
            // Get state frament that this subgraph is contained in
            let state_frag = &state[state_id];

            // Get ids of extensions of this subgraph
            let children_ids = &dag[frag_id].children;

            for child_id in children_ids {
                // Check if this extension is valid
                // i.e. the extended subgraph is contained in a fragment of state
                let child_frag = &dag[*child_id].fragment;
                let valid = child_frag.is_subset(state_frag);

                if valid {
                    let canon_id = &dag[*child_id].canon_id;
                    
                    // Add fragment to bucket
                    buckets.entry(*canon_id)
                        .and_modify(|bucket| bucket.push((*child_id, state_id)))
                        .or_insert(vec![(*child_id, state_id)]);
                }
            }
        }

        // Search through buckets and create matches
        for fragments in buckets.values() {
            let mut has_match = BitSet::with_capacity(fragments.len());
            // Loop over pairs of fragments
            for i in 0..fragments.len() {
                for j in i+1..fragments.len() {
                    let mut frag1 = (i, fragments[i]);
                    let mut frag2 = (j, fragments[j]);

                    // Order fragments
                    if frag1.1.0 < frag2.1.0 {
                        let temp = frag1;
                        frag1 = frag2;
                        frag2 = temp;
                    }

                    let frag1_bucket_idx = frag1.0;
                    let frag1_id = frag1.1.0;
                    let frag1_state_id = frag1.1.1;

                    let frag2_bucket_idx = frag2.0;
                    let frag2_id = frag2.1.0;
                    let frag2_state_id = frag2.1.1;

                    // Check for valid match
                    // TODO: try alternatives to this
                    if let Some(x) = matches.get(&(frag1_id, frag2_id)) {
                        // Only add match if it occurs later in the ordering than
                        // the last removed match
                        if *x >= last_removed {
                            valid_matches.push((frag1_id, frag2_id));

                            // If this is the first time seeing that this frag has a match,
                            // add it to the new_subgraphs list and modify mask
                            if !has_match.contains(frag1_bucket_idx) {
                                new_subgraphs.push((frag1_id, frag1_state_id));
                                used_edge_mask[frag1_state_id].union_with(&dag[frag1_id].fragment);
                                has_match.insert(frag1_bucket_idx);
                            }
                            if !has_match.contains(frag2_bucket_idx) {
                                new_subgraphs.push((frag2_id, frag2_state_id));
                                used_edge_mask[frag2_state_id].union_with(&dag[frag2_id].fragment);
                                has_match.insert(frag2_bucket_idx);
                            }
                        }
                    }
                }
            }
        }

        // Add mask to list
        masks.push(used_edge_mask);

        // Update to new subgraphs
        subgraphs = new_subgraphs;
    }

    // If there are no matches, return
    if valid_matches.len() == 0 {
        return (state_index, 1);
    }

    valid_matches.sort();
    valid_matches.reverse();

    // Let new state be only the edges which can be matched
    let mut state = Vec::new();

    for frag in masks[0].iter() {
        // Add connected components to new state
        state.extend(connected_components_under_edges(mol.graph(), &frag));
    }

    // Remove frags of size 1 or less
    state.retain(|frag| frag.len() >= 2);

    let best = best_index.load(Relaxed);
    let mut best_bound = 0;
    let mut largest_length = 2;

    // Use bounding strategy to find the largest length match
    // to try removing from this state.
    for (i, list) in masks.iter().enumerate() {
        let mut bound = 0;
        let i = i + 2;
        for (j, frag) in list.iter().enumerate() {
            let mf = masks[0][j].len();
            let mi = frag.len();
            let x = mf + (mi % i) - mi;
            bound += mf - (mi / i) - x / (i-1) - (x % (i-1) != 0) as usize;
        }
        bound -= (i as f32).log2().ceil() as usize;

        best_bound = max(best_bound, bound);

        // Check if removing at this length can give a more optimal answer
        // If yes, stop and return the largest_length to remove
        if state_index - best_bound < best {
            break;
        }

        largest_length += 1;
    }

    // Remove from valid_matches any matches that have size smaller
    // than largest_length, and thus their removal does not need to be tested
    if largest_length > 2 {
        let mut final_idx = 0;
        while dag[valid_matches[final_idx].0].fragment.len() >= largest_length {
            final_idx += 1;
        }
        valid_matches.truncate(final_idx);
    }


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
    let  recurse_on_match = |i: usize, v: (usize, usize)| {
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
                *matches.get(&v).unwrap(),
                BitSet::new(),
                {
                    let mut clone = removal_order.clone();
                    // clone.push(*matches.get(&v).unwrap());
                    clone.push(i);
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
    kernel_mode: KernelMode,
    bounds: &[Bound],
    clique: bool,
) -> (u32, u32, usize) {
    // Catch not-yet-implemented modes.

    // Enumerate non-overlapping isomorphic subgraph pairs.
    let (matches, dag) = matches(mol, canonize_mode);

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
        0,
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
