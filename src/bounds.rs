//! Prune assembly states from which the assembly index cannot improve.
//!
//! There are two types of bounds: *state bounds* and *match bounds*. State
//! bounds are applied as soon as an assembly state is reached to determine
//! whether any children of this state may possibly yield an improved assembly
//! index; if not, this state is removed from the search tree. Otherwise, match
//! bounds are applied to exclude any matches (i.e., pairs of edge-disjoint,
//! isomorphic subgraphs) whose removal from the current assembly state cannot
//! not yield an improved assembly index.
//!
//! Each state bound considers the current assembly state (i.e., set of
//! fragments) and computes an upper bound on the "savings" (in terms of number
//! of joining operations) that can possibly be obtained when constructing the
//! molecule using this state's fragments and subfragments thereof. Let
//! `state_index` be this assembly state's assembly index, `best_index` be the
//! smallest assembly index found across any assembly state so far, and `bound`
//! be the upper bound on this assembly state's possible savings. If ever
//! `state_index` - `bound` >= `best_index`, then no descendant of this
//! assembly state can possibly yield an assembly index better than
//! `best_index` and thus this assembly state can be pruned.
//!
//! Each match bound considers both the current assembly state and the matches
//! (i.e., pairs of edge-disjoint, isomorphic subgraphs) that can be removed
//! from it. It determines a size (in edges) such that the removal of any match
//! smaller than this size cannot possibly yield a better assembly index. This
//! size is then used as a filter when generating matches to remove from a
//! given assembly state in [`crate::matches::Matches::matches_to_remove`].

use bit_set::BitSet;
use clap::ValueEnum;

use crate::{
    object::{AObject, NodeBehavior, EdgeBehavior},
    state::State,
};

/// State and match bounds to prune assembly index search.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Bound {
    /// The shortest number of joining operations to create a molecule with |E|
    /// bonds is log_2(|E|), i.e., if it is possible to always join the largest
    /// fragment with itself to produce the molecule. Thus, an upper bound on
    /// a state's savings is [#fragment bonds] - log_2([#fragment bonds]); see
    /// [Jirasek et al. (2024)](https://doi.org/10.1021/acscentsci.4c00120).
    Log,
    /// An improvement over `Log` that also uses the size of the "largest
    /// duplicatable subgraph" for this state in an integer addition chain; see
    /// [Seet et al. (2024)](https://arxiv.org/abs/2410.09100).
    Int,
    /// Uses the types of bonds in the molecule to bound the number of assembly
    /// steps remaining. The first time a unique bond type is added to the
    /// graph, it could not have been part of a duplicate since that bond type
    /// has not been used yet. Thus the number of unique bond types gives
    /// information on how many more joins are required.
    VecSimple,
    /// Considers the fragments of size two in the current fragmentation. In
    /// the remaining top-down process, such fragments will require one step to
    /// remove if there is a duplicate set of two bonds in the graph.
    /// Otherwise, they will require two steps.
    VecSmallFrags,
    /// Like `Int`, but uses only matchable edges as computed in
    /// [`crate::matches::Matches::matches_to_remove`].
    MatchableEdges,
}

/// Edge information used in vector addition chain bounds.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeType {
    kind: u8,
    ends: (u8, u8),
}

/// Return `true` iff any of the given state bounds prune this assembly state.
pub(crate) fn state_bounds<T: AObject>(
    mol: &T,
    state: &State,
    best_index: usize,
    bounds: &[Bound],
) -> bool {
    let fragments = state.fragments();
    let state_index = state.index();
    let largest_removed = state.largest_removed();

    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => state_index - log_bound(fragments) >= best_index,
            Bound::Int => state_index - int_bound(fragments, largest_removed) >= best_index,
            Bound::VecSimple => {
                state_index - vec_simple_bound(fragments, largest_removed, mol) >= best_index
            }
            Bound::VecSmallFrags => {
                state_index - vec_small_frags_bound(fragments, largest_removed, mol) >= best_index
            }
            _ => false,
        };
        if exceeds {
            return true;
        }
    }
    false
}

/// Compute a lower bound on the size of the largest match whose removal from
/// the given assembly state could possibly yield an improved assembly index.
pub(crate) fn match_bounds(
    state_index: usize,
    best_index: usize,
    matchable_edge_masks: &[Vec<BitSet>],
    bounds: &[Bound],
) -> usize {
    // Test different match removal sizes.
    for removal_size in 2..matchable_edge_masks.len() + 2 {
        let mut bounded = false;

        for bound_type in bounds {
            bounded |= match bound_type {
                Bound::MatchableEdges => {
                    state_index
                        >= usable_edges_bound(matchable_edge_masks, removal_size) + best_index
                }
                _ => false,
            };
            if bounded {
                break;
            }
        }

        // If no bounds apply, the largest match size has been found.
        if !bounded {
            return removal_size;
        }
    }

    matchable_edge_masks.len() + 2
}

/// TODO
fn log_bound(fragments: &[BitSet]) -> usize {
    let mut size = 0;
    for f in fragments {
        size += f.len();
    }

    size - (size as f32).log2().ceil() as usize
}

/// TODO
fn int_bound(fragments: &[BitSet], m: usize) -> usize {
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

/// TODO
// Count number of unique edges in a fragment
// Helper function for vector bounds
fn unique_edges<T: AObject>(fragment: &BitSet, mol: &T) -> Vec<EdgeType> {
    let g = mol.graph();
    let mut nodes: Vec<u8> = Vec::new();
    for v in g.node_weights() {
        nodes.push(v.kind());
    }
    let edges: Vec<petgraph::prelude::EdgeIndex> = g.edge_indices().collect();
    let weights: Vec<T::EdgeLabel> = g.edge_weights().copied().collect();

    // types will hold an element for every unique edge type in fragment
    let mut types: Vec<EdgeType> = Vec::new();
    for idx in fragment.iter() {
        let bond = weights[idx];
        let e = edges[idx];

        let (e1, e2) = g.edge_endpoints(e).expect("bad");
        let e1 = nodes[e1.index()];
        let e2 = nodes[e2.index()];
        let ends = if e1 < e2 { (e1, e2) } else { (e2, e1) };

        let edge_type = EdgeType { kind: bond.kind(), ends };

        if types.contains(&edge_type) {
            continue;
        } else {
            types.push(edge_type);
        }
    }

    types
}

/// TODO
fn vec_simple_bound<T: AObject>(fragments: &[BitSet], m: usize, mol: &T) -> usize {
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

/// TODO
fn vec_small_frags_bound<T: AObject>(fragments: &[BitSet], m: usize, mol: &T) -> usize {
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

/// TODO
fn usable_edges_bound(matchable_edge_masks: &[Vec<BitSet>], removal_size: usize) -> usize {
    let mut bound = 0;

    for (frag_ix, frag) in matchable_edge_masks[removal_size - 2].iter().enumerate() {
        let total_removable_edges = matchable_edge_masks[0][frag_ix].len();
        let removable_edges = frag.len();
        let leftover_edges =
            (total_removable_edges - removable_edges) + (removable_edges % removal_size);
        bound += total_removable_edges
            - (removable_edges / removal_size)
            - leftover_edges.div_ceil(removal_size - 1);
    }

    bound.saturating_sub((removal_size as f32).log2().ceil() as usize)
}
