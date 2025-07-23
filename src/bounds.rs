//! Prune assembly states from which the assembly index cannot improve.
//!
//! Each bound takes information about the current assembly state (i.e., set of
//! fragments) and computes an upper bound on the "savings" (in terms of number
//! of joining operations) that can possibly be obtained when constructing the
//! molecule using this state's fragments and subfragments thereof. Let
//! `state_index` be this assembly state's assembly index, `best_index` be the
//! smallest assembly index found across any assembly state so far, and `bound`
//! be the upper bound on this assembly state's possible savings. If ever
//! `state_index` - `bound` >= `best_index`, then no descendant of this
//! assembly state can possibly yield an assembly index better than
//! `best_index` and thus this assembly state can be pruned.

use bit_set::BitSet;
use clap::ValueEnum;

use crate::molecule::{Bond, Element, Molecule};

/// Type of upper bound on the "savings" possible from an assembly state.
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
    /// A weighted independent set cover provides a bound on the size of a max.
    /// weight clique in the compatibility graph. Uses a greedy algorithm  to
    /// construct such a cover and obtain a bound. See
    /// [Lamm et al. (2019)](https://doi.org/10.1137/1.9781611975499.12) for
    /// the definition of a cover. (Note that they solve the equivalent
    /// weighted independent set problem and thus use a clique cover instead.)
    CoverNoSort,
    /// Like `CoverNoSort`, buts sorts the vertices of the compatibility
    /// graph by degree before creating the greedy independent set cover.
    CoverSort,
    /// Uses the compatibility graph to determine the largest duplicatable
    /// subraphs remaining in each fragment. Uses this to bound the best
    /// possible savings obtainable for each fragment.
    CliqueBudget,
}

/// Edge information used in vector addition chain bounds.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeType {
    bond: Bond,
    ends: (Element, Element),
}

/// Returns `true` iff any of the given bounds would prune this assembly state.
pub fn bound_exceeded(
    mol: &Molecule,
    fragments: &[BitSet],
    state_index: usize,
    best_index: usize,
    largest_remove: usize,
    bounds: &[Bound],
) -> bool {
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => state_index - log_bound(fragments) >= best_index,
            Bound::Int => state_index - int_bound(fragments, largest_remove) >= best_index,
            Bound::VecSimple => {
                state_index - vec_simple_bound(fragments, largest_remove, mol) >= best_index
            }
            Bound::VecSmallFrags => {
                state_index - vec_small_frags_bound(fragments, largest_remove, mol) >= best_index
            }
            _ => {
                panic!("One of the chosen bounds is not implemented yet!")
            }
        };
        if exceeds {
            return true;
        }
    }
    false
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

        if types.contains(&edge_type) {
            continue;
        } else {
            types.push(edge_type);
        }
    }

    types
}

/// TODO
fn vec_simple_bound(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
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
fn vec_small_frags_bound(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
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
