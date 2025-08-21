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

use crate::{
    molecule::{Bond, Element, Molecule},
    reductions::CompatGraph,
};

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
    matches: &Vec<(BitSet, BitSet)>,
    graph: &Option<CompatGraph>,
    fragments: &[BitSet],
    subgraph: &BitSet,
    state_index: usize,
    best_index: usize,
    largest_remove: usize,
    bounds: &[Bound],
    masks: &Vec<Vec<BitSet>>,
) -> bool {
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => state_index - log_bound(fragments) >= best_index,
            //Bound::Int => state_index - int_bound_seet(masks, largest_remove) >= best_index,
            Bound::Int => best_index + int_bound_new(masks) <= state_index,
            Bound::VecSimple => {
                state_index - vec_simple_bound(fragments, largest_remove, mol) >= best_index
            }
            Bound::VecSmallFrags => {
                state_index - vec_small_frags_bound(fragments, largest_remove, mol) >= best_index
            }
            Bound::CliqueBudget => best_index + clique_budget_bound(matches, subgraph, fragments) <= state_index,
            Bound::CoverNoSort => {
                if let Some(g) = graph {
                    best_index + cover_bound(matches, g, subgraph, false) <= state_index
                }
                else {
                    false
                }
            }
            Bound::CoverSort => {
                if let Some(g) = graph {
                    best_index + cover_bound(matches, g, subgraph, true) <= state_index
                }
                else {
                    false
                }
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
fn _int_bound(fragments: &[BitSet], m: usize) -> usize {
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

fn int_bound_new(masks: &Vec<Vec<BitSet>>) -> usize {
    let mut max = 0;
    for (i, list) in masks.iter().enumerate() {
        let mut val = 0;
        let i = i + 2;
        for (j, frag) in list.iter().enumerate() {
            let mf = masks[0][j].len();
            let mi = frag.len();
            let x = mf + (mi % i) - mi;
            val += mf - (mi / i) - x / (i-1) - (x % (i-1) != 0) as usize;
        }
        val -= (i as f32).log2().ceil() as usize;
        max = max.max(val);
    }

    max
}

fn _int_bound_seet(masks: &Vec<Vec<BitSet>>, largest_remove: usize) -> usize {
    let mut dup_bonds_2 = 0;
    let mut dup_bonds_total;
    let mut max;
    let mut size_lists = vec![Vec::new(); masks.len()];

    for j in 0..masks.len()
    {
        size_lists[j] = vec![0; masks[j].len()];
        for i in 0..masks[0].len()
        {
            size_lists[j][i] = masks[j][i].len();
        }
    }

    for i in 0..size_lists[0].len() {
        dup_bonds_2 += size_lists[0][i] / 2;
    }
    dup_bonds_2 -= 1;
    max = dup_bonds_2;

    for j in 3..largest_remove+1
    {
        let size_list = &size_lists[j - 2];
        let size_list2 = &size_lists[0];
        let mut adjusted_size_list = vec![0; size_list.len()];
        let mut adjusted_size_list2 = vec![0; size_list.len()];
        for i in 0..size_list.len()
        {
            adjusted_size_list[i] = size_list[i] - size_list[i] % j;
            adjusted_size_list2[i] = size_list2[i] - adjusted_size_list[i];
        }
        dup_bonds_total = 0;
        for i in 0..size_list.len()
        {
            dup_bonds_total += adjusted_size_list[i] - adjusted_size_list[i] / j;
            dup_bonds_total += adjusted_size_list2[i] - adjusted_size_list2[i] / (j - 1);
            if adjusted_size_list2[i] % (j - 1) != 0 {
                dup_bonds_total -= 1;
            }
        }
        dup_bonds_total -= (j as f32).log2().ceil() as usize;
        max = max.max(dup_bonds_total);
    }

    max
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

pub fn clique_budget_bound(matches: &Vec<(BitSet, BitSet)>, subgraph: &BitSet, fragments: &[BitSet]) -> usize {
    let total_bonds = fragments.iter().map(|x| x.len()).sum::<usize>();
    let mut bound = 0;
    let sizes = {
        let mut vec = vec![];
        let mut prev = 0;
        for s in subgraph.iter().map(|v| matches[v].0.len()) {
            if s != prev {
                vec.push(s);
                prev = s;
            }
        }

        vec
    };
    
    for i in sizes {
        let mut bound_temp = 0;
        let mut has_bonds = fragments.len();
        let mut num_bonds: Vec<usize> = fragments.iter().map(|x| x.len()).collect();
        let mut smallest_remove = i;

        for v in subgraph.iter() {
            if has_bonds == 0 {
                break;
            }
            if matches[v].0.len() > i {
                continue;
            }


            let dup = &matches[v];
            let bond = dup.1.iter().next().unwrap();
            let mut j = 0;
            while !fragments[j].contains(bond) {
                j += 1;
            }

            if num_bonds[j] > 0 {
                let remove = std::cmp::min(dup.0.len(), num_bonds[j]);
                bound_temp += 1;
                num_bonds[j] -= remove;
                smallest_remove = std::cmp::min(smallest_remove, remove);

                if num_bonds[j] == 0 {
                    has_bonds -= 1;
                }
            }
        }

        let leftover = num_bonds.iter().sum::<usize>();
        let log = {
            if leftover > 0 {
                0
            }
            else {
                (smallest_remove as f32).log2().ceil() as usize
            }
        };
        bound = std::cmp::max(bound, total_bonds - bound_temp - leftover - log);
    }

    bound
}

pub fn cover_bound(matches: &Vec<(BitSet, BitSet)>, graph: &CompatGraph, subgraph: &BitSet, sort: bool) -> usize {
    // Sort vertices
    if sort {
        let mut vertices: Vec<(usize, usize)> = Vec::with_capacity(subgraph.len());
        for v in subgraph {
            vertices.push((v, graph.degree(v, subgraph)));
        }   
        vertices.sort_by(|a, b| b.1.cmp(&a.1));
        cover_bound_helper(matches, graph, subgraph, vertices.iter().map(|(v, _)| *v))
    }
    else {
        let vertices = (0..graph.len()).rev().filter(|v| subgraph.contains(*v));
        cover_bound_helper(matches, graph, subgraph, vertices)
    }
}

fn cover_bound_helper(matches: &Vec<(BitSet, BitSet)>, graph: &CompatGraph, subgraph: &BitSet, iter: impl Iterator<Item = usize>) -> usize {
    let mut colors: Vec<Option<Vec<usize>>> = vec![None; graph.len()];
    let mut col_weights = vec![];
    let mut num_col = 0;

    for v in iter {
        let mut v_col = Vec::new();
        let mut used = vec![0; num_col];

        // Find colors used in neighborhood of v
        for u in subgraph.intersection(graph.compatible_with(v)) {
            let Some(u_col) = &colors[u] else {
                continue;
            };

            for c in u_col {
                used[*c] = 1;
            }
        }

        let mut total_weight = 0;
        let v_val = matches[v].0.len() - 1;
        // Find colors to give to v
        for c in 0..num_col {
            if used[c] == 1 {
                continue;
            }

            v_col.push(c);
            total_weight += col_weights[c];

            if total_weight >= v_val {
                break;
            }
        }

        if total_weight == 0 {
            v_col.push(num_col);
            col_weights.push(v_val);
            num_col += 1
        }
        else if total_weight < v_val {
            let mut k = num_col - 1;
            while used[k] == 1 {
                k -= 1
            }
            col_weights[k] += v_val - total_weight
        }

        colors[v] = Some(v_col);
    };

    col_weights.iter().sum()
}
