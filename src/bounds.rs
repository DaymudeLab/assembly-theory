//! TODO

use bit_set::BitSet;
use clap::ValueEnum;

use crate::{
    molecule::{Bond, Element, Molecule},
    reductions::{CompatGraph},
};

/// Bounding strategies for the search phase of assembly index calcluation.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Bound {
    /// Trivial bound of log_2(# remaining bonds/edges).
    Log,
    /// Bound using the length of the shortest integer addition chain defined
    /// using fragment sizes.
    Int,
    /// Bound using the length of the shortest vector addition chain defined
    /// using fragments' number and types of edges.
    VecSimple,
    /// Bound using the length of the shortest vector addition chain defined
    /// using information about the molecule's number of fragments of size 2.
    VecSmallFrags,
    /// TODO
    CoverSort,
    /// TODO
    CoverNoSort,
    /// TODO
    CliqueBudget,
}

/// Edge information used in vector addition chain bounds.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeType {
    bond: Bond,
    ends: (Element, Element),
}

/// TODO
pub fn bound_exceeded(
    mol: &Molecule,
    fragments: &[BitSet],
    ix: usize,
    best: usize,
    largest_remove: usize,
    bounds: &[Bound],
    matches_graph: &CompatGraph,
    subgraph: &BitSet,
) -> bool {
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::Int => ix - int_bound(fragments, largest_remove) >= best,
            Bound::VecSimple => ix - vec_simple_bound(fragments, largest_remove, mol) >= best,
            Bound::VecSmallFrags => {
                ix - vec_small_frags_bound(fragments, largest_remove, mol) >= best
            }
            Bound::CliqueBudget => ix - clique_budget_bound(matches_graph, subgraph, fragments) >= best,
            Bound::CoverNoSort => ix - cover_bound(matches_graph, subgraph, false),
            Bound::CoverSort => ix - cover_bound(matches_graph, subgraph, true),
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
pub fn log_bound(fragments: &[BitSet]) -> usize {
    let mut size = 0;
    for f in fragments {
        size += f.len();
    }

    size - (size as f32).log2().ceil() as usize
}

/// TODO
pub fn int_bound(fragments: &[BitSet], m: usize) -> usize {
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
pub fn vec_simple_bound(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
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
pub fn vec_small_frags_bound(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
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

pub fn remaining_weight_bound(graph: &CompatGraph, subgraph: &BitSet) -> usize {
    let deg_sum = subgraph.iter().map(|v| graph.degree(v, subgraph)).sum::<usize>() as f32;
    let max_clique = ((1_f32 + (4_f32 * deg_sum + 1_f32).sqrt()) / 2_f32).floor() as usize;
    let mut sum = 0;
    let mut iter = subgraph.iter();
    for _ in 0..max_clique {
        sum += graph.weight(iter.next().unwrap());
    };

    sum
}

pub fn color_bound(graph: &CompatGraph, subgraph: &BitSet) -> usize{
    // Greedy coloring
    let mut colors: Vec<i32> = vec![-1; graph.len()];
    let mut num_colors = 0;
    let mut largest: Vec<usize> = Vec::new();
    

    for v in (0..graph.len()).rev() {
        if !subgraph.contains(v) {
            continue;
        }

        let mut used: Vec<usize> = vec![0; num_colors];

        for u in subgraph.intersection(graph.compatible_with(v)) {
            if colors[u] != -1 {
                used[colors[u] as usize] = 1;
            }
        }

        let mut max = 0;
        let mut max_idx = num_colors;
        for i in 0..num_colors {
            if used[i] == 0 && largest[i] > max {
                max = largest[i];
                max_idx = i;
            }
        }

        if max_idx == num_colors {
            num_colors += 1;
            largest.push(0);
        }
        if graph.weight(v) > largest[max_idx] {
            largest[max_idx] = graph.weights(v)
        }

        colors[v] = max_idx as i32;
    }

    largest.iter().sum::<usize>()
}

pub fn cover_bound(graph: &CompatGraph, subgraph: &BitSet, sort: bool) -> usize {
    // Sort vertices
    if sort {
        let mut vertices: Vec<(usize, usize)> = Vec::with_capacity(subgraph.len());
        for v in subgraph {
            vertices.push((v, graph.degree(v, subgraph)));
        }   
        vertices.sort_by(|a, b| b.1.cmp(&a.1));
        cover_bound_helper(graph, subgraph, vertices.iter().map(|(v, _)| *v))
    }
    else {
        let vertices = (0..graph.len()).rev().filter(|v| subgraph.contains(*v));
        cover_bound_helper(graph, subgraph, vertices)
    }
}

fn cover_bound_helper(graph: &CompatGraph, subgraph: &BitSet, iter: impl Iterator<Item = usize>) -> usize {
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
        let v_val = graph.weight(v);
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

pub fn clique_budget_bound(graph: &CompatGraph, subgraph: &BitSet, fragments: &[BitSet]) -> usize {
    let total_bonds = fragments.iter().map(|x| x.len()).sum::<usize>();
    let mut bound = 0;
    let sizes = {
        let mut vec = vec![];
        let mut prev = 0;
        for s in subgraph.iter().map(|v| graph.weight(v) + 1) {
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
            if graph.weight(v) + 1 > i {
                continue;
            }


            let dup = &graph.matches(v);
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
