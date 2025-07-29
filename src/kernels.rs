//! Kernelize match-compatibility graphs to improve top-down search efficiency.
//!
//! The problem of computing the minimum assembly index of a molecule can be
//! reduced to finding the maximum weight clique in a compatibility graph over
//! matches (i.e., pairs of non-overlapping isomorphic subgraphs). Strucutral
//! properties of this graph can be used to determine match pairs (i.e., nodes)
//! that *definitely will* or *definitely won't* be used in an optimal
//! solution. We call the process of identifying these nodes *kernelization*.
//! Uses the strategies of neighborhood removal, isolated vertex removal, and
//! domination as described in Section 5.2 of [Lamm et al.
//! (2019)](https://doi.org/10.1137/1.9781611975499.12). (Note that they solve
//! the equivalent problem of weighted independent set.)

use clap::ValueEnum;
use bit_set::BitSet;
use crate::reductions::CompatGraph;

/// Graph kernelization strategy when searching using the clique reduction.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum KernelMode {
    /// Do not apply any kernelizations.
    None,
    /// Apply kernels only after the initial construction of the compatibility
    /// graph.
    Once,
    /// Apply kernels after the initial construction of the compability graph
    /// and again after any fragmentations of the full molecule.
    DepthOne,
    /// Apply kernels after every fragmentation step.
    Always,
}

/// Takes a subgraph as input and removes all nodes that will not be part of any maximum weight clique.
/// Uses the strategy of domination as described in Lamm et al.
pub fn deletion_kernel(matches: &Vec<(BitSet, BitSet)>, g: &CompatGraph, mut subgraph: BitSet) -> BitSet{
    let subgraph_copy = subgraph.clone();

    for v in subgraph_copy.iter() {
        let v_val = matches[v].0.len();
        let neighbors_v = g.neighbors(v, &subgraph);

        let Some(w1) = neighbors_v.iter().next() else {
            continue;
        };
        let Some(w2) = neighbors_v.iter().last() else {
            continue;
        };

        let mut s = subgraph.clone();
        s.intersect_with(&g.compatible_with(w1));
        for u in s.intersection(&&g.compatible_with(w2)) {
            if g.are_adjacent(v, u) || v == u {
                continue;
            }

            let u_val = matches[u].0.len();
            if v_val > u_val {
                continue;
            }

            let neighbors_u = g.neighbors(u, &subgraph);

            if neighbors_v.is_subset(&neighbors_u) {
                subgraph.remove(v);
                break;
            }
        }
    }

    subgraph
}

/// Takes a subgraph as input and returns the first vertex that will be included in some maximum weight clique.
/// Uses the strategies of neighborhood removal and isolated vertex removal from Lamm et al.
pub fn inclusion_kernel(matches: &Vec<(BitSet, BitSet)>, g: &CompatGraph, subgraph: &BitSet) -> usize {
    let tot = subgraph.iter().map(|v| matches[v].0.len() - 1).sum::<usize>();

    'outer: for v in subgraph {
        let v_val = matches[v].0.len() - 1;

        // Neighborhood removal
        let neighbors_val = g.neighbors(v, subgraph).iter().map(|u| matches[u].0.len() - 1).sum::<usize>();
        if v_val >= tot - neighbors_val - v_val { 
            return v;
        }

        // Isolated vertex removal.
        // Only makes it through this loop if the non-neighborhood of v is independent and
        // contains vertices with weight no higher than v.
        let mut neighbors: Vec<usize> = vec![];
        for u in subgraph.difference(&g.compatible_with(v)) {
            if u == v {
                continue;
            }
            if matches[u].0.len() - 1 > v_val {
                continue 'outer;
            }

            for w in neighbors.iter() {
                if g.are_adjacent(u, *w) {
                    continue 'outer;
                }   
            }

            neighbors.push(u);
        }

        return v;
    }

    matches.len()

}