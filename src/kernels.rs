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