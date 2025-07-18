//! Kernelize molecular graphs to improve top-down search efficiency.
//!
//! TODO: Longer explanation of what that means.

use clap::ValueEnum;
use bit_set::BitSet;
use crate::reductions::CompatGraph;

/// Graph kernelization strategy when searching using the clique reduction.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum KernelMode {
    /// No kernelization.
    None,
    /// Only kernelize the original molecule.
    Once,
    /// Kernelize the original molecule and the recursion's first level only.
    DepthOne,
    /// Perform kernelization at every recursive step.
    Always,
}

fn deletion_kernel(g: &CompatGraph, mut subgraph: BitSet) -> BitSet {
    let subgraph_copy = subgraph.clone();

    for v in subgraph_copy.iter() {
        let v_val = g.weight(v);
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

            let u_val = g.weight(u);
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

fn inclusion_kernel(g: &CompatGraph, subgraph: &BitSet) -> Vec<usize> {
    let mut kernel = Vec::new();
    let tot = subgraph.iter().map(|v| g.weight(v)).sum::<usize>();

    'outer: for v in subgraph {
        let v_val = g.weight(v);
        let neighbors_val = g.neighbors(v, subgraph).iter().map(|u| g.weight(u)).sum::<usize>();
        if v_val >= tot - neighbors_val - v_val { 
            kernel.push(v);
            continue;
        }

        let mut neighbors: Vec<usize> = vec![];

        for u in subgraph.difference(&g.compatible_with(v)) {
            if u == v {
                continue;
            }
            if g.weight(u) > v_val {
                continue 'outer;
            }

            for w in neighbors.iter() {
                if g.are_adjacent(u, *w) {
                    continue 'outer;
                }   
            }

            neighbors.push(u);
        }

        kernel.push(v);
    }

    kernel
}
