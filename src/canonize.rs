//! Create canonical labelings for molecular graphs.

use std::{collections::HashMap, hash::Hash};

use bit_set::BitSet;
use clap::ValueEnum;
use graph_canon::CanonLabeling;
use petgraph::{
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::molecule::{AtomOrBond, Index, Molecule};

/// Algorithm for graph canonization.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum CanonizeMode {
    /// Use the Nauty algorithm of [McKay & Piperno
    /// (2014)](https://doi.org/10.1016/j.jsc.2013.09.003).
    Nauty,
    /// Use the algorithm of
    /// [Faulon et al. (2004)](https://doi.org/10.1021/ci0341823).
    Faulon,
    /// Use a tree canonization algorithm if applicable; else use `Nauty`.
    TreeNauty,
    /// Use a tree canonization algorithm if applicable; else use `Faulon`.
    TreeFaulon,
}

/// Canonical labeling returned by our graph canonization functions.
#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub enum Labeling {
    /// Labeling returned by the `graph_canon` crate.
    Nauty(CanonLabeling<AtomOrBond>),

    /// A string labeling returned by our implementation of
    /// [Faulon et al. (2004)](https://doi.org/10.1021/ci0341823).
    // TODO: This should be a `Vec<u8>`
    Faulon(String),
}

/// Obtain a canonical labeling of the specified subgraph using the specified
/// algorithm.
pub fn canonize(mol: &Molecule, subgraph: &BitSet, mode: CanonizeMode) -> Labeling {
    match mode {
        CanonizeMode::Nauty => {
            let cgraph = subgraph_to_cgraph(mol, subgraph);
            Labeling::Nauty(CanonLabeling::new(&cgraph))
        }
        _ => {
            panic!("The chosen --canonize mode is not implemented yet!")
        }
    }
}

/// A graph representation interpretable by Nauty.
type CGraph = Graph<AtomOrBond, (), Undirected, Index>;

/// Convert the specified `subgraph` to the format expected by Nauty.
fn subgraph_to_cgraph(mol: &Molecule, subgraph: &BitSet) -> CGraph {
    let mut h = CGraph::with_capacity(subgraph.len(), 2 * subgraph.len());
    let mut vtx_map = HashMap::<NodeIndex, NodeIndex>::new();
    for e in subgraph {
        let eix = EdgeIndex::new(e);
        let (src, dst) = mol.graph().edge_endpoints(eix).unwrap();
        let src_w = mol.graph().node_weight(src).unwrap();
        let dst_w = mol.graph().node_weight(dst).unwrap();
        let e_w = mol.graph().edge_weight(eix).unwrap();

        let h_enode = h.add_node(AtomOrBond::Bond(*e_w));

        let h_src = vtx_map
            .entry(src)
            .or_insert(h.add_node(AtomOrBond::Atom(*src_w)));
        h.add_edge(*h_src, h_enode, ());

        let h_dst = vtx_map
            .entry(dst)
            .or_insert(h.add_node(AtomOrBond::Atom(*dst_w)));
        h.add_edge(*h_dst, h_enode, ());
    }
    h
}

mod tests {
    #[allow(unused_imports)]
    use super::*;

    #[allow(unused_imports)]
    use petgraph::algo::is_isomorphic_matching;

    #[test]
    fn noncanonical() {
        let mut p3_010 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = p3_010.add_node(0);
        let n1 = p3_010.add_node(1);
        let n2 = p3_010.add_node(0);
        p3_010.add_edge(n0, n1, ());
        p3_010.add_edge(n1, n2, ());

        let mut p3_001 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = p3_001.add_node(0);
        let n1 = p3_001.add_node(0);
        let n2 = p3_001.add_node(1);
        p3_001.add_edge(n0, n1, ());
        p3_001.add_edge(n1, n2, ());

        let repr_a = CanonLabeling::new(&p3_010);
        let repr_b = CanonLabeling::new(&p3_001);

        assert_ne!(repr_a, repr_b);
    }

    #[test]
    fn nonisomorphic() {
        let mut p3_010 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = p3_010.add_node(0);
        let n1 = p3_010.add_node(1);
        let n2 = p3_010.add_node(0);
        p3_010.add_edge(n0, n1, ());
        p3_010.add_edge(n1, n2, ());

        let mut p3_001 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = p3_001.add_node(0);
        let n1 = p3_001.add_node(0);
        let n2 = p3_001.add_node(1);
        p3_001.add_edge(n0, n1, ());
        p3_001.add_edge(n1, n2, ());

        assert!(!is_isomorphic_matching(
            &p3_001,
            &p3_010,
            |e0, e1| e0 == e1,
            |n0, n1| n0 == n1
        ))
    }
}
