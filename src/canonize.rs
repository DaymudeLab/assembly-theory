//! Create canonical labelings for molecular graphs.

use std::{
    collections::HashMap,
};

use bit_set::BitSet;
use clap::ValueEnum;
use graph_canon::CanonLabeling;
use petgraph::{
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::{
    molecule::{AtomOrBond, Index, Molecule},
};

/// Algorithm for graph canonization.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum CanonizeMode {
    /// Use the Nauty algorithm of McKay & Piperno (2014).
    Nauty,
    /// Use the algorithm of Faulon et al. (2004).
    Faulon,
    /// Use a fast tree canonization algorithm if applicable; else use Nauty.
    TreeNauty,
    /// Use a fast tree canonization algorithm if applicable; else use Faulon.
    TreeFaulon,
}

/// Obtain a canonical labeling of the specified `subgraph` using the
/// algorithm specified by `mode`.
// TODO: Should return a Box<dyn Hash> to flexibly extend over different
// canonization algorithms.
pub fn canonize(mol: &Molecule, subgraph: &BitSet, mode: CanonizeMode)
    -> CanonLabeling<AtomOrBond> {
    match mode {
        CanonizeMode::Nauty => {
            let cgraph = subgraph_to_cgraph(mol, subgraph);
            CanonLabeling::new(&cgraph)
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

