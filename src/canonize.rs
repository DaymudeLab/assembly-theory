//! Create canonical labelings for molecular graphs.

use std::{
    collections::{BTreeSet, HashMap},
    hash::Hash,
    iter,
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
    utils::node_count_under_edge_mask,
};

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

    Tree(Vec<usize>),
}

/// Obtain a canonical labeling of the specified subgraph using the specified
/// algorithm.
pub fn canonize(mol: &Molecule, subgraph: &BitSet, mode: CanonizeMode) -> Labeling {
    match mode {
        CanonizeMode::Nauty => {
            let cgraph = subgraph_to_cgraph(mol, subgraph);
            Labeling::Nauty(CanonLabeling::new(&cgraph))
        }
        CanonizeMode::TreeNauty => {
            if is_tree(mol, subgraph) {
                Labeling::Tree(tree_canonize(mol, subgraph))
            } else {
                let cgraph = subgraph_to_cgraph(mol, subgraph);
                Labeling::Nauty(CanonLabeling::new(&cgraph))
            }
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

fn tree_canonize(mol: &Molecule, subgraph: &BitSet) -> Vec<usize> {
    let graph = mol.graph();
    let order = subgraph.len() + 1;
    let mut adjacencies = vec![BitSet::with_capacity(order); order];
    let mut partial_canonical_sets = vec![BTreeSet::<Vec<usize>>::new(); order];
    let mut unlabeled_vertices = BitSet::with_capacity(order);
    for ix in subgraph.iter() {
        let (u, v) = graph
            .edge_endpoints(EdgeIndex::new(ix))
            .expect("malformed bitset!");

        for node in [u, v] {
            let index = node.index();
            if unlabeled_vertices.contains(index) {
                continue;
            }
            let weight = graph.node_weight(node).unwrap();
            partial_canonical_sets[index].insert(vec![weight.element().into()]);
        }

        let (u, v) = (u.index(), v.index());
        adjacencies[u].insert(v);
        adjacencies[v].insert(u);
    }

    while unlabeled_vertices.len() > 2 {
        let leaves = adjacencies
            .iter()
            .enumerate()
            .filter_map(|(i, list)| (list.len() == 1).then_some(i))
            .collect::<Vec<_>>();

        for leaf in leaves {
            let parent = adjacencies[leaf].iter().next().unwrap();
            let edge = graph
                .edges_connecting(NodeIndex::new(parent), NodeIndex::new(leaf))
                .next()
                .unwrap();
            let mut canonical_label = vec![(*edge.weight()).into()];
            canonical_label.extend(partial_canonical_sets[leaf].iter().flatten());
            partial_canonical_sets[parent].insert(canonical_label);

            adjacencies[leaf].clear();
            adjacencies[parent].remove(leaf);
            unlabeled_vertices.remove(leaf);
        }
    }

    if unlabeled_vertices.len() == 2 {
        let mut iter = unlabeled_vertices.iter();
        let (u, v) = (iter.next().unwrap(), iter.next().unwrap());
        let edge = graph
            .edges_connecting(NodeIndex::new(u), NodeIndex::new(v))
            .next()
            .unwrap();
        let u_label = std::mem::take(&mut partial_canonical_sets[u])
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        let v_label = std::mem::take(&mut partial_canonical_sets[v])
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        let (first, second) = if u_label < v_label {
            (u_label, v_label)
        } else {
            (v_label, u_label)
        };
        iter::once((*edge.weight()).into())
            .chain(first)
            .chain(second)
            .collect()
    } else {
        let canonical_root = unlabeled_vertices.iter().next().unwrap();
        let canonical_set = std::mem::take(&mut partial_canonical_sets[canonical_root]);
        canonical_set.into_iter().flatten().collect()
    }
}

// Assuming subgraph represents a connected graph, check if it induces a tree
fn is_tree(mol: &Molecule, subgraph: &BitSet) -> bool {
    node_count_under_edge_mask(mol.graph(), subgraph) == subgraph.len() - 1
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
