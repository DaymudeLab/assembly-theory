//! Create canonical labelings for molecular graphs.

use std::{collections::HashMap, hash::Hash, iter, mem};

use bit_set::BitSet;
use clap::ValueEnum;
use petgraph::{
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::{
    molecule::{AtomOrBond, Index, Molecule},
    nauty::CanonLabeling,
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
#[derive(Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Labeling {
    /// Labeling returned by the `graph_canon` crate.
    Nauty(CanonLabeling<AtomOrBond>),

    /// A string labeling returned by our implementation of
    /// [Faulon et al. (2004)](https://doi.org/10.1021/ci0341823).
    // TODO: This should be a `Vec<u8>`
    Faulon(String),

    /// Labeling returned by our implementation of tree canonization.
    Tree(Vec<u8>),
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

/// Returns `True` iff the *connected* `subgraph` induces a tree.
fn is_tree(mol: &Molecule, subgraph: &BitSet) -> bool {
    node_count_under_edge_mask(mol.graph(), subgraph) == subgraph.len() + 1
}

/// Wrap a bytevec with 0..vec..255.
fn wrap_with_delimiters(data: Vec<u8>) -> impl Iterator<Item = u8> {
    iter::once(u8::MIN).chain(data).chain(iter::once(u8::MAX))
}

/// Sort subvectors and collapse into delimited sequence. For small n, vec+sort
/// is faster than delimiters.
// TODO: Swap with a radix sort.
fn collapse_set(mut set: Vec<Vec<u8>>) -> Vec<u8> {
    set.sort_unstable();
    set.into_iter().flat_map(wrap_with_delimiters).collect()
}

/// Obtain a canonical labeling of a `subgraph` inducing a tree using an
/// algorithm inspired by Aho, Hopcroft, and Ullman; see
/// [Read (1972)](https://doi.org/10.1016/B978-1-4832-3187-7.50017-9) and
/// [Ingels (2024)](https://arxiv.org/abs/2309.14441).
///
/// Two isomorphic trees will have the same canonical labeling. Assumes
/// `subgraph` induces a tree; this should be checked before calling this
/// function with `is_tree`.
fn tree_canonize(mol: &Molecule, subgraph: &BitSet) -> Vec<u8> {
    let graph = mol.graph();
    let order = graph.node_count();

    // An adjacency list of the current state of the algorithm. At each stage,
    // leaf vertices are pruned and removed from the adjacency list.
    let mut adjacencies = vec![BitSet::with_capacity(order); order];

    // Sets of canonical labels of the current node and all its children. A
    // canonical set is collected and moved to the parent set when its
    // corresponding node is refined into a leaf.
    let mut partial_canonical_sets = vec![Vec::<Vec<u8>>::new(); order];

    // Nodes not yet given a complete label and pruned.
    let mut unlabeled_vertices = BitSet::with_capacity(order);

    // Initialize each partial canonical set as a singleton containing the
    // node's own label, and set up adjacency list.
    for ix in subgraph.iter() {
        let (u, v) = graph
            .edge_endpoints(EdgeIndex::new(ix))
            .expect("malformed bitset!");

        for node in [u, v] {
            let index = node.index();
            if unlabeled_vertices.contains(index) {
                continue;
            }
            unlabeled_vertices.insert(index);
            let weight = graph.node_weight(node).unwrap();
            partial_canonical_sets[index].push(vec![weight.element().repr()]);
        }

        let (u, v) = (u.index(), v.index());
        adjacencies[u].insert(v);
        adjacencies[v].insert(u);
    }

    // Leaf pruning proceeds until the tree is an isolated edge or node.
    while unlabeled_vertices.len() > 2 {
        let leaves = unlabeled_vertices
            .iter()
            .filter(|&i| adjacencies[i].len() == 1)
            .collect::<Vec<_>>();

        for leaf in leaves {
            let parent = adjacencies[leaf].iter().next().unwrap();
            let edge = graph
                .edges_connecting(NodeIndex::new(parent), NodeIndex::new(leaf))
                .next()
                .unwrap();

            // Collapse parent-leaf edge + partial canonical set into canonical
            // label vector. Then move canonical label into parent's partial
            // canonical set.
            let mut canonical_label = vec![(*edge.weight()).repr()];
            canonical_label.extend(collapse_set(mem::take(&mut partial_canonical_sets[leaf])));
            partial_canonical_sets[parent].push(canonical_label);

            // Proceed as though the pruned leaf no longer exists.
            adjacencies[leaf].clear();
            adjacencies[parent].remove(leaf);
            unlabeled_vertices.remove(leaf);
        }
    }

    if unlabeled_vertices.len() == 2 {
        // Case 1: Tree collapses into isolated edge. Obtain node labels,
        // lexicographically sort, and then glue together into canonical label
        // alongside edge weight.
        let mut iter = unlabeled_vertices.iter();
        let (u, v) = (iter.next().unwrap(), iter.next().unwrap());
        let edge = graph
            .edges_connecting(NodeIndex::new(u), NodeIndex::new(v))
            .next()
            .unwrap();

        let u = collapse_set(mem::take(&mut partial_canonical_sets[u]));
        let v = collapse_set(mem::take(&mut partial_canonical_sets[v]));

        let (first, second) = if u < v { (u, v) } else { (v, u) };

        [(*edge.weight()).repr()]
            .into_iter()
            .chain(wrap_with_delimiters(first))
            .chain(wrap_with_delimiters(second))
            .collect()
    } else {
        // Case 2: Tree collapses into isolated node. Use node label.
        let canonical_root = unlabeled_vertices.iter().next().unwrap();
        let canonical_set = mem::take(&mut partial_canonical_sets[canonical_root]);
        collapse_set(canonical_set)
    }
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
