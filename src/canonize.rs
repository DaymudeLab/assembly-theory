//! Create canonical labelings for molecular graphs.

use std::{collections::HashMap, cmp::Ordering, hash::Hash};

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
    /// Use the Nauty algorithm of McKay & Piperno (2014).
    Nauty,
    /// Use the algorithm of Faulon et al. (2004).
    Faulon,
    /// Use a fast tree canonization algorithm if applicable; else use Nauty.
    TreeNauty,
    /// Use a fast tree canonization algorithm if applicable; else use Faulon.
    TreeFaulon,
}

/// Canonical label returned by our graph canonization functions.
#[derive(Hash, PartialEq, Eq, Debug)]
pub enum Labeling {
    /// The label returned by the graph_canon crate.
    Nauty(CanonLabeling<AtomOrBond>),

    /// A string label returned by our implementation of Faulon et al. (2004).
    // TODO: This should be a `Vec<u8>`
    Faulon(String),
}

/// Obtain a canonical labeling of the specified `subgraph` using the
/// algorithm specified by `mode`.
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

// First checks if the subgraph forms a tree (given the assumption that subgraph is always connected)
// returns None if not a tree
// Otherwise returns canon string for the (labeled) tree
pub fn canonize_subtree(molecule: &Molecule, subgraph: &BitSet) -> Option<Vec<u8>> {
    let mgraph = molecule.graph();
    let mut vtx_set = BitSet::with_capacity(mgraph.node_count());
    let mut subgraph_adj = vec![BitSet::with_capacity(mgraph.node_count()); mgraph.node_count()];
    let mut vtx_strs: Vec<Vec<Vec<u8>>> = vec![vec![]; mgraph.node_count()];
    let mut vtx_label: Vec<Vec<u8>> = vec![vec![]; mgraph.node_count()];

    // count nodes in the subgraph tree
    for subgraph_bond_idx in subgraph {
        let bond_idx = EdgeIndex::new(subgraph_bond_idx);
        let (start_atom_idx, end_atom_idx) = mgraph.edge_endpoints(bond_idx).unwrap();

        let inserted_start = vtx_set.insert(start_atom_idx.index());
        let inserted_end = vtx_set.insert(end_atom_idx.index());

        // update adj matrix for subgraph
        subgraph_adj[start_atom_idx.index()].insert(end_atom_idx.index());
        subgraph_adj[end_atom_idx.index()].insert(start_atom_idx.index());

        // add empty strings for all the nodes
        if inserted_start {
            vtx_label[start_atom_idx.index()].push(b'(');
            vtx_label[start_atom_idx.index()].push(
                mgraph
                    .node_weight(start_atom_idx)
                    .unwrap()
                    .element()
                    .to_string()
                    .as_bytes()[0],
            );
            vtx_label[start_atom_idx.index()].push(b')');
        }
        if inserted_end {
            vtx_label[end_atom_idx.index()].push(b'(');
            vtx_label[end_atom_idx.index()].push(
                mgraph
                    .node_weight(end_atom_idx)
                    .unwrap()
                    .element()
                    .to_string()
                    .as_bytes()[0],
            );
            vtx_label[end_atom_idx.index()].push(b')');
        }
    }

    // check if the subgraph is a tree or not (given the assumption that subgraph is connected graph)
    // e = n - 1

    if subgraph.len() == (vtx_set.len() - 1) {
        // keep merging the labels until 1/2 vertices remain
        while vtx_set.len() > 2 {
            // locate all the leaves
            let leaves: Vec<usize> = subgraph_adj
                .iter()
                .enumerate()
                .filter_map(|(i, adj_list)| if adj_list.len() == 1 { Some(i) } else { None })
                .collect();
            let mut parents = BitSet::with_capacity(mgraph.node_count());
            // add strings of the leaves to the parent's string array
            leaves.iter().for_each(|leaf_id| {
                let parent = subgraph_adj[*leaf_id].iter().next().unwrap();
                let leaf_str = &vtx_label[*leaf_id];
                vtx_strs[parent].push(leaf_str.clone());

                // add parent to parents list
                parents.insert(parent);

                // remove leaf node from adj matrix
                subgraph_adj[*leaf_id].clear();
                subgraph_adj[parent].remove(*leaf_id);
                // remove leaf node from bit-set
                vtx_set.remove(*leaf_id);
            });

            // merge leaf strings in parent's primary string
            parents.iter().for_each(|parent| {
                // only do this once the parent has seen all its children and becomes a leaf in next iteration
                if subgraph_adj[parent].len() <= 1 {
                    vtx_strs[parent].sort();
                    let parent_str = vtx_strs[parent].join(&b',');

                    vtx_label[parent].pop(); // remove the closing bracket in starting parent label: "(X)"
                    vtx_label[parent].push(b',');
                    vtx_label[parent].extend_from_slice(&parent_str);
                    vtx_label[parent].push(b')');
                }
            });
        }

        let vtx_1 = vtx_set.iter().next().unwrap();
        let vtx_1_str = &vtx_label[vtx_1];

        // when 2 vertices are left, merge their labels
        if vtx_set.len() == 2 {
            let vtx_2 = vtx_set.iter().next().unwrap();
            let vtx_2_str = &vtx_label[vtx_2];
            let mut return_str_vec: Vec<u8> = vec![];
            if vtx_1_str.cmp(vtx_2_str) == Ordering::Less {
                return_str_vec.push(b'(');
                return_str_vec.extend(vtx_1_str);
                return_str_vec.push(b',');
                return_str_vec.extend(vtx_2_str);
                return_str_vec.push(b')');
                Some(return_str_vec)
            } else {
                return_str_vec.push(b'(');
                return_str_vec.extend(vtx_2_str);
                return_str_vec.push(b',');
                return_str_vec.extend(vtx_1_str);
                return_str_vec.push(b')');
                Some(return_str_vec)
            }
        } else {
            // Some(vtx_1_str.as_bytes().to_vec())
            Some(vtx_1_str.clone())
        }
    } else {
        None
    }
}

mod tests {
    #[allow(unused_imports)]
    use super::*;

    use std::fs;
    use std::path::PathBuf;

    #[allow(unused_imports)]
    use petgraph::algo::is_isomorphic_matching;

    use crate::loader;

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

    #[test]
    fn canonize_benzene() {
        let path = PathBuf::from(format!("./data/checks/benzene.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let mut subgraph = BitSet::from_iter(molecule.graph().edge_indices().map(|e| e.index()));
        subgraph.remove(molecule.graph().edge_indices().next().unwrap().index());
        let canonical_repr = canonize_subtree(&molecule, &subgraph).unwrap();

        assert_eq!(
            String::from_utf8(canonical_repr).unwrap(),
            "((C,(C,(C))),(C,(C,(C))))"
        )
    }

    #[test]
    fn canonize_anthracene() {
        let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let mut subgraph = BitSet::from_iter(molecule.graph().edge_indices().map(|e| e.index()));
        subgraph.remove(0);
        subgraph.remove(6);
        subgraph.remove(9);
        let canonical_repr = canonize_subtree(&molecule, &subgraph).unwrap();

        assert_eq!(
            String::from_utf8(canonical_repr).unwrap(),
            "((C,(C,(C),(C,(C,(C,(C)))))),(C,(C,(C),(C,(C,(C,(C)))))))"
        )
    }

    #[test]
    fn canonize_aspirin() {
        let path = PathBuf::from(format!("./data/checks/aspirin.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let mut subgraph = BitSet::from_iter(molecule.graph().edge_indices().map(|e| e.index()));
        subgraph.remove(molecule.graph().edge_indices().next().unwrap().index());
        let canonical_repr = canonize_subtree(&molecule, &subgraph).unwrap();

        assert_eq!(
            String::from_utf8(canonical_repr).unwrap(),
            "(C,(C,(C,(C,(C,(O),(O))))),(C,(C,(O,(C,(C),(O))))))"
        )
    }
}
