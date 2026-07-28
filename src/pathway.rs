//! Reconstruct minimum assembly pathways from recursive search information.

use std::fmt;

use bit_set::BitSet;
use petgraph::{
    dot::{Config, Dot},
    graph::{DiGraph, EdgeIndex, NodeIndex},
};

use crate::{
    canonize::{canonize, CanonizeMode},
    matches::Matches,
    molecule::Molecule,
};

/// Assembly pathway of pairwise joining operations yielding the target graph.
pub struct Pathway {
    /// A directed, acyclic graph representation of the assembly pathway; see
    /// [`Pathway::dag`] for details.
    dag: DiGraph<BitSet, BitSet>,
}

impl Pathway {
    /// Construct the assembly [`Pathway`] of the given molecule corresponding
    /// to the given order of removed matches. Efficiently implements the
    /// duplicate/remnant algorithm described in the Supporting Information
    /// of [Seet et al. (2025)](https://doi.org/10.1021/acs.jcim.5c01964).
    pub fn new(mol: &Molecule, matches: &Matches, removal_order: &Vec<usize>) -> Self {
        // Create a fragment representing the complete molecule.
        let mut mol_frag = BitSet::new();
        mol_frag.extend(mol.graph().edge_indices().map(|ix| ix.index()));

        // Get all duplicate fragments from removed matches in reverse order so
        // that smaller fragments can be used to construct larger ones.
        let duplicates: Vec<(BitSet, BitSet)> = removal_order
            .iter()
            .rev()
            .map(|&match_ix| {
                let (h1, h2) = matches.match_fragments(match_ix);
                (h1.clone(), h2.clone())
            })
            .collect();

        // Initialize the assembly pathway's bag of "pieces" as all remnants,
        // i.e., singleton edges comprising the fragments remaining after all
        // duplicate fragments are removed. As in [`assembly::fragments`], we
        // retain each match's first duplicate fragment and remove the second.
        let mut remnants = mol_frag.clone();
        for (_, h2) in &duplicates {
            remnants.difference_with(h2);
        }
        let mut pieces: Vec<BitSet> = remnants
            .iter()
            .map(|edge_ix| {
                let mut piece = BitSet::with_capacity(mol.graph().edge_count());
                piece.insert(edge_ix);
                piece
            })
            .collect();

        // Use the bag of pieces to construct all duplicate fragments across
        // matches, starting with the smallest matches and building up.
        let mut dag = DiGraph::<BitSet, BitSet>::new();
        for (h1, h2) in &duplicates {
            // If h1 is not already in the bag (as it could be if it is used in
            // multiple matches), construct it by joining pieces in the bag.
            if pieces.iter().find(|&p| p == h1) == None {
                // Take all pieces that are part of h1; this algorithm ensures
                // that such pieces are an edge-disjoint partition of h1.
                let mut h1_pieces: Vec<BitSet> =
                    pieces.extract_if(.., |p| p.is_subset(h1)).collect();

                // Join those pieces iteratively into h1 and add it to the bag.
                join_pieces(mol, &mut h1_pieces, &mut dag);
                pieces.push(h1.clone());
            }

            // Having constructed/found h1, add its duplicate h2 to the bag.
            pieces.push(h2.clone());
        }

        // Join whatever fragments are left in the bag until the original
        // molecular graph is recovered.
        join_pieces(mol, &mut pieces, &mut dag);

        Self { dag }
    }

    /// Get the pathway's directed, acyclic (multi)graph representation whose
    /// nodes represent fragments and edges represent joining operations.
    ///
    /// In detail, each node is labeled by a BitSet indexing the bonds its
    /// fragment contains. An edge `(u1, v)` labeled `l1` and an edge `(u2, v)`
    /// labeled `l2` exist iff the fragment BitSet `l1` is isomorphic to `u1`,
    /// the fragment BitSet `l2` is isomorphic to `u2`, and `l1` and `l2` are
    /// joined to produce a fragment isomorphic to `v`.
    ///
    /// Note that an edge `(u, v)` may be in the DAG multiple times, e.g., if
    /// two fragments isomorphic to `u` are joined to produce a fragment
    /// isomorphic to `v`. It may even be in the DAG multiple times *with the
    /// same label*, e.g., if some fragment isomorphic to `u` is duplicated and
    /// joined to itself to produce a fragment isomorphic to `v`.
    pub fn dag(&self) -> &DiGraph<BitSet, BitSet> {
        &self.dag
    }
}

impl fmt::Display for Pathway {
    /// Format the assembly pathway DAG as a DOT string whose node and edge
    /// labels are the BitSets of their respective fragments' bonds.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:?}",
            Dot::with_attr_getters(
                self.dag(),
                &[Config::EdgeNoLabel, Config::NodeNoLabel], // Suppress default labels.
                &|_, e| format!("label = \"{}\" ", e.weight()),
                &|_, (_, n_weight)| format!("label = \"{}\" ", n_weight)
            )
        )
    }
}

/// Helper function for [`Pathway::new`] that takes as input an edge-disjoint
/// set of pieces and repeatedly joins all compatible pairs of pieces together,
/// tracking these join operations in the given assembly pathway.
fn join_pieces(mol: &Molecule, pieces: &mut Vec<BitSet>, pathway: &mut DiGraph<BitSet, BitSet>) {
    // For each piece, identify the set of nodes it contains.
    let mut piece_nodes = Vec::<BitSet>::new();
    for piece in pieces.iter() {
        let mut nodes = BitSet::with_capacity(mol.graph().node_count());
        for e_ix in piece.iter() {
            let (end1, end2) = mol.graph().edge_endpoints(EdgeIndex::new(e_ix)).unwrap();
            nodes.insert(end1.index());
            nodes.insert(end2.index());
        }
        piece_nodes.push(nodes);
    }

    // Join all compatible pieces together.
    let mut i = 0;
    while i < pieces.len() {
        let mut j = i + 1;
        let mut joined_piece_i = false;
        while j < pieces.len() {
            // If pieces i and j share a node, join them into one piece.
            if !piece_nodes[i].is_disjoint(&piece_nodes[j]) {
                // Join piece j into piece i and update the nodes contained in
                // the newly formed piece.
                let piece_i = pieces[i].clone();
                let piece_j = pieces.swap_remove(j);
                pieces[i].union_with(&piece_j);
                let piece_nodes_j = piece_nodes.swap_remove(j);
                piece_nodes[i].union_with(&piece_nodes_j);

                // Ensure pieces i and j exist as nodes in the pathway DAG.
                let i_ix = find_or_add_dag_node(mol, &piece_i, pathway);
                let j_ix = find_or_add_dag_node(mol, &piece_j, pathway);

                // Add the joined piece and its incoming edges to the DAG.
                let joined_ix = pathway.add_node(pieces[i].clone());
                pathway.add_edge(i_ix, joined_ix, piece_i);
                pathway.add_edge(j_ix, joined_ix, piece_j);

                joined_piece_i = true;
            } else {
                // Piece j is not compatible with piece i; check the next one.
                j += 1
            }
        }

        // If piece i was not joined to any later pieces, it is a "final"
        // structure in this assembly (e.g., one connected component in a
        // disconnected input molecular graph) and can now be set aside.
        if !joined_piece_i {
            i += 1
        }
    }
}

// Helper function for [`join_pieces`] that finds the pathway node isomorphic
// to the given piece or adds a new node if one doesn't exist.
fn find_or_add_dag_node(
    mol: &Molecule,
    piece: &BitSet,
    pathway: &mut DiGraph<BitSet, BitSet>,
) -> NodeIndex {
    let piece_label = canonize(mol, piece, CanonizeMode::TreeNauty);
    match pathway
        .raw_nodes()
        .iter()
        .position(|v| canonize(mol, &v.weight, CanonizeMode::TreeNauty) == piece_label)
    {
        Some(v_ix) => NodeIndex::new(v_ix),
        _ => pathway.add_node(piece.clone()),
    }
}
