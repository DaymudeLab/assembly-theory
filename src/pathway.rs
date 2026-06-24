//! Reconstruct minimum assembly pathways from recursive search information.

use std::fmt;

use bit_set::BitSet;
use petgraph::{
    dot::{Config, Dot},
    graph::{DiGraph, EdgeIndex, NodeIndex},
};

use crate::{matches::Matches, molecule::Molecule, utils::subgraph_from_edge_mask};

/// Assembly pathway of pairwise joining operations yielding the target graph.
pub struct Pathway {
    /// A directed, acyclic graph representation of the assembly pathway; see
    /// [`Pathway::dag`] for details.
    dag: DiGraph<String, ()>,
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
        let mut dag = DiGraph::<String, ()>::new();
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

    /// Get the pathway's directed, acyclic graph representation.
    ///
    /// Nodes are fragments (represented as DOT strings) and edges `(u1, v)`
    /// and `(u2, v)` exist if fragments `u1` and `u2` are joined to produce
    /// fragment `v`. If `u` is duplicated and joined to itself to produce `v`,
    /// then the edge `(u, v)` is included twice.
    pub fn dag(&self) -> &DiGraph<String, ()> {
        &self.dag
    }
}

impl fmt::Display for Pathway {
    /// Format the assembly pathway DAG as a DOT string whose node labels are
    /// the DOT strings of their respective fragments.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:?}",
            Dot::with_config(self.dag(), &[Config::EdgeNoLabel])
        )
    }
}

/// Helper function for [`Pathway::new`] that takes as input an edge-disjoint
/// set of pieces and repeatedly joins all compatible pairs of pieces together,
/// tracking these join operations in the given assembly pathway.
fn join_pieces(mol: &Molecule, pieces: &mut Vec<BitSet>, pathway: &mut DiGraph<String, ()>) {
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
                // Get the DOT string representations of pieces i and j.
                let dot_i = format!(
                    "{:?}",
                    Dot::new(&subgraph_from_edge_mask(mol.graph(), &pieces[i]))
                );
                let dot_j = format!(
                    "{:?}",
                    Dot::new(&subgraph_from_edge_mask(mol.graph(), &pieces[j]))
                );

                // Ensure pieces i and j exist as nodes in the pathway DAG.
                let i_ix = match pathway.raw_nodes().iter().position(|v| v.weight == dot_i) {
                    Some(v_ix) => NodeIndex::new(v_ix),
                    _ => pathway.add_node(dot_i),
                };
                let j_ix = match pathway.raw_nodes().iter().position(|v| v.weight == dot_j) {
                    Some(v_ix) => NodeIndex::new(v_ix),
                    _ => pathway.add_node(dot_j),
                };

                // Join piece j into piece i and update the nodes contained in
                // the newly formed piece.
                let piece_j = pieces.swap_remove(j);
                pieces[i].union_with(&piece_j);
                let piece_nodes_j = piece_nodes.swap_remove(j);
                piece_nodes[i].union_with(&piece_nodes_j);

                // Get the DOT string representation of the new piece and add
                // it as a node to the pathway DAG.
                let dot_joined = format!(
                    "{:?}",
                    Dot::new(&subgraph_from_edge_mask(mol.graph(), &pieces[i]))
                );
                let joined_ix = pathway.add_node(dot_joined);

                // Add edges in the pathway DAG representing the join.
                pathway.add_edge(i_ix, joined_ix, ());
                pathway.add_edge(j_ix, joined_ix, ());

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
