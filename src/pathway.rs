//! Assembly pathway reconstruction from top-down decomposition results.
//!
//! Implements the `Generate_Pathway` algorithm from the supplementary
//! information of [Seet et al. (2025)](https://doi.org/10.1021/acs.jcim.5c01964).

use std::fmt;

use bit_set::BitSet;
use petgraph::graph::EdgeIndex;

use crate::{matches::Matches, molecule::Molecule};

/// Format an edge set as a list of `(src, dst)` vertex pairs.
pub fn edges_as_vertex_pairs(mol: &Molecule, edges: &BitSet) -> Vec<(usize, usize)> {
    edges
        .iter()
        .map(|e| {
            let (src, dst) = mol.graph().edge_endpoints(EdgeIndex::new(e)).unwrap();
            (src.index(), dst.index())
        })
        .collect()
}

/// A single step in an assembly pathway: joining two pieces into one.
#[derive(Debug, Clone)]
pub struct PathwayStep {
    /// First piece being joined (edge set).
    pub piece_a: BitSet,
    /// Second piece being joined (edge set).
    pub piece_b: BitSet,
    /// The result of joining the two pieces (union of edge sets).
    pub result: BitSet,
}

impl fmt::Display for PathwayStep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let a: Vec<usize> = self.piece_a.iter().collect();
        let b: Vec<usize> = self.piece_b.iter().collect();
        let r: Vec<usize> = self.result.iter().collect();
        write!(f, "{:?} + {:?} -> {:?}", a, b, r)
    }
}

/// Check whether two edge sets share at least one node (vertex) in the
/// molecular graph.
fn shares_node(mol: &Molecule, a: &BitSet, b: &BitSet) -> bool {
    let mut nodes_a = BitSet::with_capacity(mol.graph().node_count());
    for e in a {
        let (src, dst) = mol.graph().edge_endpoints(EdgeIndex::new(e)).unwrap();
        nodes_a.insert(src.index());
        nodes_a.insert(dst.index());
    }
    for e in b {
        let (src, dst) = mol.graph().edge_endpoints(EdgeIndex::new(e)).unwrap();
        if nodes_a.contains(src.index()) || nodes_a.contains(dst.index()) {
            return true;
        }
    }
    false
}

/// `Consistent_Join`: merge any two pieces sharing a vertex, recording each
/// join as a pathway step. Returns `true` if any join was performed.
fn consistent_join(mol: &Molecule, pieces: &mut Vec<BitSet>, steps: &mut Vec<PathwayStep>) -> bool {
    let mut joined = false;
    let mut i = 0;
    while i < pieces.len() {
        let mut j = i + 1;
        let mut did_join = false;
        while j < pieces.len() {
            if shares_node(mol, &pieces[i], &pieces[j]) {
                // Join pieces[j] into pieces[i].
                let piece_a = pieces[i].clone();
                let piece_b = pieces.swap_remove(j);
                pieces[i].union_with(&piece_b);
                steps.push(PathwayStep {
                    piece_a,
                    piece_b,
                    result: pieces[i].clone(),
                });
                joined = true;
                did_join = true;
                // Don't increment j — a new element may have been swapped in.
            } else {
                j += 1;
            }
        }
        if !did_join {
            i += 1;
        }
        // If we joined, re-check pieces[i] against remaining pieces.
    }
    joined
}

/// Find pieces that overlap with a given fragment (by shared nodes).
fn filter_pieces(mol: &Molecule, pieces: &[BitSet], fragment: &BitSet) -> Vec<usize> {
    pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| shares_node(mol, p, fragment))
        .map(|(i, _)| i)
        .collect()
}

/// `Duplicate_Construction`: process duplicates (matched fragment pairs) from
/// smallest to largest, adding copies or first reassembling containing pieces.
fn duplicate_construction(
    mol: &Molecule,
    duplicates: &[(BitSet, BitSet)],
    pieces: &mut Vec<BitSet>,
    steps: &mut Vec<PathwayStep>,
) {
    for (h1, h2) in duplicates {
        // Check if h2 (the "template") is already exactly one of the pieces.
        let exact_match = pieces.iter().position(|p| p == h2);
        if let Some(_) = exact_match {
            // The template is already a piece; just add the copy.
            pieces.push(h1.clone());
        } else {
            // Filter pieces that overlap with h2.
            let filtered_ixs = filter_pieces(mol, pieces, h2);
            if filtered_ixs.is_empty() {
                // Template edges are already part of a larger assembled piece;
                // just add the copy.
                pieces.push(h1.clone());
                continue;
            }

            // Extract filtered pieces and remove them from the main list.
            // Remove in reverse index order to avoid invalidating indices.
            let mut filtered: Vec<BitSet> = Vec::new();
            let mut sorted_ixs = filtered_ixs;
            sorted_ixs.sort_unstable_by(|a, b| b.cmp(a));
            for ix in &sorted_ixs {
                filtered.push(pieces.swap_remove(*ix));
            }

            // Join the filtered pieces until they form one piece.
            while filtered.len() > 1 {
                if !consistent_join(mol, &mut filtered, steps) {
                    break;
                }
            }

            // Add both the reassembled piece and the duplicate copy.
            pieces.extend(filtered);
            pieces.push(h1.clone());
        }
    }
}

/// Reconstruct the bottom-up assembly pathway from the top-down decomposition
/// results.
///
/// # Arguments
/// - `mol`: The molecule.
/// - `matches`: The matches structure used during the search.
/// - `match_sequence`: The sequence of global `match_ix` values from root to
///   leaf of the optimal decomposition path.
/// - `remnants`: The leaf state's fragments (the pieces left after all
///   duplicates have been removed). Note: these may exclude singleton (single-
///   edge) fragments that were dropped during the search.
///
/// # Returns
/// A vector of [`PathwayStep`]s representing the assembly pathway from basic
/// fragments to the complete molecule.
pub fn generate_pathway(
    mol: &Molecule,
    matches: &Matches,
    match_sequence: &[usize],
    _remnants: &[BitSet],
) -> Vec<PathwayStep> {
    let mut steps: Vec<PathwayStep> = Vec::new();

    // Build the list of duplicates from the match sequence. Each duplicate is
    // (copy, template): the search kept h1 (template) and discarded h2 (copy).
    // For bottom-up reconstruction, we reverse so smallest fragments come first.
    let mut duplicates: Vec<(BitSet, BitSet)> = match_sequence
        .iter()
        .map(|&mix| {
            let (h1, h2) = matches.match_fragments(mix);
            (h2.clone(), h1.clone()) // (copy, template)
        })
        .collect();
    duplicates.reverse();

    // Compute the edges present at the leaf level: all molecule edges minus
    // the net-removed edges (h2 from each decomposition step).
    let mut present_edges = BitSet::with_capacity(mol.graph().edge_count());
    present_edges.extend(mol.graph().edge_indices().map(|ix| ix.index()));
    for &mix in match_sequence.iter() {
        let (_, h2) = matches.match_fragments(mix);
        present_edges.difference_with(h2);
    }

    // Start pieces as individual edges (basic units).
    let mut pieces: Vec<BitSet> = present_edges
        .iter()
        .map(|e| {
            let mut bs = BitSet::with_capacity(mol.graph().edge_count());
            bs.insert(e);
            bs
        })
        .collect();

    // Process duplicates: for each one, reassemble its template from basic
    // pieces if needed, then add the copy.
    duplicate_construction(mol, &duplicates, &mut pieces, &mut steps);

    // Iteratively join remaining pieces that share vertices.
    loop {
        if !consistent_join(mol, &mut pieces, &mut steps) {
            break;
        }
    }

    steps
}
