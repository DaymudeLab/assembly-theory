//! Enumerate connected, non-induced subgraphs of a molecular graph.

use std::collections::HashSet;

use bit_set::BitSet;
use clap::ValueEnum;
use petgraph::graph::EdgeIndex;

use crate::molecule::Molecule;

/// Strategy for enumerating connected, non-induced subgraphs.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum EnumerateMode {
    /// Grow connected subgraphs from each edge using BFS.
    Bfs,
    /// Like Bfs, but at each level of the BFS, prune any subgraphs that do not
    /// have isomorphic components since these will not be useful later.
    BfsPrune,
    /// From a subgraph, choose an edge from its boundary and recursively grow
    /// it by adding this edge or erode its remainder by discarding the edge.
    GrowErode,
    /// An iterative (memory-efficient) implementation of GrowErode.
    GrowErodeIterative,
}

/// Return an interator over all connected, non-induced subgraphs of the
/// molecular graph `mol` using the algorithm specified by `mode`.
pub fn enumerate_subgraphs(mol: &Molecule, mode: EnumerateMode) -> impl Iterator<Item = BitSet> {
    match mode {
        EnumerateMode::GrowErode => grow_erode(mol).into_iter(),
        _ => {
            panic!("The chosen --enumerate mode is not implemented yet!");
        }
    }
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges; at
/// each recursive step, choose one edge from the current subgraph's boundary
/// and either add it to the subgraph or discard it from the remainder.
fn grow_erode(mol: &Molecule) -> HashSet<BitSet> {
    // Initialize the current subset of edges and its union with its edge
    // boundary (in this algorithm, "subsetplus") as empty.
    let subset = BitSet::new();
    let subsetplus = BitSet::new();

    // The remainder is all edges not in the current subset; initially, this is
    // everything.
    let remainder = BitSet::from_iter(mol.graph().edge_indices().map(|ix| ix.index()));

    // Set up a set of subgraphs enumerated so far.
    let mut subgraphs = HashSet::new();

    // Recurse, and ultimately return the final set of enumerated subgraphs.
    grow_erode_recurse(mol, subset, subsetplus, remainder, &mut subgraphs);
    subgraphs
}

fn grow_erode_recurse(
    mol: &Molecule,
    mut subset: BitSet,
    mut subsetplus: BitSet,
    mut remainder: BitSet,
    subgraphs: &mut HashSet<BitSet>,
) {
    // Get the next edge from the current subset's boundary or, if the subset
    // is empty, from the remainder.
    let candidate = if subset.is_empty() {
        remainder.iter().next()
    } else {
        remainder.intersection(&subsetplus).next()
    };

    if let Some(e) = candidate {
        // In the first recursive branch, discard the candidate edge entirely.
        remainder.remove(e);
        grow_erode_recurse(
            mol,
            subset.clone(),
            subsetplus.clone(),
            remainder.clone(),
            subgraphs,
        );

        // The other recursive branch will add the candidate edge to the
        // current subset and update the boundary accordingly. However, since
        // we ultimately only care about connected, non-induced subgraphs that
        // may be part of a non-overlapping isomorphic pair, we need not
        // recurse if adding the edge exceeds |E|/2 edges.
        if subset.len() < mol.graph().edge_count() / 2 {
            // Add the candidate edge to the current subset.
            subset.insert(e);

            // Grow the boundary to include edges incident to the candidate.
            let (src, dst) = mol
                .graph()
                .edge_endpoints(EdgeIndex::new(e))
                .expect("malformed input");
            subsetplus.extend(
                mol.graph()
                    .neighbors(src)
                    .filter_map(|n| mol.graph().find_edge(src, n).map(|ix| ix.index())),
            );
            subsetplus.extend(
                mol.graph()
                    .neighbors(dst)
                    .filter_map(|n| mol.graph().find_edge(dst, n).map(|ix| ix.index())),
            );

            // Recurse.
            grow_erode_recurse(mol, subset, subsetplus, remainder, subgraphs);
        }
    } else if subset.len() > 1 {
        // When all candidate edges have been exhausted, add this subset as a
        // new subgraph if it is nonempty.
        subgraphs.insert(subset);
    }
}
