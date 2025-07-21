//! Enumerate connected, non-induced subgraphs of a molecular graph.
//!
//! Specifically, for a molecule with |E| edges, enumerate all connected,
//! non-induced subgraphs with at most |E|/2 edges. Any larger subgraphs
//! cannot be "duplicatable" (i.e., in a pair of non-overlapping, isomorphic
//! subgraphs), so we don't need them.

use std::collections::HashSet;

use bit_set::BitSet;
use clap::ValueEnum;
use petgraph::graph::EdgeIndex;

use crate::{molecule::Molecule, utils::edge_neighbors};

/// Strategy for enumerating connected, non-induced subgraphs.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum EnumerateMode {
    /// Grow connected subgraphs from each edge using iterative extension.
    Extend,
    /// From a subgraph, choose an edge from its boundary and either grow it by
    /// adding this edge or erode its remainder by discarding the edge.
    GrowErode,
}

/// Return an interator over all connected, non-induced subgraphs of the
/// molecular graph `mol` using the algorithm specified by `mode`.
pub fn enumerate_subgraphs(mol: &Molecule, mode: EnumerateMode) -> HashSet<BitSet> {
    match mode {
        EnumerateMode::Extend => extend(mol),
        EnumerateMode::GrowErode => grow_erode(mol),
    }
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges using
/// a process of iterative extension starting from each individual edge.
fn extend(mol: &Molecule) -> HashSet<BitSet> {
    // Maintain a vector of sets of subgraphs at each level of the process,
    // starting with all edges individually at the first level.
    let mut subgraphs: Vec<HashSet<BitSet>> =
        vec![HashSet::from_iter(mol.graph().edge_indices().map(|ix| {
            let mut subgraph = BitSet::new();
            subgraph.insert(ix.index());
            subgraph
        }))];

    // At each level, collect and deduplicate all ways of extending subgraphs
    // by one neighboring edge.
    for level in 0..(mol.graph().edge_count() / 2) {
        let mut extended_subgraphs = HashSet::new();
        for subgraph in &subgraphs[level] {
            // Find all "frontier" edges incident to this subgraph (contains
            // both this subgraph's edges and its edge boundary).
            let frontier =
                BitSet::from_iter(subgraph.iter().flat_map(|i| {
                    edge_neighbors(mol.graph(), EdgeIndex::new(i)).map(|ix| ix.index())
                }));

            // Collect and deduplicate all subgraphs obtained by extending the
            // current subgraph using one edge from its boundary.
            for edge in frontier.difference(subgraph) {
                let mut extended_subgraph = subgraph.clone();
                extended_subgraph.insert(edge);
                extended_subgraphs.insert(extended_subgraph);
            }
        }

        subgraphs.push(extended_subgraphs);
    }

    // Return an iterator over subgraphs, skipping singleton edges.
    subgraphs
        .into_iter()
        .skip(1)
        .flatten()
        .collect::<HashSet<_>>()
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges; at
/// each step, choose one edge from the current subgraph's boundary and either
/// add it to the subgraph or discard it from the remainder. See:
/// - https://stackoverflow.com/a/15722579
/// - https://stackoverflow.com/a/15658245
fn grow_erode(mol: &Molecule) -> HashSet<BitSet> {
    // Initialize the current subgraph and its "frontier" (i.e., the union of
    // its edges and its edge boundary) as well as the "remainder", which is
    // all edges not in the current subgraph.
    let subgraph = BitSet::new();
    let frontier = BitSet::new();
    let remainder = BitSet::from_iter(mol.graph().edge_indices().map(|ix| ix.index()));

    // Set up a set of subgraphs enumerated so far.
    let mut subgraphs = HashSet::new();

    // Maintain a stack of subgraph instances to extend.
    let mut stack = vec![(subgraph, frontier, remainder)];
    while let Some((mut subgraph, mut frontier, mut remainder)) = stack.pop() {
        // Get the next edge from the subgraph's edge boundary or, if the
        // subgraph is empty, from the remainder.
        let candidate = if subgraph.is_empty() {
            remainder.iter().next()
        } else {
            remainder.intersection(&frontier).next()
        };

        if let Some(e) = candidate {
            // Make a new instance by discarding the candidate edge entirely.
            remainder.remove(e);
            stack.push((subgraph.clone(), frontier.clone(), remainder.clone()));

            // Make another instance by adding the candidate edge to the
            // subgraph and updating the frontier accordingly if the new
            // subgraph was not previously enumerated and is not too large to
            // be part of a non-overlapping isomorphic pair.
            subgraph.insert(e);
            if !subgraphs.contains(&subgraph) && subgraph.len() <= mol.graph().edge_count() / 2 {
                frontier
                    .extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|ix| ix.index()));
                stack.push((subgraph, frontier, remainder));
            }
        } else if subgraph.len() > 1 {
            // When all candidate edges are exhausted, collect this subgraph
            // unless it is just a singleton edge (basic unit).
            subgraphs.insert(subgraph);
        }
    }

    subgraphs
}
