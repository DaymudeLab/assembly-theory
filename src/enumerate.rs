//! Enumerate connected, non-induced subgraphs of a molecular graph.

use std::collections::{HashMap, HashSet};

use bit_set::BitSet;
use clap::ValueEnum;
use graph_canon::CanonLabeling;
use petgraph::graph::EdgeIndex;

use crate::{
    canonize::subgraph_to_cgraph,
    molecule::{AtomOrBond, Molecule},
    utils::edge_neighbors,
};

/// Strategy for enumerating connected, non-induced subgraphs.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
pub enum EnumerateMode {
    /// Grow connected subgraphs from each edge using iterative extension.
    Extend,
    /// Like Extend, but at each level, prune any subgraphs that do not have
    /// isomorphic components since these will not be useful later.
    ExtendPrune,
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
        EnumerateMode::Extend => extend_iterative(mol).into_iter(),
        EnumerateMode::GrowErode => grow_erode(mol).into_iter(),
        _ => {
            panic!("The chosen --enumerate mode is not implemented yet!");
        }
    }
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges using
/// a process of iterative extension starting from each individual edge.
fn extend_iterative(mol: &Molecule) -> HashSet<BitSet> {
    // Maintain a vector of sets of subgraphs at each level of the process,
    // starting with all edges individually at the first level.
    let mut subgraphs: Vec<HashSet<BitSet>> =
        vec![HashSet::from_iter(mol.graph().edge_indices().map(|ix| {
            let mut set = BitSet::new();
            set.insert(ix.index());
            set
        }))];

    // At each level, collect and deduplicate all ways of extending subgraphs
    // by one neighboring edge.
    for level in 0..(mol.graph().edge_count() / 2) {
        let mut extended_subgraphs = HashSet::new();
        for subgraph in &subgraphs[level] {
            // Find all "frontier" edges incident to this subgraph (this
            // contains both this subgraph's edges and its edge boundary).
            let frontier = BitSet::from_iter(subgraph
                .iter()
                .map(|i| edge_neighbors(&mol.graph(), EdgeIndex::new(i)).map(|ix| ix.index()))
                .flatten(),
            );
            
            // Collect and deduplicate all subgraphs obtained by extending the
            // current subgraph using one edge from its boundary.
            for edge in frontier.difference(&subgraph) {
                let mut extended_subgraph = subgraph.clone();
                extended_subgraph.insert(edge);
                extended_subgraphs.insert(extended_subgraph);
            }
        }

        subgraphs.push(extended_subgraphs);
    }

    // Return an iterator over subgraphs, skipping singleton edges.
    subgraphs.into_iter().skip(1).flatten().collect::<HashSet<_>>()
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges; at
/// each recursive step, choose one edge from the current subgraph's boundary
/// and either add it to the subgraph or discard it from the remainder. See:
/// - https://stackoverflow.com/a/15722579
/// - https://stackoverflow.com/a/15658245
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
    recurse_grow_erode(mol, subset, subsetplus, remainder, &mut subgraphs);
    subgraphs
}

fn recurse_grow_erode(
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
        recurse_grow_erode(
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
            subsetplus.extend(edge_neighbors(&mol.graph(), EdgeIndex::new(e)).map(|ix| ix.index()));

            // Recurse.
            recurse_grow_erode(mol, subset, subsetplus, remainder, subgraphs);
        }
    } else if subset.len() > 1 {
        // When all candidate edges have been exhausted, add this subset as a
        // new subgraph if it is more than just a singleton edge (basic unit).
        subgraphs.insert(subset);
    }
}

// Iterative version of grow_erode_recurse
fn grow_erode_iterative(mol: &Molecule) -> impl Iterator<Item = BitSet> {
    let remainder = BitSet::from_iter(mol.graph().edge_indices().map(|ix| ix.index()));
    let subset = BitSet::new();
    let neighbors = BitSet::new();
    let mut stack = vec![(subset, neighbors, remainder)];

    let mut solutions = HashSet::new();

    while let Some((mut subset, mut neighbors, mut remainder)) = stack.pop() {
        let candidate = if subset.is_empty() {
            remainder.iter().next()
        } else {
            remainder.intersection(&neighbors).next()
        };

        if let Some(e) = candidate {
            remainder.remove(e);
            stack.push((subset.clone(), neighbors.clone(), remainder.clone()));

            subset.insert(e);
            if solutions.contains(&subset) || subset.len() > mol.graph().edge_count() / 2 {
                continue;
            }

            neighbors.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|ix| ix.index()));
            stack.push((subset, neighbors, remainder));
        } else if subset.len() > 1 {
            solutions.insert(subset);
        }
    }
    solutions.into_iter()
}

// TODO: matches_by_iterative_expansion and naive_matches_by_iterative_expansion are big, ugly
// functions that need to be split up and cleaned.
pub fn matches_by_iterative_expansion(mol: &Molecule) -> impl Iterator<Item = (BitSet, BitSet)> {
    let mut solutions: HashMap<BitSet, BitSet> =
        HashMap::from_iter(mol.graph().edge_indices().map(|ix| {
            let mut set = BitSet::new();
            set.insert(ix.index());
            let neighborhood =
                BitSet::from_iter(edge_neighbors(&mol.graph(), ix).map(|e| e.index()));
            (set, neighborhood)
        }));

    let mut matches = Vec::new();

    for _ in 0..(mol.graph().edge_count() / 2) {
        let mut next_set = HashMap::new();
        for (subgraph, neighborhood) in solutions {
            for neighbor in neighborhood.difference(&subgraph) {
                if subgraph.contains(neighbor) {
                    continue;
                }

                let mut next = subgraph.clone();
                next.insert(neighbor);

                if next_set.contains_key(&next) {
                    continue;
                }

                let mut next_neighborhood = neighborhood.clone();
                next_neighborhood.extend(
                    edge_neighbors(&mol.graph(), EdgeIndex::new(neighbor)).map(|e| e.index()),
                );

                next_set.insert(next, next_neighborhood);
            }
        }

        let mut local_isomorphic_map = HashMap::<CanonLabeling<AtomOrBond>, Vec<BitSet>>::new();
        for (subgraph, _) in &next_set {
            let cgraph = subgraph_to_cgraph(mol, &subgraph);
            let repr = CanonLabeling::new(&cgraph);

            local_isomorphic_map
                .entry(repr)
                .and_modify(|bucket| bucket.push(subgraph.clone()))
                .or_insert(vec![subgraph.clone()]);
        }

        for (_, sets) in &local_isomorphic_map {
            if sets.len() == 1 {
                next_set.remove(&sets[0]);
            }
        }

        solutions = next_set;
        for bucket in local_isomorphic_map.values().filter(|v| v.len() > 1) {
            for (i, first) in bucket.iter().enumerate() {
                for second in &bucket[i..] {
                    if first.is_disjoint(second) {
                        matches.push((first.clone(), second.clone()));
                    }
                }
            }
        }
    }

    matches.into_iter()
}

