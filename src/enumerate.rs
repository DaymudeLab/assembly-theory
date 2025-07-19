//! Enumerate connected, non-induced subgraphs of a molecular graph.
//!
//! Specifically, for a molecule with |E| edges, enumerate all connected,
//! non-induced subgraphs with at most |E|/2 edges. Any larger subgraphs
//! cannot be "duplicatable" (i.e., in a pair of non-overlapping, isomorphic
//! subgraphs), so we don't need them.

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
    /// Like Extend, but only extends subgraphs that are isomoprhic to at least
    /// one other subgraph.
    ExtendIsomorphic,
    /// From a subgraph, choose an edge from its boundary and either grow it by
    /// adding this edge or erode its remainder by discarding the edge.
    GrowErode,
}

/// Return an interator over all connected, non-induced subgraphs of the
/// molecular graph `mol` using the algorithm specified by `mode`.
pub fn enumerate_subgraphs(mol: &Molecule, mode: EnumerateMode) -> impl Iterator<Item = BitSet> {
    match mode {
        EnumerateMode::Extend => extend(mol).into_iter(),
        EnumerateMode::GrowErode => grow_erode(mol).into_iter(),
        _ => {
            panic!("The chosen --enumerate mode is not implemented yet!");
        }
    }
}

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges using
/// a process of iterative extension starting from each individual edge.
fn extend(mol: &Molecule) -> HashSet<BitSet> {
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
            // Find all "frontier" edges incident to this subgraph (contains
            // both this subgraph's edges and its edge boundary).
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

/// Enumerate connected, non-induced subgraphs with at most |E|/2 edges which
/// are isomorphic to at least one other subgraph (i.e., the subgraphs that
/// have a chance of being in a non-overlapping, isomorphic subgraph pair later
/// on). Uses a similar iterative extension process to extend_iterative, but at
/// each level, removes any subgraphs in singleton isomorphism classes. 
fn extend_isomorphic(mol: &Molecule) -> impl Iterator<Item = (BitSet, BitSet)> {
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
                frontier.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|ix| ix.index()));
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
