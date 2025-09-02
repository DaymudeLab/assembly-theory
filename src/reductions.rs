//! Reduces the optimal molecular assembly index problem into maximum weighted clique.
//! 
//! Removing a duplicatable subgraph pair (h1, h2) gives a savings of |h1| - 1.
//! The problem of finding the optimal molecular assembly index can thus be thought of as finding
//! a set of duplicatable subgraph removals that are all compatible with each other (i.e. can be removed
//! at the same time) with the largest savings possible. Thus, this problem is close to that of finding
//! a maximum clique in a graph with duplicatable subgraph pairs for vertices, vertex weights |h1| - 1, and
//! edges representing compatibility. However, the order of duplicatable subgraph removals can affect
//! compatibility, and so edges must be directed. By ordering the vertices in a nice way, we can ignore
//! edge directionality.

use bit_set::BitSet;
use std::time::Instant;
use std::collections::HashMap;

/// Graph representation of the compatibility of duplicatable subgraph pairs.
pub struct CompatGraph {
    /// The graph is implemented as an adjacency matrix. The ith element of graph
    /// has a 1 at position j if {i, j} is an edge.
    graph: Vec<BitSet>,
}

impl CompatGraph {
    /// Constructs a compatibility graph given a set of matches.
    pub fn new(init_matches: &Vec<(usize, usize)>, dag: &HashMap<usize, (BitSet, Vec<usize>)>) -> Self {
        let start = Instant::now();
        let size = init_matches.len();

        // Initialize weights and empty graph
        let mut init_graph: Vec<BitSet> = Vec::with_capacity(size);
        for _ in 0..init_matches.len() {
            init_graph.push(BitSet::with_capacity(size));
        }

        // Populate graph
        for (idx1, (h1, h2)) in init_matches.iter().enumerate() {
            for (idx2, (h1p, h2p)) in init_matches[idx1 + 1..].iter().enumerate() {
                let idx2 = idx2 + idx1 + 1;
                let h1 = &dag.get(h1).unwrap().0;
                let h1p = &dag.get(h1p).unwrap().0;
                let h2 = &dag.get(h2).unwrap().0;
                let h2p = &dag.get(h2p).unwrap().0;

                let forward_compatible = {
                    h2.is_disjoint(h1p) && 
                    h2.is_disjoint(h2p) &&
                    (h1.is_disjoint(h1p) || h1.is_superset(h1p)) &&
                    (h1.is_disjoint(h2p) || h1.is_superset(h2p))
                };

                if forward_compatible {
                    init_graph[idx1].insert(idx2);
                    init_graph[idx2].insert(idx1);
                }
            }
        }

        println!("Graph Time: {:?}", start.elapsed());

        Self {
            graph: init_graph,
        }
    }

    /// Returns the number of vertices in the graph.
    pub fn len(&self) -> usize {
        self.graph.len()
    }

    /// Returns the degree vertex v in the subgraph
    pub fn degree(&self, v: usize, subgraph: &BitSet) -> usize {
        self.graph[v].intersection(subgraph).count()
    }

    /// Returns the set of matches compatible with v. I.e. the set neighbors
    /// of v.
    pub fn compatible_with(&self, v: usize) -> &BitSet {
        &self.graph[v]
    }

    /// Returns the neighbors of v in the given subgraph.
    /// The returned BitSet is a clone, and thus can be safely modified by the calling function
    /// without changing the graph.
    pub fn neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        
        neighbors
    }

    /// Returns the neighbors of v that occur after v in the matches list.
    /// The returned BitSet is a clone, and thus can be safely modified by the calling function
    /// without changing the graph.
    pub fn forward_neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        let mut to_remove = vec![];
        for u in neighbors.iter() {
            if u <= v {
                to_remove.push(u);
            }
            if u > v {
                break;
            }
        }
        for u in to_remove {
            neighbors.remove(u);
        }

        neighbors
    }

    /// Returns true if vertices v and u are adjacent.
    pub fn are_adjacent(&self, v: usize, u: usize) -> bool {
        self.graph[v].contains(u)
    }
}