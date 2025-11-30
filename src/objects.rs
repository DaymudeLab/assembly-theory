//! Trait definition for AssemblyObject.

use petgraph::{
    graph::{EdgeIndex, Graph, NodeIndex},
    dot::Dot,
    Undirected,
};
use std::collections::{BTreeSet, HashSet};
use crate::utils::{edge_induced_subgraph, is_subset_connected};


pub(crate) type Index = u32;
pub(crate) type MGraph = Graph<u32, u32, Undirected, Index>;
type EdgeSet = BTreeSet<EdgeIndex<Index>>;

type NodeSet = BTreeSet<NodeIndex<Index>>;

/// Trait for objects that operate on assembly graphs.
pub trait AssemblyObject {
    /// Construct an object from an existing graph.
    fn from_graph(graph: MGraph) -> Self;

    /// Return a reference to the underlying graph.
    fn graph(&self) -> &MGraph;

    /// Return a pretty-printable representation of this object.
    fn info(&self) -> String {
        let dot = Dot::new(&self.graph);
        format!("{dot:?}")
    }
    /// Check if the graph contains self-loops or multiple edges.
    fn is_malformed(&self) -> bool {
        let mut uniq = HashSet::new();
        self.graph().edge_indices().all(|ix| {
            uniq.insert(ix)
                && self
                    .graph()
                    .edge_endpoints(ix)
                    .is_some_and(|(src, dst)| src != dst)
        })
    }

    /// Check if the graph is a basic unit (one edge, two nodes).
    fn is_basic_unit(&self) -> bool {
        self.graph().edge_count() == 1 && self.graph().node_count() == 2
    }

    /// Return all valid partitions of the graph.
    fn partitions(&self) -> Option<Vec<(Self, Self)>>
    where
        Self: Sized,
    {
        let mut solutions = HashSet::new();
        let remaining_edges = self.graph().edge_indices().collect();
        self.backtrack(remaining_edges, BTreeSet::new(), BTreeSet::new(), &mut solutions);
        Some(
            solutions
                .into_iter()
                .map(|(left, right)| {
                    (
                        Self::from_graph(edge_induced_subgraph(self.graph().clone(), &left)),
                        Self::from_graph(edge_induced_subgraph(self.graph().clone(), &right)),
                    )
                })
                .collect(),
        )
    }

    /// Helper method for recursive partitioning.
    fn backtrack(
        &self,
        mut remaining_edges: Vec<EdgeIndex<Index>>,
        left: EdgeSet,
        right: EdgeSet,
        solutions: &mut HashSet<(EdgeSet, EdgeSet)>,
    ) {
        if let Some(suffix) = remaining_edges.pop() {
            let mut lc = left.clone();
            lc.insert(suffix);

            let mut rc = right.clone();
            rc.insert(suffix);

            self.backtrack(remaining_edges.clone(), lc, right, solutions);
            self.backtrack(remaining_edges, left, rc, solutions);
        } else if self.is_valid_partition(&left, &right) {
            solutions.insert((left, right));
        }
    }

    /// Check if a partition is valid.
    fn is_valid_partition(&self, left: &EdgeSet, right: &EdgeSet) -> bool {
        !left.is_empty()
            && !right.is_empty()
            && is_subset_connected(self.graph(), left)
            && is_subset_connected(self.graph(), right)
    }
}
