//! Define the abstract object trait and associated types.
//! 
use petgraph::{
     graph::{Graph},
    Undirected,
    dot::Dot
};
use std::collections::HashSet;
use std::fmt::Debug;
use std::marker::{Copy, Sync};
use std::hash::Hash;

pub(crate) type Index = u32;

/// Methods required of nodes
pub trait NodeBehavior {
    fn kind(&self) -> u8;
}
/// Methods required of edges
pub trait EdgeBehavior {
    fn kind(&self) -> u8;
}

/// Abstract representation of an object.
pub trait AObject : Sync + Send {
    type NodeLabel: NodeBehavior + Copy + Sync + Debug + Eq + Hash;
    type EdgeLabel: EdgeBehavior + Copy + Sync + Debug + Eq + Hash;

    fn from_graph(g: Graph<Self::NodeLabel, Self::EdgeLabel, Undirected, Index>) -> Self;

    /// Return a representation of this molecule as an undirected Graph 
    /// with node and edge labels.
    fn graph(&self) -> &Graph<Self::NodeLabel, Self::EdgeLabel, Undirected, Index>;

    /// Return a pretty-printable representation of this object.
    fn info(&self) -> String {
        let dot = Dot::new(self.graph());
        format!("{dot:?}")
    }

    /// Return `true` iff this object contains self-loops or multiple edges
    /// between any pair of nodes. TODO: Check the logic here, need better testing
    fn is_malformed(&self) -> bool {
        let mut uniq = HashSet::new();
        !self.graph().edge_indices().all(|ix| {
            uniq.insert(ix)
                && self
                    .graph()
                    .edge_endpoints(ix)
                    .is_some_and(|(src, dst)| src != dst)
        })
    }

    /// Return `true` iff this object comprises only one edge (of any type).
    fn is_basic_unit(&self) -> bool {
        self.graph().edge_count() == 1 && self.graph().node_count() == 2
    }

}