use petgraph::graph::{Graph,Undirected};
use std::fmt::Debug;
use std::hash::Hash;

pub trait AObject {
    type NodeLabel: Debug + Clone + Eq;
    type EdgeLabel: Debug + Clone + Eq;

}