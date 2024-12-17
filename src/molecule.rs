use std::collections::{HashMap, HashSet};

use petgraph::{
    data::{Build, DataMap},
    graph::{Graph, NodeIndex},
    Undirected,
};

type Index = u32;

#[derive(Debug, Clone)]
pub enum Element {
    Hydrogen,
    Carbon,
    Nitrogen,
    Oxygen,
}

#[derive(Debug, Clone)]
pub struct Atom {
    element: Element,
    capacity: u32,
}

#[derive(Debug, Clone)]
pub enum Bond {
    Single,
    Double,
}

#[derive(Debug, Clone)]
pub struct Molecule {
    graph: Graph<Atom, Bond, Undirected, Index>,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            capacity: 0,
        }
    }
}

impl Molecule {
    pub fn join(
        &self,
        other: &Molecule,
        on: impl IntoIterator<Item = (NodeIndex, NodeIndex)>,
    ) -> Option<Molecule> {
        let mut output_graph = self.clone();

        let mut u_set = HashSet::<NodeIndex>::new();
        let mut v_set = HashSet::<NodeIndex>::new();
        let mut io_map = HashMap::<NodeIndex, NodeIndex>::new();

        for (u, v) in on.into_iter() {
            u_set.insert(u);
            v_set.insert(v);
            io_map.insert(v, u);
        }

        for ix in other.graph.node_indices() {
            if !v_set.contains(&ix) {
                let out = output_graph
                    .graph
                    .add_node(other.graph.node_weight(ix).unwrap().clone());
                io_map.insert(ix, out);
            }
        }

        for ix in other.graph.edge_indices() {
            let (u, v) = other.graph.edge_endpoints(ix).unwrap();
            let w = other.graph.edge_weight(ix).unwrap().clone();
            let um = io_map.get(&u).unwrap();
            let vm = io_map.get(&v).unwrap();

            output_graph.graph.add_edge(*um, *vm, w);
        }

        Some(output_graph)
    }

    // Iterator over every joinable set of n vertices in self
    // pub fn n_join_points(&self, n: usize) -> impl Iterator<Item = impl Iterator<Item = Index>> {}

    pub fn from_graph(g: Graph<Atom, Bond, Undirected, Index>) -> Self {
        Self { graph: g }
    }
}
