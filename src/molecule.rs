use std::collections::{HashMap, HashSet};

use petgraph::{
    graph::{Graph, NodeIndex}, Undirected
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

        for ix in self.graph.node_indices() {
            if !u_set.contains(&ix) {
                io_map.insert(ix, ix);
            }
        }

        for ix in other.graph.node_indices() {
            if !v_set.contains(&ix) {
                io_map.insert(ix, ix);
            }
        }

        println!("{:?}", io_map);

        for ix in self.graph.node_indices() {
            let w = self.graph.node_weight(ix);
        }

        None
    }

    // Iterator over every joinable set of n vertices in self
    // pub fn n_join_points(&self, n: usize) -> impl Iterator<Item = impl Iterator<Item = Index>> {}

    pub fn from_graph(g: Graph<Atom, Bond, Undirected, Index>) -> Self {
        Self { graph: g }
    }
}
