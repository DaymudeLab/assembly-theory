use std::collections::{BTreeSet, HashMap, HashSet};

use petgraph::{
    algo::{is_isomorphic, is_isomorphic_subgraph, subgraph_isomorphisms_iter},
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::utils::{edge_induced_subgraph, is_subset_connected};

type Index = u32;
type MGraph = Graph<Atom, Bond, Undirected, Index>;

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
    graph: MGraph,
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
    /// Join self with other on `on`
    pub fn join(
        &self,
        other: &Molecule,
        on: impl IntoIterator<Item = (NodeIndex, NodeIndex)>,
    ) -> Option<Molecule> {
        let mut output_graph = self.clone();

        let mut v_set = HashSet::<NodeIndex>::new();
        let mut io_map = HashMap::<NodeIndex, NodeIndex>::new();

        for (u, v) in on.into_iter() {
            v_set.insert(v);
            io_map.insert(v, u);
        }

        for ix in other.graph.node_indices() {
            if !v_set.contains(&ix) {
                let w = other.graph.node_weight(ix).unwrap().clone();
                let out = output_graph.graph.add_node(w);
                io_map.insert(ix, out);
            }
        }

        for ix in other.graph.edge_indices() {
            let (u, v) = other.graph.edge_endpoints(ix).unwrap();
            let um = io_map.get(&u).unwrap();
            let vm = io_map.get(&v).unwrap();
            let w = other.graph.edge_weight(ix).unwrap().clone();

            output_graph.graph.add_edge(*um, *vm, w);
        }

        Some(output_graph)
    }

    pub fn is_isomorphic_to(&self, other: &Molecule) -> bool {
        is_isomorphic(&self.graph, &other.graph)
    }

    pub fn is_subgraph_of(&self, other: &Molecule) -> bool {
        is_isomorphic_subgraph(&self.graph, &other.graph)
    }

    pub fn enumerate_subgraphs(&self) -> impl Iterator<Item = BTreeSet<NodeIndex>> {
        let mut solutions = HashSet::new();
        let remainder = BTreeSet::from_iter(self.graph.node_indices());
        self.generate_connected_subgraphs(
            remainder,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        solutions.into_iter()
    }

    // From
    // https://stackoverflow.com/a/15722579
    // https://stackoverflow.com/a/15658245
    fn generate_connected_subgraphs(
        &self,
        mut remainder: BTreeSet<NodeIndex>,
        mut subset: BTreeSet<NodeIndex>,
        mut neighbors: BTreeSet<NodeIndex>,
        solutions: &mut HashSet<BTreeSet<NodeIndex>>,
    ) {
        let candidates = if subset.is_empty() {
            remainder.clone()
        } else {
            remainder.intersection(&neighbors).cloned().collect()
        };

        if candidates.is_empty() {
            solutions.insert(subset);
        } else {
            let v = candidates.first().unwrap();
            remainder.remove(v);

            self.generate_connected_subgraphs(
                remainder.clone(),
                subset.clone(),
                neighbors.clone(),
                solutions,
            );

            subset.insert(*v);
            neighbors.extend(self.graph.neighbors(*v));
            self.generate_connected_subgraphs(remainder, subset, neighbors, solutions);
        }
    }

    pub fn partitions(&self) -> Option<impl Iterator<Item = (Molecule, Molecule)> + '_> {
        let mut solutions = HashSet::new();
        let remaining_edges = self.graph.edge_indices().collect();
        self.backtrack(
            remaining_edges,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        Some(solutions.into_iter().map(|(left, right)| {
            (
                Molecule {
                    graph: edge_induced_subgraph(self.graph.clone(), &left),
                },
                Molecule {
                    graph: edge_induced_subgraph(self.graph.clone(), &right),
                },
            )
        }))
    }

    fn backtrack(
        &self,
        mut remaining_edges: Vec<EdgeIndex>,
        left: BTreeSet<EdgeIndex>,
        right: BTreeSet<EdgeIndex>,
        solutions: &mut HashSet<(BTreeSet<EdgeIndex>, BTreeSet<EdgeIndex>)>,
    ) {
        if remaining_edges.is_empty() {
            if self.is_valid_partition(&left, &right) {
                solutions.insert((left, right));
            }
            return;
        }

        let suffix = remaining_edges.pop().unwrap();
        let mut lc = left.clone();
        lc.insert(suffix);

        let mut rc = right.clone();
        rc.insert(suffix);

        self.backtrack(remaining_edges.clone(), lc, right, solutions);
        self.backtrack(remaining_edges, left, rc, solutions);
    }

    fn is_valid_partition(&self, left: &BTreeSet<EdgeIndex>, right: &BTreeSet<EdgeIndex>) -> bool {
        !left.is_empty()
            && !right.is_empty()
            && is_subset_connected(&self.graph, left)
            && is_subset_connected(&self.graph, right)
    }

    pub fn from_graph(g: MGraph) -> Self {
        Self { graph: g }
    }

    pub fn single_bond() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        g.add_edge(u, v, Bond::Single);
        Self { graph: g }
    }

    pub fn double_bond() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        g.add_edge(u, v, Bond::Double);
        Self { graph: g }
    }

    pub fn carbonyl() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Double);
        Self { graph: g }
    }

    pub fn hydroxyl() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Single);
        Self { graph: g }
    }

    pub fn is_basic_unit(&self) -> bool {
        self.is_isomorphic_to(&Molecule::hydroxyl())
            || self.is_isomorphic_to(&Molecule::carbonyl())
            || self.is_isomorphic_to(&Molecule::single_bond())
            || self.is_isomorphic_to(&Molecule::double_bond())
    }
}

mod tests {
    use super::*;

    #[test]
    fn join_carbonyl_hydroxyl() {
        let c = Molecule::carbonyl();
        let h = Molecule::hydroxyl();
        let j = c
            .join(&h, [(NodeIndex::new(0), NodeIndex::new(0))])
            .unwrap();

        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        let w = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Single);
        g.add_edge(u, w, Bond::Double);

        assert!(j.is_isomorphic_to(&Molecule { graph: g }));
    }
}
