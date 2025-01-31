use std::collections::{BTreeSet, HashMap, HashSet};

use bit_set::BitSet;
use petgraph::{
    algo::{is_isomorphic, is_isomorphic_subgraph, subgraph_isomorphisms_iter},
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::utils::{edge_induced_subgraph, edges_contained_within, is_subset_connected};

pub type Index = u32;
pub type MGraph = Graph<Atom, Bond, Undirected, Index>;
type MSubgraph = Graph<Atom, Option<Bond>, Undirected, Index>;
type EdgeSet = BTreeSet<EdgeIndex<Index>>;
type NodeSet = BTreeSet<NodeIndex<Index>>;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Element {
    Hydrogen,
    Carbon,
    Nitrogen,
    Oxygen,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Atom {
    element: Element,
    capacity: u32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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
        on: impl IntoIterator<Item = (NodeIndex<Index>, NodeIndex<Index>)>,
    ) -> Option<Molecule> {
        let mut output_graph = self.clone();

        let mut v_set = NodeSet::new();
        let mut io_map = HashMap::<NodeIndex<Index>, NodeIndex<Index>>::new();

        for (u, v) in on.into_iter() {
            v_set.insert(v);
            io_map.insert(v, u);
        }

        for ix in other.graph.node_indices() {
            if !v_set.contains(&ix) {
                let w = *other.graph.node_weight(ix)?;
                let out = output_graph.graph.add_node(w);
                io_map.insert(ix, out);
            }
        }

        for ix in other.graph.edge_indices() {
            let (u, v) = other.graph.edge_endpoints(ix)?;
            let um = io_map.get(&u)?;
            let vm = io_map.get(&v)?;
            let w = *other.graph.edge_weight(ix)?;

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

    pub fn enumerate_subgraphs(&self) -> impl Iterator<Item = NodeSet> {
        let mut solutions = HashSet::new();
        let remainder = BTreeSet::from_iter(self.graph.node_indices());
        self.generate_connected_subgraphs(
            remainder,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        solutions.into_iter().filter(|s| !s.is_empty())
    }

    // From
    // https://stackoverflow.com/a/15722579
    // https://stackoverflow.com/a/15658245
    fn generate_connected_subgraphs(
        &self,
        mut remainder: NodeSet,
        mut subset: NodeSet,
        mut neighbors: NodeSet,
        solutions: &mut HashSet<NodeSet>,
    ) {
        let candidates = if subset.is_empty() {
            remainder.clone()
        } else {
            remainder.intersection(&neighbors).cloned().collect()
        };

        if let Some(v) = candidates.first() {
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
        } else if subset.len() > 2 {
            solutions.insert(subset);
        }
    }

    pub fn matches(&self) -> impl Iterator<Item = (BitSet, BitSet)> {
        let mut matches = BTreeSet::new();
        for subgraph in self.enumerate_subgraphs() {
            let mut h = self.graph().clone();
            h.retain_nodes(|_, n| subgraph.contains(&n));

            let h_prime = self.graph().map(
                |_, n| *n,
                |i, e| {
                    let (src, dst) = self.graph.edge_endpoints(i).unwrap();
                    (!subgraph.contains(&src) || !subgraph.contains(&dst)).then_some(*e)
                },
            );

            for cert in isomorphic_subgraphs_of(&h, &h_prime) {
                let cert = BTreeSet::from_iter(cert);

                let mut c = BitSet::new();
                c.extend(edges_contained_within(&self.graph, &cert).map(|e| e.index()));

                let mut h = BitSet::new();
                h.extend(edges_contained_within(&self.graph, &subgraph).map(|e| e.index()));

                matches.insert(if c < h { (c, h) } else { (h, c) });
            }
        }

        matches.into_iter()
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

    fn is_valid_partition(
        &self,
        left: &EdgeSet,
        right: &EdgeSet,
    ) -> bool {
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

    pub fn graph(&self) -> &MGraph {
        &self.graph
    }
}

// This is bad, BIG TODO: replace with handwritten vf3
pub fn isomorphic_subgraphs_of(pattern: &MGraph, target: &MSubgraph) -> Vec<Vec<NodeIndex<Index>>> {
    if let Some(iter) =
        subgraph_isomorphisms_iter(&pattern, &target, &mut |n0, n1| n1 == n0, &mut |e0, e1| {
            e1.is_some_and(|e| *e0 == e)
        })
    {
        iter.map(|v| v.into_iter().map(NodeIndex::new).collect())
            .collect()
    } else {
        Vec::new()
    }
}

mod tests {
    #![allow(unused_imports)]
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
