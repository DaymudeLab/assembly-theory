use crate::nauty::dense::DenseGraph;
use bitvec::prelude::*;
use nauty_Traces_sys::{densenauty, empty_graph, optionblk, statsblk};
use petgraph::{EdgeType, Graph};
use std::{
    hash::{Hash, Hasher},
    os::raw::c_int,
};

#[derive(Eq, Debug, Clone)]
pub struct CanonLabeling<N> {
    pub g: Vec<u64>,
    pub e: usize,
    pub n: usize,
    dense: DenseGraph<N>,
}
impl<N> Hash for CanonLabeling<N>
where
    N: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.g.hash(state);
        self.e.hash(state);
        self.n.hash(state);
        self.dense.nodes.ptn.hash(state);
        self.dense.nodes.weights.hash(state);
    }
}

impl<N> PartialEq for CanonLabeling<N>
where
    N: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.g == other.g
            && self.e == other.e
            && self.n == other.n
            && self.dense.nodes.ptn == other.dense.nodes.ptn
            && self.dense.nodes.weights == other.dense.nodes.weights
    }
}

impl<N> PartialOrd for CanonLabeling<N>
where
    N: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(
            self.n
                .cmp(&other.n)
                .then(self.e.cmp(&other.e))
                .then(self.g.cmp(&other.g))
                .then(self.dense.nodes.ptn.partial_cmp(&other.dense.nodes.ptn)?)
                .then(
                    self.dense
                        .nodes
                        .weights
                        .partial_cmp(&other.dense.nodes.weights)?,
                ),
        )
    }
}

impl<N> Ord for CanonLabeling<N>
where
    N: Ord,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.n
            .cmp(&other.n)
            .then(self.e.cmp(&other.e))
            .then(self.g.cmp(&other.g))
            .then(self.dense.nodes.ptn.cmp(&other.dense.nodes.ptn))
            .then(self.dense.nodes.weights.cmp(&other.dense.nodes.weights))
    }
}

impl<N> CanonLabeling<N> {
    pub fn new<E, Ty>(graph: &Graph<N, E, Ty>) -> Self
    where
        Ty: EdgeType,
        N: Ord + Eq + Clone,
    {
        let mut dg = DenseGraph::<N>::from_petgraph(graph);
        let mut opt = canon_opts(graph.is_directed());
        let mut stat = statsblk::default();
        let mut cg = empty_graph(dg.m, dg.n);
        unsafe {
            densenauty(
                dg.g.as_mut_ptr(),
                dg.nodes.lab.as_mut_ptr(),
                dg.nodes.ptn.as_mut_ptr(),
                dg.nodes.orbits.as_mut_ptr(),
                &mut opt,
                &mut stat,
                dg.m as c_int,
                dg.n as c_int,
                cg.as_mut_ptr(),
            )
        }
        Self {
            g: cg,
            e: dg.e,
            n: dg.n,
            dense: dg,
        }
    }

    /// Returns the adjacency matrix as a flat vector
    pub fn flat_adjacency(&self) -> Vec<usize> {
        let mut bit_vector = Vec::with_capacity(self.n * self.n);
        for num in self.g.iter() {
            let bv = num.view_bits::<Msb0>();
            for bit in bv.iter().take(self.n) {
                if *bit {
                    bit_vector.push(1);
                } else {
                    bit_vector.push(0);
                }
            }
        }
        bit_vector
    }

    pub fn orbits(&self) -> &[i32] {
        self.dense.orbits()
    }
}

/// Creates a `optionblk` struct for canonization
///
/// # Arguments
/// * `is_directed` - Whether the graph is directed
pub fn canon_opts(is_directed: bool) -> optionblk {
    optionblk {
        getcanon: 1,
        digraph: is_directed.into(),
        defaultptn: 0,
        ..Default::default()
    }
}

/// Creates an edge list from a bit adjacency matrix
///
/// # Arguments
/// * `adj` - A bit adjacency matrix
/// * `e` - The number of edges
/// * `n` - The number of nodes
pub fn bit_adj_to_edgelist(adj: &[u64], e: usize, n: usize) -> Vec<(u32, u32)> {
    let mut edges = Vec::with_capacity(e);
    for (idx, num) in adj.iter().enumerate() {
        let bv = num.view_bits::<Msb0>();
        for (jdx, bit) in bv.iter().enumerate().take(n) {
            if *bit {
                edges.push((idx as u32, jdx as u32));
            }
        }
    }
    edges
}

/// Creates a graph from a bit adjacency matrix
///
/// # Arguments
/// * `adj` - A bit adjacency matrix
/// * `e` - The number of edges
/// * `n` - The number of nodes
pub fn bit_adj_to_graph<Ty>(adj: &[u64], e: usize, n: usize) -> Graph<(), (), Ty>
where
    Ty: EdgeType,
{
    let edges = bit_adj_to_edgelist(adj, e, n);
    Graph::from_edges(&edges)
}

#[cfg(test)]
mod testing {
    use petgraph::{Directed, Graph, Undirected};

    #[test]
    fn test_equivalent_digraph() {
        let e1 = vec![(0, 1), (0, 2), (1, 2)];
        let e2 = vec![(1, 0), (1, 2), (0, 2)];

        let g1: Graph<(), (), Directed> = Graph::from_edges(&e1);
        let g2: Graph<(), (), Directed> = Graph::from_edges(&e2);

        let l1 = super::CanonLabeling::new(&g1);
        let l2 = super::CanonLabeling::new(&g2);

        assert_eq!(l1, l2);
    }

    #[test]
    fn test_unequal_digraph() {
        let e1 = vec![(0, 1), (0, 2), (1, 2)];
        let e2 = vec![(1, 0), (1, 2), (2, 1)];

        let g1: Graph<(), (), Directed> = Graph::from_edges(&e1);
        let g2: Graph<(), (), Directed> = Graph::from_edges(&e2);

        let l1 = super::CanonLabeling::new(&g1);
        let l2 = super::CanonLabeling::new(&g2);

        assert_ne!(l1, l2);
    }

    #[test]
    fn test_equal_ungraph() {
        let e1 = vec![(0, 1), (0, 2), (1, 2)];
        let e2 = vec![(1, 0), (1, 2), (0, 2)];

        let g1: Graph<(), (), Undirected> = Graph::from_edges(&e1);
        let g2: Graph<(), (), Undirected> = Graph::from_edges(&e2);

        let l1 = super::CanonLabeling::new(&g1);
        let l2 = super::CanonLabeling::new(&g2);

        assert_eq!(l1, l2);
    }

    #[test]
    fn test_unequal_ungraph() {
        let e1 = vec![(0, 1), (0, 2), (1, 2)];
        let e2 = vec![(1, 0), (1, 2)];

        let g1: Graph<(), (), Undirected> = Graph::from_edges(&e1);
        let g2: Graph<(), (), Undirected> = Graph::from_edges(&e2);

        let l1 = super::CanonLabeling::new(&g1);
        let l2 = super::CanonLabeling::new(&g2);

        assert_ne!(l1, l2);
    }

    #[test]
    fn test_unequal_labeled_graph() {
        let mut g1 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g1.add_node(0);
        let n1 = g1.add_node(1);
        let n2 = g1.add_node(0);
        g1.add_edge(n0, n1, ());
        g1.add_edge(n1, n2, ());

        let mut g2 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g2.add_node(0);
        let n1 = g2.add_node(0);
        let n2 = g2.add_node(1);
        g2.add_edge(n0, n1, ());
        g2.add_edge(n1, n2, ());

        let repr_a = super::CanonLabeling::new(&g1);
        let repr_b = super::CanonLabeling::new(&g2);

        assert_ne!(repr_a, repr_b);
    }

    #[test]
    fn test_unequal_labeled_graph_again() {
        let mut g1 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g1.add_node(0);
        let n1 = g1.add_node(1);
        let n2 = g1.add_node(0);
        g1.add_edge(n0, n1, ());
        g1.add_edge(n1, n2, ());

        let mut g2 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g2.add_node(0);
        let n1 = g2.add_node(0);
        let n2 = g2.add_node(0);
        g2.add_edge(n0, n1, ());
        g2.add_edge(n1, n2, ());

        let repr_a = super::CanonLabeling::new(&g1);
        let repr_b = super::CanonLabeling::new(&g2);

        assert_ne!(repr_a, repr_b);
    }

    #[test]
    fn test_equal_labeled_graph() {
        let mut g1 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g1.add_node(1);
        let n1 = g1.add_node(0);
        let n2 = g1.add_node(0);
        g1.add_edge(n0, n1, ());
        g1.add_edge(n1, n2, ());

        let mut g2 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = g2.add_node(0);
        let n1 = g2.add_node(0);
        let n2 = g2.add_node(1);
        g2.add_edge(n0, n1, ());
        g2.add_edge(n1, n2, ());

        let repr_a = super::CanonLabeling::new(&g1);
        let repr_b = super::CanonLabeling::new(&g2);

        assert_eq!(repr_a, repr_b);
    }

    #[test]
    fn test_label() {
        let edges = vec![(0, 1), (0, 2), (1, 2)];
        let graph: Graph<(), (), Directed> = Graph::from_edges(&edges);
        let canon = super::CanonLabeling::new(&graph);
        assert_eq!(canon.g, vec![0, 9223372036854775808, 13835058055282163712]);
        assert_eq!(canon.e, 3);
        assert_eq!(canon.n, 3);
    }

    #[test]
    fn test_flat_adj_directed() {
        let edges = vec![(0, 1), (0, 2), (1, 2)];
        let graph: Graph<(), (), Directed> = Graph::from_edges(&edges);
        let canon = super::CanonLabeling::<()>::new(&graph);
        let flat_adj = canon.flat_adjacency();
        assert_eq!(flat_adj, vec![0, 0, 0, 1, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_flat_adj_undirected() {
        let edges = vec![(0, 1), (0, 2), (1, 2)];
        let graph: Graph<(), (), Undirected> = Graph::from_edges(&edges);
        let canon = super::CanonLabeling::new(&graph);
        let flat_adj = canon.flat_adjacency();
        assert_eq!(flat_adj, vec![0, 1, 1, 1, 0, 1, 1, 1, 0]);
    }
}
