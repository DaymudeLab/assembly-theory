use nauty_Traces_sys::{empty_graph, ADDONEARC, SETWORDSNEEDED};
use petgraph::{visit::EdgeRef, EdgeType, Graph};
use std::{collections::BTreeMap, ffi::c_int, hash::Hash};

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct DenseGraph<N> {
    pub g: Vec<u64>,
    pub n: usize,
    pub e: usize,
    pub m: usize,
    pub nodes: Nodes<N>,
}
impl<N> DenseGraph<N> {
    pub fn from_petgraph<E, Ty>(graph: &Graph<N, E, Ty>) -> Self
    where
        Ty: EdgeType,
        N: Ord + Eq + Clone,
    {
        let n = graph.node_count();
        let e = graph.edge_count();
        let m = SETWORDSNEEDED(n);
        let nodes = Nodes::new(graph);
        let mut g = empty_graph(m, n);
        for edge in graph.edge_references() {
            let src = edge.source().index();
            let dst = edge.target().index();
            ADDONEARC(&mut g, src, dst, m);
            if !graph.is_directed() {
                ADDONEARC(&mut g, dst, src, m);
            }
        }
        Self { g, n, e, m, nodes }
    }

    pub fn orbits(&self) -> &[i32] {
        &self.nodes.orbits
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct Nodes<N> {
    pub lab: Vec<c_int>,
    pub ptn: Vec<c_int>,
    pub orbits: Vec<c_int>,
    pub weights: Vec<N>,
}

impl<N> Nodes<N> {
    pub fn new<E, Ty>(graph: &Graph<N, E, Ty>) -> Self
    where
        Ty: EdgeType,
        N: Ord + Eq + Clone,
    {
        let mut buckets = BTreeMap::new();
        for (ix, weight) in graph.node_weights().enumerate() {
            let bucket = buckets.entry(weight.clone()).or_insert(vec![]);
            bucket.push(ix as i32);
        }

        let mut lab = vec![];
        let mut ptn = vec![];
        let mut weights = vec![];
        for (weight, bucket) in buckets {
            ptn.extend(vec![1; bucket.len() - 1]);
            ptn.push(0);
            lab.extend(bucket);
            weights.push(weight)
        }
        Self {
            lab,
            ptn,
            orbits: vec![0; graph.node_count()],
            weights,
        }
    }
}
