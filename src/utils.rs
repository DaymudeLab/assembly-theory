use std::collections::{BTreeSet, VecDeque};

use petgraph::{
    graph::{EdgeIndex, Graph},
    EdgeType,
};

pub fn is_subset_connected<N, E, Ix>(g: &Graph<N, E, Ix>, s: &BTreeSet<EdgeIndex>) -> bool
where
    Ix: EdgeType,
{
    let mut visited = 0;
    let mut queue = VecDeque::from([*s.iter().next().unwrap()]);
    while let Some(e) = queue.pop_front() {
        let (src, dst) = g.edge_endpoints(e).unwrap();
        queue.extend(
            g.neighbors(src)
                .filter_map(|n| g.find_edge(src, n).filter(|f| *f != e && s.contains(&e))),
        );
        queue.extend(
            g.neighbors(dst)
                .filter_map(|n| g.find_edge(src, n).filter(|f| *f != e && s.contains(&e))),
        );
        visited += 1;
    }
    return visited == g.node_count();
}

pub fn edge_induced_subgraph<N, E, Ix>(g: Graph<N, E, Ix>, s: &BTreeSet<EdgeIndex>) -> Graph<N, E, Ix> {
    g
}
