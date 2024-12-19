use std::collections::{BTreeSet, VecDeque};

use petgraph::{
    graph::{EdgeIndex, Graph},
    EdgeType,
};

pub fn is_subset_connected<N, E, Ix>(g: &Graph<N, E, Ix>, s: &BTreeSet<EdgeIndex>) -> bool
where
    Ix: EdgeType,
{
    let mut visited = BTreeSet::new();
    let mut queue = BTreeSet::from([*s.iter().next().unwrap()]);
    while let Some(e) = queue.pop_first() {
        visited.insert(e);
        let (src, dst) = g.edge_endpoints(e).unwrap();
        let nl = g.neighbors(src).filter_map(|n| {
            g.find_edge(src, n)
                .filter(|f| *f != e && s.contains(&f) && !visited.contains(f))
        });

        let nr = g.neighbors(dst).filter_map(|n| {
            g.find_edge(dst, n)
                .filter(|f| *f != e && s.contains(&f) && !visited.contains(f))
        });
        queue.extend(nl);
        queue.extend(nr);
    }

    return visited.len() == s.len();
}

pub fn edge_induced_subgraph<N, E, Ix>(
    mut g: Graph<N, E, Ix>,
    s: &BTreeSet<EdgeIndex>,
) -> Graph<N, E, Ix>
where
    Ix: EdgeType,
{
    g.retain_edges(|_, e| s.contains(&e));
    g.retain_nodes(|f, n| f.neighbors(n).count() != 0);
    g
}
