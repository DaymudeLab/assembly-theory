use std::collections::BTreeSet;

use petgraph::{
    graph::{EdgeIndex, Graph, NodeIndex},
    EdgeType,
};

pub fn is_subset_connected<N, E, Ty>(g: &Graph<N, E, Ty>, s: &BTreeSet<EdgeIndex>) -> bool
where
    Ty: EdgeType,
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

pub fn edge_induced_subgraph<N, E, Ty>(
    mut g: Graph<N, E, Ty>,
    s: &BTreeSet<EdgeIndex>,
) -> Graph<N, E, Ty>
where
    Ty: EdgeType,
{
    g.retain_edges(|_, e| s.contains(&e));
    g.retain_nodes(|f, n| f.neighbors(n).count() != 0);
    g
}

pub fn node_induced_subgraph<N, E, Ty>(
    mut g: Graph<N, E, Ty>,
    s: &BTreeSet<NodeIndex>,
) -> Graph<N, E, Ty>
where
    Ty: EdgeType,
{
    g.retain_nodes(|_, n| s.contains(&n));
    g
}

pub fn edge_induced_cosubgraph<N, E, Ty>(
    mut g: Graph<N, E, Ty>,
    s: &BTreeSet<EdgeIndex>,
) -> Graph<N, E, Ty>
where
    Ty: EdgeType,
{
    g.retain_edges(|_, e| !s.contains(&e));
    g.retain_nodes(|f, n| f.neighbors(n).count() != 0);
    g
}
