use std::collections::BTreeSet;

use petgraph::{
    graph::{Edge, EdgeIndex, Graph, NodeIndex},
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

pub fn connected_components_under<N, E, Ty>(
    g: &Graph<N, E, Ty>,
    s: &BTreeSet<NodeIndex>,
) -> impl Iterator<Item = BTreeSet<NodeIndex>>
where
    Ty: EdgeType,
{
    let mut remainder = s.clone();
    let mut components = Vec::new();
    while !remainder.is_empty() {
        let mut visited = BTreeSet::new();
        let mut queue = BTreeSet::from([*remainder.iter().next().unwrap()]);
        while let Some(v) = queue.pop_first() {
            visited.insert(v);
            let neighbors = g
                .neighbors(v)
                .filter(|n| !visited.contains(n) && s.contains(n));
            queue.extend(neighbors)
        }
        remainder = remainder.difference(&visited).cloned().collect();
        components.push(visited);
    }
    components.into_iter()
}

pub fn connected_components_under_edges<N, E, Ty>(
    g: &Graph<N, E, Ty>,
    s: &BTreeSet<EdgeIndex>,
) -> impl Iterator<Item = BTreeSet<EdgeIndex>>
where
    Ty: EdgeType,
{
    let mut remainder = s.clone();
    let mut components = Vec::new();
    while !remainder.is_empty() {
        let mut visited = BTreeSet::new();
        let mut queue = BTreeSet::from([*remainder.iter().next().unwrap()]);
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
        remainder = remainder.difference(&visited).cloned().collect();
        components.push(visited);
    }
    components.into_iter()
}

pub fn edge_seperator<N, E, Ty>(
    g: &Graph<N, E, Ty>,
    s: &BTreeSet<NodeIndex>,
) -> (BTreeSet<EdgeIndex>, BTreeSet<EdgeIndex>)
where
    Ty: EdgeType,
{
    let left = g.edge_indices().filter(|e| {
        let (src, dst) = g.edge_endpoints(*e).unwrap();
        s.contains(&src) && s.contains(&dst)
    });

    let right = g.edge_indices().filter(|e| {
        let (src, dst) = g.edge_endpoints(*e).unwrap();
        !s.contains(&src) && !s.contains(&dst)
    });
    (BTreeSet::from_iter(left), BTreeSet::from_iter(right))
}

pub fn edges_contained_within<'a, N, E, Ty>(
    g: &'a Graph<N, E, Ty>,
    s: &'a BTreeSet<NodeIndex>,
) -> impl Iterator<Item = EdgeIndex> + 'a
where
    Ty: EdgeType,
{
    g.edge_indices().filter(|e| {
        let (src, dst) = g.edge_endpoints(*e).unwrap();
        s.contains(&src) && s.contains(&dst)
    })
}

pub fn edges_incident_to<'a, N, E, Ty>(
    g: &'a Graph<N, E, Ty>,
    s: &'a BTreeSet<NodeIndex>,
) -> impl Iterator<Item = EdgeIndex> + 'a
where
    Ty: EdgeType,
{
    g.edge_indices().filter(|e| {
        let (src, dst) = g.edge_endpoints(*e).unwrap();
        s.contains(&src) || s.contains(&dst)
    })
}
