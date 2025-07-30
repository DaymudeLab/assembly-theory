#![allow(dead_code)]

use std::collections::BTreeSet;

use bit_set::BitSet;
use petgraph::{
    graph::{EdgeIndex, Graph, IndexType, NodeIndex},
    EdgeType,
};

pub fn is_subset_connected<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    s: &BTreeSet<EdgeIndex<Ix>>,
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut visited = BTreeSet::new();
    let mut queue = BTreeSet::from([*s.iter().next().unwrap()]);
    while let Some(e) = queue.pop_first() {
        visited.insert(e);
        let (src, dst) = g.edge_endpoints(e).unwrap();
        let nl = g.neighbors(src).filter_map(|n| {
            g.find_edge(src, n)
                .filter(|f| *f != e && s.contains(f) && !visited.contains(f))
        });

        let nr = g.neighbors(dst).filter_map(|n| {
            g.find_edge(dst, n)
                .filter(|f| *f != e && s.contains(f) && !visited.contains(f))
        });
        queue.extend(nl);
        queue.extend(nr);
    }

    visited.len() == s.len()
}

pub fn edge_induced_subgraph<N, E, Ty, Ix>(
    mut g: Graph<N, E, Ty, Ix>,
    s: &BTreeSet<EdgeIndex<Ix>>,
) -> Graph<N, E, Ty, Ix>
where
    Ty: EdgeType,
    Ix: IndexType,
{
    g.retain_edges(|_, e| s.contains(&e));
    g.retain_nodes(|f, n| f.neighbors(n).count() != 0);
    g
}

pub fn subgraph_from_edge_mask<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    s: &BitSet,
) -> Graph<N, E, Ty, Ix>
where
    Ty: EdgeType,
    Ix: IndexType,
    N: Clone,
    E: Clone,
{
    let mut g = g.clone();
    g.retain_edges(|_, e| s.contains(e.index()));
    g.retain_nodes(|f, n| f.neighbors(n).count() != 0);
    g
}

pub fn node_count_under_edge_mask<N, E, Ty, Ix>(g: &Graph<N, E, Ty, Ix>, s: &BitSet) -> usize
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut node_set = BitSet::with_capacity(g.node_count());
    for ix in s.into_iter().map(|ix| EdgeIndex::new(ix)) {
        let (src, dst) = g.edge_endpoints(ix).expect("malformed bitset!");
        node_set.insert(src.index());
        node_set.insert(dst.index());
    }
    node_set.len()
}

pub fn node_induced_subgraph<N, E, Ty, Ix>(
    mut g: Graph<N, E, Ty, Ix>,
    s: &BTreeSet<NodeIndex<Ix>>,
) -> Graph<N, E, Ty, Ix>
where
    Ty: EdgeType,
    Ix: IndexType,
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

pub fn connected_components_under_edges<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    s: &BitSet,
) -> impl Iterator<Item = BitSet>
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut remainder = s.clone();
    let mut components = Vec::new();
    while let Some(c) = remainder.iter().next() {
        let mut queue = BitSet::new();
        queue.reserve_len(s.len());
        queue.insert(c);

        let mut visited = BitSet::new();
        visited.reserve_len(s.len());

        while let Some(e) = queue.iter().next() {
            queue.remove(e);
            visited.insert(e);
            let (src, dst) = g.edge_endpoints(EdgeIndex::new(e)).unwrap();

            let nl = g.neighbors(src).filter_map(|n| {
                g.find_edge(src, n).and_then(|f| {
                    let ix = f.index();
                    (ix != e && s.contains(ix) && !visited.contains(ix)).then_some(ix)
                })
            });

            let nr = g.neighbors(dst).filter_map(|n| {
                g.find_edge(dst, n).and_then(|f| {
                    let ix = f.index();
                    (ix != e && s.contains(ix) && !visited.contains(ix)).then_some(ix)
                })
            });

            queue.extend(nl);
            queue.extend(nr);
        }
        remainder.difference_with(&visited);
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

pub fn edges_contained_within<'a, N, E, Ty, Ix>(
    g: &'a Graph<N, E, Ty, Ix>,
    s: &'a BTreeSet<NodeIndex<Ix>>,
) -> impl Iterator<Item = EdgeIndex<Ix>> + 'a
where
    Ty: EdgeType,
    Ix: IndexType,
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

pub fn edge_neighbors<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    e: EdgeIndex<Ix>,
) -> impl Iterator<Item = EdgeIndex<Ix>> + '_
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let (src, dst) = g.edge_endpoints(e).unwrap();
    let src_neighbors = g.neighbors(src).map(move |n| g.find_edge(src, n));
    let dst_neighbors = g.neighbors(dst).map(move |n| g.find_edge(dst, n));
    src_neighbors
        .chain(dst_neighbors)
        .filter_map(move |w| w.filter(|i| e != *i))
}

pub fn node_weight_between<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    left: EdgeIndex<Ix>,
    right: EdgeIndex<Ix>,
) -> Option<&N>
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let (lsrc, ldst) = g.edge_endpoints(left)?;
    let (rsrc, rdst) = g.edge_endpoints(right)?;

    if lsrc == rsrc || lsrc == rdst {
        Some(g.node_weight(lsrc)?)
    } else if ldst == rsrc || ldst == rdst {
        Some(g.node_weight(ldst)?)
    } else {
        None
    }
}

pub fn node_between<N, E, Ty, Ix>(
    g: &Graph<N, E, Ty, Ix>,
    left: EdgeIndex<Ix>,
    right: EdgeIndex<Ix>,
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let Some((lsrc, ldst)) = g.edge_endpoints(left) else {
        return false;
    };
    let Some((rsrc, rdst)) = g.edge_endpoints(right) else {
        return false;
    };

    lsrc == rsrc || lsrc == rdst || ldst == rsrc || ldst == rdst
}
