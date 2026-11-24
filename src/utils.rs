use std::collections::BTreeSet;

use bit_set::BitSet;
use petgraph::{
    graph::{EdgeIndex, Graph, IndexType},
    EdgeType,
};
use smallvec::SmallVec;

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
    node_set.count()
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
    let mut components = SmallVec::<[BitSet; 8]>::new();

    while let Some(c) = remainder.iter().next() {
        let mut queue = BitSet::with_capacity(s.count());
        let mut visited = BitSet::with_capacity(s.count());
        queue.insert(c);

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
