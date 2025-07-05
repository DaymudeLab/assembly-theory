#![allow(dead_code)]

use std::collections::HashSet;

use bit_set::BitSet;
use petgraph::{
    graph::{EdgeIndex, Graph, IndexType},
    Undirected,
};

use crate::utils::{edge_neighbors, node_weight_between};

struct VF3State<'a, N, E1, E2, Ty, Ix, F> {
    pattern: &'a Graph<N, E1, Ty, Ix>,
    target: &'a Graph<N, E2, Ty, Ix>,
    pattern_map: Vec<Option<EdgeIndex<Ix>>>,
    target_map: Vec<Option<EdgeIndex<Ix>>>,
    pattern_depths: Vec<Option<usize>>,
    target_depths: Vec<Option<usize>>,
    depth: usize,
    edge_comparator: F,
}

impl<'a, N, E1, E2, Ix, F> VF3State<'a, N, E1, E2, Undirected, Ix, F>
where
    E1: PartialEq,
    E2: PartialEq,
    N: PartialEq,
    Ix: IndexType,
    F: Fn(&E1, &E2) -> bool,
{
    fn new(
        pattern: &'a Graph<N, E1, Undirected, Ix>,
        target: &'a Graph<N, E2, Undirected, Ix>,
        edge_comparator: F,
    ) -> Self {
        VF3State {
            pattern_map: vec![None; pattern.edge_count()],
            target_map: vec![None; target.edge_count()],
            pattern_depths: vec![None; pattern.edge_count()],
            target_depths: vec![None; target.edge_count()],
            depth: 0,
            pattern,
            target,
            edge_comparator,
        }
    }

    fn is_consistent(&self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) -> bool {
        self.semantic_rule(pattern_edge, target_edge)
            && self.core_rule(pattern_edge, target_edge)
            && self.frontier_rule(pattern_edge, target_edge)
            && self.remainder_rule(pattern_edge, target_edge)
    }

    fn core_rule(&self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) -> bool {
        for neighbor in edge_neighbors(self.pattern, pattern_edge) {
            if let Some(neighbor_in_target) = self.pattern_map[neighbor.index()] {
                if node_weight_between(self.pattern, pattern_edge, neighbor)
                    != node_weight_between(self.target, target_edge, neighbor_in_target)
                {
                    return false;
                }
            }
        }

        for neighbor in edge_neighbors(self.target, target_edge) {
            if let Some(neighbor_in_pattern) = self.target_map[neighbor.index()] {
                if node_weight_between(self.target, target_edge, neighbor)
                    != node_weight_between(self.pattern, pattern_edge, neighbor_in_pattern)
                {
                    return false;
                }
            }
        }

        true
    }

    fn frontier_rule(&self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) -> bool {
        let card_pattern = edge_neighbors(self.pattern, pattern_edge)
            .filter(|e| {
                self.pattern_depths[e.index()].is_some() && self.pattern_map[e.index()].is_none()
            })
            .count();

        let card_target = edge_neighbors(self.target, target_edge)
            .filter(|e| {
                self.target_depths[e.index()].is_some() && self.target_map[e.index()].is_none()
            })
            .count();

        card_target >= card_pattern
    }

    fn remainder_rule(&self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) -> bool {
        let card_pattern = edge_neighbors(self.pattern, pattern_edge)
            .filter(|e| self.pattern_map[e.index()].is_none())
            .count();

        let card_target = edge_neighbors(self.target, target_edge)
            .filter(|e| self.target_map[e.index()].is_none())
            .count();

        card_target >= card_pattern
    }

    fn semantic_rule(&self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) -> bool {
        let edge_match = (self.edge_comparator)(
            self.pattern.edge_weight(pattern_edge).unwrap(),
            self.target.edge_weight(target_edge).unwrap(),
        );

        let (pattern_src, pattern_dst) = self.pattern.edge_endpoints(pattern_edge).unwrap();
        let (target_src, target_dst) = self.target.edge_endpoints(target_edge).unwrap();

        let pattern_src = self.pattern.node_weight(pattern_src);
        let pattern_dst = self.pattern.node_weight(pattern_dst);
        let target_src = self.target.node_weight(target_src);
        let target_dst = self.target.node_weight(target_dst);

        let node_match = (pattern_src == target_src && pattern_dst == target_dst)
            || (pattern_src == target_dst && pattern_dst == target_src);
        edge_match && node_match
    }

    fn pop_mapping(&mut self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) {
        self.pattern_map[pattern_edge.index()] = None;
        self.target_map[target_edge.index()] = None;
        for i in 0..self.pattern_depths.len() {
            if self.pattern_depths[i].is_some_and(|depth| depth >= self.depth) {
                self.pattern_depths[i] = None
            }
        }
        for i in 0..self.target_depths.len() {
            if self.target_depths[i].is_some_and(|depth| depth >= self.depth) {
                self.target_depths[i] = None
            }
        }
        self.depth -= 1;
    }

    fn push_mapping(&mut self, pattern_edge: EdgeIndex<Ix>, target_edge: EdgeIndex<Ix>) {
        self.pattern_map[pattern_edge.index()] = Some(target_edge);
        self.target_map[target_edge.index()] = Some(pattern_edge);
        self.depth += 1;

        if self.pattern_depths[pattern_edge.index()].is_none() {
            self.pattern_depths[pattern_edge.index()] = Some(self.depth);
        }

        if self.target_depths[target_edge.index()].is_none() {
            self.target_depths[target_edge.index()] = Some(self.depth);
        }

        for i in 0..self.pattern_map.len() {
            let neighbors = edge_neighbors(self.pattern, EdgeIndex::new(i)).map(|e| e.index());
            for neighbor in neighbors {
                if self.pattern_map[neighbor].is_none() && self.pattern_depths[neighbor].is_none() {
                    self.pattern_depths[neighbor] = Some(self.depth);
                }
            }
        }

        for i in 0..self.target_map.len() {
            let neighbors = edge_neighbors(self.target, EdgeIndex::new(i)).map(|e| e.index());
            for neighbor in neighbors {
                if self.target_map[neighbor].is_none() && self.target_depths[neighbor].is_none() {
                    self.target_depths[neighbor] = Some(self.depth);
                }
            }
        }
    }

    fn generate_pairs(&mut self) -> Vec<(EdgeIndex<Ix>, EdgeIndex<Ix>)> {
        let mut target_frontier = (0..self.target.edge_count())
            .filter_map(|i| {
                (self.target_map[i].is_none() && self.target_depths[i].is_some())
                    .then_some(EdgeIndex::new(i))
            })
            .peekable();

        let pattern_frontier = (0..self.pattern.edge_count()).filter_map(|i| {
            (self.pattern_map[i].is_none() && self.pattern_depths[i].is_some())
                .then_some(EdgeIndex::new(i))
        });

        if let (Some(u), Some(_)) = (pattern_frontier.min(), target_frontier.peek()) {
            target_frontier.map(|t| (u, t)).collect()
        } else {
            let u = (0..self.pattern.edge_count())
                .find(|i| self.pattern_map[*i].is_none())
                .unwrap();
            (0..self.target.edge_count())
                .filter_map(|i| {
                    self.target_map[i]
                        .is_none()
                        .then_some((EdgeIndex::new(u), EdgeIndex::new(i)))
                })
                .collect()
        }
    }

    pub fn bitset_from_current_mapping(&self) -> BitSet {
        BitSet::from_iter(
            self.target_map
                .iter()
                .enumerate()
                .filter_map(|(ix, e)| e.map(|_| ix)),
        )
    }

    pub fn all_subgraphs(&mut self) -> HashSet<BitSet> {
        let mut isomorphisms = HashSet::new();
        if self.depth == self.pattern.edge_count() {
            isomorphisms.insert(self.bitset_from_current_mapping());
        } else {
            for (pattern_edge, target_edge) in self.generate_pairs() {
                if self.is_consistent(pattern_edge, target_edge) {
                    self.push_mapping(pattern_edge, target_edge);
                    isomorphisms.extend(&mut self.all_subgraphs().into_iter());
                    self.pop_mapping(pattern_edge, target_edge)
                }
            }
        }
        isomorphisms
    }
}

pub fn noninduced_subgraph_isomorphism_iter<N, E1, E2, F, Ix>(
    pattern: &Graph<N, E1, Undirected, Ix>,
    target: &Graph<N, E2, Undirected, Ix>,
    edge_comparator: F,
) -> impl Iterator<Item = BitSet>
where
    N: PartialEq,
    E1: PartialEq,
    E2: PartialEq,
    F: Fn(&E1, &E2) -> bool,
    Ix: IndexType,
{
    let mut state = VF3State::new(pattern, target, edge_comparator);
    state.all_subgraphs().into_iter()
}

mod tests {
    #![allow(unused_imports)]
    use super::*;

    #[test]
    fn p3_is_subgraph_of_c4() {
        let mut p3 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = p3.add_node(0);
        let n1 = p3.add_node(1);
        let n2 = p3.add_node(0);
        p3.add_edge(n0, n1, ());
        p3.add_edge(n1, n2, ());

        let mut c4 = Graph::<u8, (), Undirected>::new_undirected();
        let n0 = c4.add_node(0);
        let n1 = c4.add_node(1);
        let n2 = c4.add_node(0);
        let n3 = c4.add_node(1);
        c4.add_edge(n0, n1, ());
        c4.add_edge(n1, n2, ());
        c4.add_edge(n2, n3, ());
        c4.add_edge(n3, n0, ());

        let count = noninduced_subgraph_isomorphism_iter(&p3, &c4, |e1, e2| e1 == e2)
            .inspect(|g| println!("{g:?}"))
            .count();
        assert_eq!(count, 2)
    }
}
