use bit_set::BitSet;
use petgraph::graph::{EdgeIndex, Graph};

use crate::utils::{edge_neighbors, node_between};

struct VF3State<N, E> {
    pattern: Graph<N, E>,
    target: Graph<N, E>,
    pattern_map: Vec<Option<EdgeIndex>>,
    target_map: Vec<Option<EdgeIndex>>,
    pattern_depths: Vec<Option<usize>>,
    target_depths: Vec<Option<usize>>,
    depth: usize,
}

impl<N, E> VF3State<N, E>
where
    E: PartialEq,
    N: PartialEq,
{
    fn new(pattern: Graph<N, E>, target: Graph<N, E>) -> Self {
        VF3State {
            pattern_map: vec![None; pattern.edge_count()],
            target_map: vec![None; target.edge_count()],
            pattern_depths: vec![None; pattern.edge_count()],
            target_depths: vec![None; target.edge_count()],
            depth: 0,
            pattern,
            target,
        }
    }

    fn is_consistent(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {
        self.semantic_rule(pattern_edge, target_edge)
            && self.core_rule(pattern_edge, target_edge)
            && self.frontier_rule(pattern_edge, target_edge)
            && self.remainder_rule(pattern_edge, target_edge)
    }

    fn core_rule(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {
        for neighbor in edge_neighbors(&self.pattern, pattern_edge) {
            let Some(neighbor_in_target) = self.pattern_map[neighbor.index()] else {
                return false;
            };
            if node_between(&self.pattern, pattern_edge, neighbor)
                != node_between(&self.target, target_edge, neighbor_in_target)
            {
                return false;
            }
        }

        for neighbor in edge_neighbors(&self.target, pattern_edge) {
            let Some(neighbor_in_pattern) = self.pattern_map[neighbor.index()] else {
                return false;
            };
            if node_between(&self.target, pattern_edge, neighbor)
                != node_between(&self.pattern, pattern_edge, neighbor_in_pattern)
            {
                return false;
            }
        }

        true
    }

    fn frontier_rule(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {
        let card_pattern = edge_neighbors(&self.pattern, pattern_edge)
            .filter(|e| {
                self.pattern_depths[e.index()].is_some() && self.pattern_map[e.index()].is_none()
            })
            .count();

        let card_target = edge_neighbors(&self.target, target_edge)
            .filter(|e| {
                self.target_depths[e.index()].is_some() && self.target_map[e.index()].is_none()
            })
            .count();

        card_target >= card_pattern
    }

    fn remainder_rule(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {
        let card_pattern = edge_neighbors(&self.pattern, pattern_edge)
            .filter(|e| self.pattern_map[e.index()].is_none())
            .count();

        let card_target = edge_neighbors(&self.target, target_edge)
            .filter(|e| self.target_map[e.index()].is_none())
            .count();

        card_target >= card_pattern
    }

    fn semantic_rule(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {
        let edge_match =
            self.pattern.edge_weight(pattern_edge) == self.target.edge_weight(target_edge);

        let (pattern_src, pattern_dst) = self.pattern.edge_endpoints(pattern_edge).unwrap();
        let (target_src, target_dst) = self.pattern.edge_endpoints(pattern_edge).unwrap();

        let pattern_src = self.pattern.node_weight(pattern_src);
        let pattern_dst = self.pattern.node_weight(pattern_dst);
        let target_src = self.target.node_weight(target_src);
        let target_dst = self.target.node_weight(target_dst);

        let node_match = (pattern_src == target_src && pattern_dst == target_dst)
            || (pattern_src == target_dst && pattern_dst == target_src);
        edge_match && node_match
    }

    fn pop_mapping(&mut self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) {
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

    fn push_mapping(&mut self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) {
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
            let neighbors = edge_neighbors(&self.pattern, EdgeIndex::new(i)).map(|e| e.index());
            for neighbor in neighbors {
                if self.pattern_map[neighbor].is_none() && self.pattern_depths[neighbor].is_none() {
                    self.pattern_depths[neighbor] = Some(self.depth);
                }
            }
        }

        for i in 0..self.target_map.len() {
            let neighbors = edge_neighbors(&self.target, EdgeIndex::new(i)).map(|e| e.index());
            for neighbor in neighbors {
                if self.target_map[neighbor].is_none() && self.target_depths[neighbor].is_none() {
                    self.target_depths[neighbor] = Some(self.depth);
                }
            }
        }
    }

    fn generate_pairs(&mut self) -> Vec<(EdgeIndex, EdgeIndex)> {
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

    pub fn search(&mut self) -> Option<Vec<Option<EdgeIndex>>> {
        if self.depth == self.pattern.edge_count() {
            return Some(self.pattern_map.clone());
        } else {
            for (pattern_edge, target_edge) in self.generate_pairs() {
                if self.is_consistent(pattern_edge, target_edge) {
                    self.push_mapping(pattern_edge, target_edge);
                    if let Some(mapping) = self.search() {
                        return Some(mapping);
                    }
                    self.pop_mapping(pattern_edge, target_edge)
                }
            }
        }
        None
    }

    pub fn bitset_from_current_mapping(&self) -> BitSet {
        BitSet::from_iter(
            self.target_map
                .iter()
                .enumerate()
                .filter_map(|(ix, e)| e.map(|_| ix)),
        )
    }

    pub fn all_subgraphs(&mut self) -> Vec<BitSet> {
        let mut isomorphisms = vec![];
        if self.depth == self.pattern.edge_count() {
            isomorphisms.push(self.bitset_from_current_mapping());
        } else {
            for (pattern_edge, target_edge) in self.generate_pairs() {
                if self.is_consistent(pattern_edge, target_edge) {
                    self.push_mapping(pattern_edge, target_edge);
                    isomorphisms.append(&mut self.all_subgraphs());
                    self.pop_mapping(pattern_edge, target_edge)
                }
            }
        }
        isomorphisms
    }
}

pub fn noninduced_subgraph_isomorphism_iter<N, E>(
    pattern: Graph<N, E>,
    target: Graph<N, E>,
) -> impl Iterator<Item = BitSet>
where
    N: PartialEq,
    E: PartialEq,
{
    let mut state = VF3State::new(pattern, target);
    state.all_subgraphs().into_iter()
}
