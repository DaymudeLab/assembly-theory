use petgraph::{
    csr::Neighbors,
    graph::{EdgeIndex, Graph},
    visit::{EdgeCount, IntoEdges, IntoNeighbors},
};

use crate::utils::edge_neighbors;

struct VF3State<N, E> {
    pattern: Graph<N, E>,
    target: Graph<N, E>,
    pattern_map: Vec<Option<EdgeIndex>>,
    target_map: Vec<Option<EdgeIndex>>,
    pattern_depths: Vec<Option<usize>>,
    target_depths: Vec<Option<usize>>,
    depth: usize,
}

impl<N, E> VF3State<N, E> {
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
        todo!()
    }

    fn rule_pred_succ(&self, pattern_edge: EdgeIndex, target_edge: EdgeIndex) -> bool {}

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

    fn search(&mut self) -> Option<Vec<Option<EdgeIndex>>> {
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

    fn iter(&mut self) -> Vec<Vec<Option<EdgeIndex>>> {
        let mut isomorphisms = vec![];
        if self.depth == self.pattern.edge_count() {
            isomorphisms.push(self.pattern_map.clone());
        } else {
            for (pattern_edge, target_edge) in self.generate_pairs() {
                if self.is_consistent(pattern_edge, target_edge) {
                    self.push_mapping(pattern_edge, target_edge);
                    isomorphisms.append(&mut self.iter());
                    self.pop_mapping(pattern_edge, target_edge)
                }
            }
        }
        isomorphisms
    }
}

pub fn noninduced_subgraph_isomorphism_iter<N, E>(pattern: Graph<N, E>, target: Graph<N, E>) {}
