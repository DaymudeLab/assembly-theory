use bit_set::BitSet;

pub struct CompatGraph {
    graph: Vec<BitSet>,
}

impl CompatGraph {
    pub fn new(init_matches: &Vec<(BitSet, BitSet)>) -> Self {
        let size = init_matches.len();

        // Initialize weights and empty graph
        let mut init_graph: Vec<BitSet> = Vec::with_capacity(size);
        for _ in 0..init_matches.len() {
            init_graph.push(BitSet::with_capacity(size));
        }

        // Populate graph
        for (idx1, (h1, h2)) in init_matches.iter().enumerate() {
            for (idx2, (h1p, h2p)) in init_matches[idx1 + 1..].iter().enumerate() {
                let idx2 = idx2 + idx1 + 1;

                let forward_compatible = {
                    h2.is_disjoint(h1p) && 
                    h2.is_disjoint(h2p) &&
                    (h1.is_disjoint(h1p) || h1.is_superset(h1p)) &&
                    (h1.is_disjoint(h2p) || h1.is_superset(h2p))
                };

                if forward_compatible {
                    init_graph[idx1].insert(idx2);
                    init_graph[idx2].insert(idx1);
                }
            }
        }

        Self {
            graph: init_graph,
        }
    }

    pub fn len(&self) -> usize {
        self.graph.len()
    }

    pub fn degree(&self, v: usize, subgraph: &BitSet) -> usize {
        self.graph[v].intersection(subgraph).count()
    }

    pub fn compatible_with(&self, v: usize) -> &BitSet {
        &self.graph[v]
    }

    pub fn neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        
        neighbors
    }

    pub fn forward_neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        let mut to_remove = vec![];
        for u in neighbors.iter() {
            if u <= v {
                to_remove.push(u);
            }
            if u > v {
                break;
            }
        }
        for u in to_remove {
            neighbors.remove(u);
        }

        neighbors
    }

    pub fn are_adjacent(&self, v: usize, u: usize) -> bool {
        self.graph[v].contains(u)
    }
}