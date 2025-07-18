use bit_set::BitSet;

pub struct CompatGraph {
    graph: Vec<BitSet>,
    weights: Vec<usize>,
    matches: Vec<(BitSet, BitSet)>
}

impl CompatGraph {
    pub fn new(init_matches: &Vec<(BitSet, BitSet)>) -> Self{
        let size = init_matches.len();

        // Initialize weights and empty graph
        let mut init_graph: Vec<BitSet> = Vec::with_capacity(size);
        let mut init_weights: Vec<usize> = Vec::with_capacity(size);
        for m in init_matches.iter() {
            init_graph.push(BitSet::with_capacity(size));
            init_weights.push(m.0.len() - 1);
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
            weights: init_weights,
            matches: Vec::new(),
        }
    }

    /*pub fn savings_ground_truth(&self, subgraph: &BitSet) -> usize {
        self.savings_ground_truth_recurse(0, 0, subgraph)
    }

    fn savings_ground_truth_recurse(&self, ix: usize, mut best: usize, subgraph: &BitSet) -> usize {
        if subgraph.len() == 0 {
            return ix;
        }
        let mut cx = ix;

        /*if ix + subgraph.iter().count() <= best && ix + self.remaining_weight_bound(&subgraph) <= best {
            return ix;
        }*/
        /*if ix + self.color_bound(&subgraph) <= best {
            return ix;
        }*/
        /*if ix + self.cover_bound(&subgraph, true) <= best{
            return ix;
        }*/

        // Search for duplicatable fragment
        for v in subgraph.iter() {
            let subgraph_clone = self.forward_neighbors(v, &subgraph);

            cx = cx.max(self.savings_ground_truth_recurse(
                ix + self.weights[v],
                best,
                &subgraph_clone,
            ));
            best = best.max(cx);
        }

        cx
    }*/

    pub fn len(&self) -> usize {
        self.graph.len()
    }

    pub fn degree(&self, v: usize, subgraph: &BitSet) -> usize {
        self.graph[v].intersection(subgraph).count()
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

    pub fn remaining_weight_bound(&self, subgraph: &BitSet) -> usize {
        let deg_sum = subgraph.iter().map(|v| self.degree(v, subgraph)).sum::<usize>() as f32;
        let max_clique = ((1_f32 + (4_f32 * deg_sum + 1_f32).sqrt()) / 2_f32).floor() as usize;
        let mut sum = 0;
        let mut iter = subgraph.iter();
        for _ in 0..max_clique {
            sum += self.weights[iter.next().unwrap()];
        };

        sum
    }

    pub fn color_bound(&self, subgraph: &BitSet) -> usize{
        // Greedy coloring
        let mut colors: Vec<i32> = vec![-1; self.len()];
        let mut num_colors = 0;
        let mut largest: Vec<usize> = Vec::new();
        

        for v in (0..self.matches.len()).rev() {
            if !subgraph.contains(v) {
                continue;
            }

            let mut used: Vec<usize> = vec![0; num_colors];

            for u in subgraph.intersection(&self.graph[v]) {
                if colors[u] != -1 {
                    used[colors[u] as usize] = 1;
                }
            }

            let mut max = 0;
            let mut max_idx = num_colors;
            for i in 0..num_colors {
                if used[i] == 0 && largest[i] > max {
                    max = largest[i];
                    max_idx = i;
                }
            }

            if max_idx == num_colors {
                num_colors += 1;
                largest.push(0);
            }
            if self.weights[v] > largest[max_idx] {
                largest[max_idx] = self.weights[v]
            }

            colors[v] = max_idx as i32;
        }

        largest.iter().sum::<usize>()
    }

    pub fn cover_bound(&self, subgraph: &BitSet, sort: bool) -> usize {
        // Sort vertices
        if sort {
            let mut vertices: Vec<(usize, usize)> = Vec::with_capacity(subgraph.len());
            for v in subgraph {
                vertices.push((v, self.degree(v, subgraph)));
            }   
            vertices.sort_by(|a, b| b.1.cmp(&a.1));
            self.cover_bound_helper(subgraph, vertices.iter().map(|(v, _)| *v))
        }
        else {
            let vertices = (0..self.matches.len()).rev().filter(|v| subgraph.contains(*v));
            self.cover_bound_helper(subgraph, vertices)
        }
    }

    fn cover_bound_helper(&self, subgraph: &BitSet, iter: impl Iterator<Item = usize>) -> usize {
        let mut colors: Vec<Option<Vec<usize>>> = vec![None; self.len()];
        let mut col_weights = vec![];
        let mut num_col = 0;

        for v in iter {
            let mut v_col = Vec::new();
            let mut used = vec![0; num_col];

            // Find colors used in neighborhood of v
            for u in subgraph.intersection(&self.graph[v]) {
                let Some(u_col) = &colors[u] else {
                    continue;
                };

                for c in u_col {
                    used[*c] = 1;
                }
            }

            let mut total_weight = 0;
            let v_val = self.weights[v];
            // Find colors to give to v
            for c in 0..num_col {
                if used[c] == 1 {
                    continue;
                }

                v_col.push(c);
                total_weight += col_weights[c];

                if total_weight >= v_val {
                    break;
                }
            }

            if total_weight == 0 {
                v_col.push(num_col);
                col_weights.push(v_val);
                num_col += 1
            }
            else if total_weight < v_val {
                let mut k = num_col - 1;
                while used[k] == 1 {
                    k -= 1
                }
                col_weights[k] += v_val - total_weight
            }

            colors[v] = Some(v_col);
        };

        col_weights.iter().sum()
    }

    pub fn frag_bound(&self, subgraph: &BitSet, fragments: &[BitSet]) -> usize {
        let total_bonds = fragments.iter().map(|x| x.len()).sum::<usize>();
        let mut bound = 0;
        let sizes = {
            let mut vec = vec![];
            let mut prev = 0;
            for s in subgraph.iter().map(|v| self.weights[v] + 1) {
                if s != prev {
                    vec.push(s);
                    prev = s;
                }
            }

            vec
        };
        
        for i in sizes {
            let mut bound_temp = 0;
            let mut has_bonds = fragments.len();
            let mut num_bonds: Vec<usize> = fragments.iter().map(|x| x.len()).collect();
            let mut smallest_remove = i;

            for v in subgraph.iter() {
                if has_bonds == 0 {
                    break;
                }
                if self.weights[v] + 1 > i {
                    continue;
                }


                let dup = &self.matches[v];
                let bond = dup.1.iter().next().unwrap();
                let mut j = 0;
                while !fragments[j].contains(bond) {
                    j += 1;
                }

                if num_bonds[j] > 0 {
                    let remove = std::cmp::min(dup.0.len(), num_bonds[j]);
                    bound_temp += 1;
                    num_bonds[j] -= remove;
                    smallest_remove = std::cmp::min(smallest_remove, remove);

                    if num_bonds[j] == 0 {
                        has_bonds -= 1;
                    }
                }
            }

            let leftover = num_bonds.iter().sum::<usize>();
            let log = {
                if leftover > 0 {
                    0
                }
                else {
                    (smallest_remove as f32).log2().ceil() as usize
                }
            };
            bound = std::cmp::max(bound, total_bonds - bound_temp - leftover - log);
        }

        bound
    }

}