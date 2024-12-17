use std::collections::{HashSet, VecDeque};

use petgraph::graph::{Graph, NodeIndex};

pub fn is_subset_connected<N, E>(g: &Graph<N, E>, s: &HashSet<NodeIndex>) -> bool {
    let mut visited = 0;
    let mut queue = VecDeque::from([*s.iter().next().unwrap()]);
    while let Some(v) = queue.pop_front() {
        queue.extend(g.neighbors(v).filter(|n| s.contains(n)));
        visited += 1;
    }
    return visited == g.node_count();
}
