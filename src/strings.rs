//! Implementation of NodeBehavior, EdgeBehavior, and AObject for strings.

use crate::object::{AObject, NodeBehavior, EdgeBehavior};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::Undirected;

// Implement NodeBehavior for char
impl NodeBehavior for char {
    fn kind(&self) -> u8 {
        *self as u8 // Use the ASCII value of the character as its kind
    }
}

// Implement EdgeBehavior for a constant edge label
impl EdgeBehavior for u8 {
    fn kind(&self) -> u8 {
        0
    }
}

pub struct StringGraph {
    string: String,
    graph: Graph<char, u8, Undirected, u32>,
}

impl StringGraph {
    pub fn new(string: String) -> Self {
        let mut graph = Graph::new_undirected();
        let mut prev_node = None;
        let edge_label = 0; // Use a constant edge label

        for ch in string.chars() {
            let current_node = graph.add_node(ch);
            if let Some(prev) = prev_node {
                // Only add an edge if the character is not a space
                if ch != ' ' {
                    graph.add_edge(prev, current_node, edge_label);
                }
            }
            prev_node = if ch == ' ' { None } else { Some(current_node) };
        }

        Self { string, graph }
    }
}

impl AObject for StringGraph {
    type NodeLabel = char;
    type EdgeLabel = u8;

    fn from_graph(g: Graph<Self::NodeLabel, Self::EdgeLabel, Undirected, u32>) -> Self {
        let mut nodes: Vec<(NodeIndex<u32>, char)> = g
            .node_indices()
            .map(|idx| (idx, *g.node_weight(idx).expect("Node weight missing")))
            .collect();

        nodes.sort_by_key(|(idx, _)| idx.index());
        let string: String = nodes.into_iter().map(|(_, c)| c).collect();

        Self { string, graph: g }
    }

    fn graph(&self) -> &Graph<Self::NodeLabel, Self::EdgeLabel, Undirected, u32> {
        &self.graph
    }

    fn info(&self) -> String {
        format!("String as graph: {}", self.string)
    }

    fn is_malformed(&self) -> bool {
        self.string.chars().any(|ch| !ch.is_ascii())
    }

    fn is_basic_unit(&self) -> bool {
        self.string.len() == 1
    }
}
