use std::collections::{HashMap, VecDeque};

use crate::molecule::Molecule;
use petgraph::{dot::{Config, Dot}, graph::NodeIndex, visit::DfsPostOrder, Directed, Graph};

// Compute the assembly index of a molecule
pub fn canonize(m: &Molecule) -> u32 {
    let mol_graph = m.get_graph();
    for root in mol_graph.node_indices() {
        // for each node in the molecule graph create a signature

        /*
        1. create a DAG from each start node
         */
        let mut DAG = Graph::<String, &str>::new();

        let mut seen_edges_cache: HashMap<(NodeIndex,NodeIndex), u32> = HashMap::new();

        let mut DAG_vertex_map: HashMap<String, NodeIndex> = HashMap::new();

        let mut visited: VecDeque<(NodeIndex,u32)> = VecDeque::new();
        visited.push_back((root, 0));

        seen_edges_cache.insert((root, root), 0);
        DAG_vertex_map.insert(format!("{:?}-{}", root, 0), DAG.add_node(format!("{:?}-{}", root, 0)));
        

        loop {

            let (curr, level) = visited.pop_front().unwrap();

            println!("Current:{:?}", curr);
            
            for neigh in mol_graph.neighbors(curr) {
                println!("Neighbor:{:?}", neigh);
                
                let mut add_node_to_dag = false;
                
                //check if curr -> neigh or neigh -> curr already exists 
                match seen_edges_cache.get(&(curr, neigh)) {
                    Some(seen_at_level) => {
                        // edge already exists at a level above
                        println!("Already Seen: {:?} <--> {:?} at level {}", curr, neigh, *seen_at_level);
                        println!("Current Level: {}", level);
                        if *seen_at_level < (level+1) {
                            continue;
                        }
                        else {
                            //edge at the same level
                            add_node_to_dag = true;
                        }
                    },
                    None => {
                        //No edge found
                        add_node_to_dag = true;
                    }
                }

                if add_node_to_dag {

                    //check if a atom has already been processed during this current level's processing ?
                    match DAG_vertex_map.get(&format!("{:?}-{}", neigh, (level+1))) {
                        Some(present_node_idx) => {
                            seen_edges_cache.insert((curr, neigh), level+1);
                            seen_edges_cache.insert((neigh, curr), level+1);
                            //get parent node's NodeIndex
                            match DAG_vertex_map.get(&format!("{:?}-{}", curr, level)) {
                                Some(parent_node_idx) => {
                                    DAG.add_edge(*parent_node_idx, *present_node_idx, "");
                                }
                                None => {}
                            }
                            //skip rest of the processing for the atom
                            continue;
                        }
                        None => {}
                    }

                    seen_edges_cache.insert((curr, neigh), level+1);
                    seen_edges_cache.insert((neigh, curr), level+1);
                    let child_node_idx = DAG.add_node(format!("{:?}-{}", neigh, level+1));
                    DAG_vertex_map.insert(format!("{:?}-{}", neigh, level+1), child_node_idx);
                    visited.push_back((neigh, level+1));
                    //get parent node's NodeIndex
                    match DAG_vertex_map.get(&format!("{:?}-{}", curr, level)) {
                        Some(parent_node_idx) => {
                            DAG.add_edge(*parent_node_idx, child_node_idx, "");
                        }
                        None => {}
                    }
                }
            }

            if visited.len() == 0 {
                break;
            }
        }
        println!("{:?}", Dot::with_config(&DAG, &[Config::EdgeNoLabel]));
        break;
    }
    0
}
