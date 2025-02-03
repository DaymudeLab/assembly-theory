use std::collections::{btree_set::Range, HashMap, HashSet, VecDeque};
use lexical_sort::{StringSort, lexical_only_alnum_cmp};

use crate::molecule::Molecule;
use petgraph::{adj::Neighbors, dot::{Config, Dot}, graph::NodeIndex, visit::DfsPostOrder, Directed, Direction::{self, Incoming, Outgoing}, Graph};

#[derive(Debug, Clone)]
struct DAGVert {
    atom_idx: NodeIndex,
    inv: u32,
    order: String,
    parents: Vec<NodeIndex>,
    level: u32
}

impl DAGVert {
    pub fn new(atom_idx: NodeIndex, parents: Vec<NodeIndex>, level: u32) -> Self {
        DAGVert {
            atom_idx,
            inv: 0,
            parents,
            level,
            order: String::new()
        }
    }
}

#[derive(Debug, Clone)]
struct MolAtomNode {
    color: u32,
    inv: u32,
    order: String
}

impl MolAtomNode {
    pub fn new(color: u32, inv: u32, order: String) -> Self {
        MolAtomNode {color, inv, order}
    }
}

// Compute the assembly index of a molecule
pub fn canonize(m: &Molecule) -> u32 {
    let mol_graph = m.get_graph().clone();
    for root in mol_graph.node_indices() {
        // for each node in the molecule graph create a signature
        /*
        1. create a DAG from each start node
         */
        let mut DAG = Graph::<DAGVert, &str>::new();

        let mut DAG_vertex_map: HashMap<String, NodeIndex> = HashMap::new();

        let mut mol_g_dag_vertex_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();

        let mut dag_level_list: HashMap<u32, Vec<NodeIndex>> = HashMap::new();

        let mut max_level: u32 = 0;

        {
        let mut seen_edges_cache: HashMap<(NodeIndex,NodeIndex), u32> = HashMap::new();

        let mut visited: VecDeque<(NodeIndex,u32)> = VecDeque::new();
        
        visited.push_back((root, 0));

        seen_edges_cache.insert((root, root), 0);
        // DAG_vertex_map.insert(format!("{:?}-{}", root, 0), DAG.add_node(format!("{:?}-{}", root, 0)));
        let root_vertex_id = DAG.add_node(DAGVert::new(root, [].to_vec(), 0));
        DAG_vertex_map.insert(format!("{:?}-{}", root, 0), root_vertex_id);
        dag_level_list.insert(0, [root_vertex_id].to_vec());
        mol_g_dag_vertex_map.insert(root, root_vertex_id);

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
                                    // add as parent in the DAGvert
                                    (&mut DAG[*present_node_idx]).parents.push(*parent_node_idx);
                                }
                                None => {}
                            }
                            //skip rest of the processing for the atom
                            continue;
                        }
                        None => {}
                    }
                    
                    // haven't seen the atom before so add it to DAG
                    max_level = level + 1;
                    seen_edges_cache.insert((curr, neigh), level+1);
                    seen_edges_cache.insert((neigh, curr), level+1);
                    // let child_node_idx = DAG.add_node(format!("{:?}-{}", neigh, level+1));
                    let child_node_idx = DAG.add_node(DAGVert::new(neigh, [].to_vec(), level+1));
                    DAG_vertex_map.insert(format!("{:?}-{}", neigh, level+1), child_node_idx);
                    mol_g_dag_vertex_map.insert(neigh, child_node_idx);

                    // Insert into a level by level hashmap of dag nodes
                    dag_level_list.entry((level+1)).and_modify(|level_vert_list| level_vert_list.push(child_node_idx)).or_insert([child_node_idx].to_vec());
                    
                    visited.push_back((neigh, level+1));
                    //get parent node's NodeIndex
                    match DAG_vertex_map.get(&format!("{:?}-{}", curr, level)) {
                        Some(parent_node_idx) => {
                            DAG.add_edge(*parent_node_idx, child_node_idx, "");
                            // add as parent in the DAGvert
                            (&mut DAG[child_node_idx]).parents.push(*parent_node_idx);
                        }
                        None => {}
                    }
                }
            }

            if visited.len() == 0 {
                break;
            }
        }
        println!("DAG before invariants");
        println!("{:?}", Dot::with_config(&DAG, &[Config::EdgeNoLabel]));
    }
    
        // println!("DAG Lvl by lvl list");
        // println!("{:?}", dag_level_list);
        // println!("Max lvl: {}",max_level);
        /*
        2.1. Initialize the molecule graph with color = 0 and invariant no. for each atom from (atom_type,#parents in DAG)
        2.2. Do lexicographical ordering of the (atom_type, #parents in DAG)
         */
        let mut extended_molg_atom_map: HashMap<NodeIndex, MolAtomNode> = HashMap::new();
        let mut order_str_set: HashSet<String> = HashSet::new();
        for atom_node in mol_graph.node_indices() {
            let dag_vert = mol_g_dag_vertex_map.get(&atom_node).unwrap();
            let parent_len = DAG[*dag_vert].parents.len();
            let atom_str = mol_graph[atom_node].get_element();
            let atom_order_str = format!("{}{}", atom_str, parent_len);
            order_str_set.insert(atom_order_str.clone());
            extended_molg_atom_map.insert(atom_node, MolAtomNode::new(0, 0, atom_order_str));
        }

        // sort
        let mut ordered_vec: Vec<_> = order_str_set.into_iter().collect();
        ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);

        // println!("{:?}", ordered_vec);

        let mut order_idx:HashMap<String, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), idx as u32);
        }

        for atom_node in mol_graph.node_indices() {
            extended_molg_atom_map.entry(atom_node).and_modify(|ext_molg_atom| ext_molg_atom.inv = *order_idx.get(&ext_molg_atom.order).unwrap());
        }

        // println!("Extended Atom Entry");
        // println!("{:?}", extended_molg_atom_map);

        /*
        3. Generate Invariant for Atoms
         */

        loop {
            // let start_invariants_atoms;
            
            let start_inv_atoms = mol_graph.node_indices()
                .into_iter()
                .map(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv.to_string())
                .collect::<String>();

            println!("Begin Atom Invariants: {}", start_inv_atoms);

            /*
            3.1 Generate Invariants for DAG vertex
             */
            invariant_dag(&mut DAG, &extended_molg_atom_map, &dag_level_list, max_level, true);

            // println!("DAG after bottom invariants");
            // println!("{:?}", Dot::with_config(&DAG, &[Config::EdgeNoLabel]));

            invariant_dag(&mut DAG, &extended_molg_atom_map, &dag_level_list, max_level, false);

            // println!("DAG after top invariants");
            // println!("{:?}", Dot::with_config(&DAG, &[Config::EdgeNoLabel]));

            let mut order_map_vert_atom: HashMap<NodeIndex, Vec<u32>> = HashMap::new();
            for atom in mol_graph.node_indices() {
                order_map_vert_atom.insert(atom, vec![0;(max_level+1).try_into().unwrap()]);
            }

            for vert in DAG.node_indices() {
                order_map_vert_atom.entry(DAG[vert].atom_idx).and_modify(|order_vec| order_vec[DAG[vert].level as usize] = DAG[vert].inv);
            }

            let mut order_to_atom: HashMap<String, Vec<NodeIndex>> = HashMap::new();

            for atom in mol_graph.node_indices() {
                let order_str = order_map_vert_atom.get(&atom).unwrap().into_iter().map(|i| i.to_string()).collect::<String>();
                order_to_atom.entry(order_str).and_modify(|atom_list| atom_list.push(atom)).or_insert([atom].to_vec());
            }

            // sort in descending order
            let mut atom_ordered_vec: Vec<_> = order_to_atom.keys().into_iter().collect();
            atom_ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);
            atom_ordered_vec.reverse();

            for (idx, order) in atom_ordered_vec.iter().enumerate() {
                for atom in order_to_atom.get(*order).unwrap() {
                    extended_molg_atom_map.entry(*atom).and_modify(|atom_node| atom_node.inv = idx as u32);
                }
            }

            let end_inv_atoms = mol_graph.node_indices()
                .into_iter()
                .map(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv.to_string())
                .collect::<String>();

            println!("End Atom Invariants: {}", end_inv_atoms);

            if start_inv_atoms == end_inv_atoms {break}
        }
    break;
    }
    0
}

fn invariant_dag(
    DAG: &mut Graph::<DAGVert, &str>,
    extended_molg_atom_map: &HashMap<NodeIndex, MolAtomNode>,
    dag_level_list: &HashMap<u32, Vec<NodeIndex>>,
    max_level: u32,
    bottom: bool) {
        // let mut vert_order_map:HashMap<NodeIndex, String> = HashMap::new();
        let mut curr_lvl_range = if bottom {max_level} else {0};
        loop {
            let mut order_str_set: HashSet<String> = HashSet::new();
            for vert in dag_level_list.get(&curr_lvl_range).unwrap() {
                let atom_idx_for_vert = DAG[*vert].atom_idx;
                let atom_node = extended_molg_atom_map.get(&atom_idx_for_vert).unwrap();
                let (atom_color, atom_inv) = (atom_node.color, atom_node.inv);
                let mut vert_order = format!("{}{}", atom_color, atom_inv);
                // let vert_neighbor_list
                if bottom {
                    for vert_neigh in DAG.neighbors_directed(*vert, Outgoing) {
                        vert_order.push_str(&format!("{}", DAG[vert_neigh].inv));
                    }
                }
                else {
                    for vert_neigh in DAG.neighbors_directed(*vert, Incoming) {
                        vert_order.push_str(&format!("{}", DAG[vert_neigh].inv));
                    }
                }
                
                DAG[*vert].order = vert_order.clone();
                order_str_set.insert(vert_order);
            }

            let mut ordered_vec: Vec<_> = order_str_set.into_iter().collect();
            ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);
            ordered_vec.reverse();

            println!("{:?}", ordered_vec);

            let mut order_idx:HashMap<String, u32> = HashMap::new();

            for (idx, order_str) in ordered_vec.iter().enumerate() {
                order_idx.insert(order_str.clone(), idx as u32);
            }
    
            for vert in dag_level_list.get(&curr_lvl_range).unwrap() {
                DAG[*vert].inv = (*order_idx.get(&DAG[*vert].order).unwrap())+1;
            }

            if bottom {
                if curr_lvl_range == 0 {break};
                curr_lvl_range -= 1;
            }
            else {
                if curr_lvl_range == max_level {break};
                curr_lvl_range += 1;
            }
        }
}
