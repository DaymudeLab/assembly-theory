use std::collections::{btree_set::Range, HashMap, HashSet, VecDeque};
use lexical_sort::{cmp, lexical_cmp, lexical_only_alnum_cmp, natural_cmp, StringSort};

use crate::molecule::{Atom, Bond, Molecule, Element};
use petgraph::{adj::Neighbors, dot::{Config, Dot}, graph::NodeIndex, visit::DfsPostOrder, Directed, Direction::{self, Incoming, Outgoing}, Graph, Undirected};

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
pub fn canonize(m: &Molecule) -> String {
    let mol_graph = m.get_graph();

    /*
    // Dummy Molecule for testing

    let mut mol_graph: Graph<Atom, Bond, Undirected> = Graph::<Atom, Bond, Undirected>::new_undirected();

    let mut vec_nodes: Vec<NodeIndex> = Vec::new();
    vec_nodes.push(NodeIndex::new(99999));
    for i in  0..16 {
        vec_nodes.push(mol_graph.add_node(Atom::new(Element::Carbon))); 
    }
    mol_graph.add_edge(vec_nodes[1],vec_nodes[2],Bond::Single);
    mol_graph.add_edge(vec_nodes[2],vec_nodes[3],Bond::Single);
    mol_graph.add_edge(vec_nodes[3],vec_nodes[4],Bond::Single);
    mol_graph.add_edge(vec_nodes[4],vec_nodes[1],Bond::Single);

    mol_graph.add_edge(vec_nodes[1],vec_nodes[5],Bond::Single);
    mol_graph.add_edge(vec_nodes[2],vec_nodes[7],Bond::Single);
    mol_graph.add_edge(vec_nodes[3],vec_nodes[9],Bond::Single);
    mol_graph.add_edge(vec_nodes[4],vec_nodes[11],Bond::Single);

    mol_graph.add_edge(vec_nodes[13],vec_nodes[14],Bond::Single);
    mol_graph.add_edge(vec_nodes[14],vec_nodes[15],Bond::Single);
    mol_graph.add_edge(vec_nodes[15],vec_nodes[16],Bond::Single);
    mol_graph.add_edge(vec_nodes[16],vec_nodes[13],Bond::Single);

    mol_graph.add_edge(vec_nodes[13],vec_nodes[6],Bond::Single);
    mol_graph.add_edge(vec_nodes[14],vec_nodes[8],Bond::Single);
    mol_graph.add_edge(vec_nodes[15],vec_nodes[10],Bond::Single);
    mol_graph.add_edge(vec_nodes[16],vec_nodes[12],Bond::Single);

    mol_graph.add_edge(vec_nodes[12],vec_nodes[5],Bond::Single);
    mol_graph.add_edge(vec_nodes[6],vec_nodes[7],Bond::Single);
    mol_graph.add_edge(vec_nodes[8],vec_nodes[9],Bond::Single);
    mol_graph.add_edge(vec_nodes[10],vec_nodes[11],Bond::Single);

    mol_graph.add_edge(vec_nodes[12],vec_nodes[11],Bond::Single);
    mol_graph.add_edge(vec_nodes[6],vec_nodes[5],Bond::Single);
    mol_graph.add_edge(vec_nodes[8],vec_nodes[7],Bond::Single);
    mol_graph.add_edge(vec_nodes[10],vec_nodes[9],Bond::Single);
     */

    let mut max_string = String::new();
    for root in mol_graph.node_indices() {
        // let root = vec_nodes[1];
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

            // println!("Current:{:?}", curr);
            
            for neigh in mol_graph.neighbors(curr) {
                // println!("Neighbor:{:?}", neigh);
                
                let mut add_node_to_dag = false;
                
                //check if curr -> neigh or neigh -> curr already exists 
                match seen_edges_cache.get(&(curr, neigh)) {
                    Some(seen_at_level) => {
                        // edge already exists at a level above
                        // println!("Already Seen: {:?} <--> {:?} at level {}", curr, neigh, *seen_at_level);
                        // println!("Current Level: {}", level);
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
        // println!("DAG before invariants");
        // println!("{:?}", Dot::with_config(&DAG, &[Config::EdgeNoLabel]));
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

        // lexico-sort
        let mut ordered_vec: Vec<_> = order_str_set.into_iter().collect();
        ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);

        let mut order_idx:HashMap<String, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), (idx as u32)+1);
        }

        // update the molecule graph invariant based on order idx of lexico-sort of (atom_type,#parents in DAG)
        for atom_node in mol_graph.node_indices() {
            extended_molg_atom_map.entry(atom_node).and_modify(|ext_molg_atom| ext_molg_atom.inv = *order_idx.get(&ext_molg_atom.order).unwrap());
        }

        // get the canonized string for current root atom
        let canon_string = canonize_signature(&mol_graph, &mut DAG, &mut extended_molg_atom_map, &dag_level_list, &mol_g_dag_vertex_map,  max_level, 1, "".to_string());

        // lexico-compare strings to save the max one.
        if lexical_cmp(&max_string, &canon_string).is_lt() {
            max_string = canon_string
        }
    }
    return max_string;
}

fn canonize_signature(
    mol_graph: &Graph<Atom, Bond, Undirected>,
    mut DAG: &mut Graph::<DAGVert, &str>,
    mut extended_molg_atom_map: &mut HashMap<NodeIndex, MolAtomNode>,
    dag_level_list: &HashMap<u32, Vec<NodeIndex>>,
    mol_g_dag_vertex_map: &HashMap<NodeIndex, NodeIndex>,
    max_level: u32,
    color_c: u32,
    s_max: String,
) -> String {
    // 1. get the invariants for each atom
    invariant_atom(&mol_graph, &mut DAG, &mut extended_molg_atom_map, &dag_level_list, max_level);

    // 2. generate orbits based on atom's invariant values
    let mut orbits: HashMap<u32, Vec<NodeIndex>> = HashMap::new();

    for atom in mol_graph.node_indices() {
        let atom_inv = extended_molg_atom_map.get(&atom).unwrap().inv;
        let parent_len = DAG[*mol_g_dag_vertex_map.get(&atom).unwrap()].parents.len() as u32;
        // only add atoms which have 2 or more parents in DAG
        if parent_len >= 2 {
            orbits.entry(atom_inv).and_modify(|atom_list| atom_list.push(atom)).or_insert([atom].to_vec());
        }
    }

    // 3. max length of any orbit
    let mut max_orbit_len = 0;
    orbits.values().for_each(|orbit| if orbit.len() > max_orbit_len {max_orbit_len = orbit.len()});

    if max_orbit_len >= 2 {
        // find the orbits with max len of atoms
        let max_orbits = orbits.keys().filter(|orbit| orbits.get(&orbit).unwrap().len() == max_orbit_len).collect::<Vec<&u32>>();
        //  if multiple then use orbit with min value
        let min_orbit = (if (&max_orbits.len()).clone() > 1 {max_orbits.iter().min()} else {max_orbits.first()}).unwrap();
        
        let mut local_smax = String::new();
        // recurse further for each of the atom in such a orbit and generate a canonized signature by diff. the atoms in same orbit
        for atom in orbits.get(&min_orbit).unwrap() {
            extended_molg_atom_map.entry(*atom).and_modify(|atom_node| atom_node.color = color_c as u32);
            local_smax = canonize_signature(&mol_graph, &mut DAG, &mut extended_molg_atom_map, &dag_level_list, &mol_g_dag_vertex_map,  max_level, color_c+1, s_max.clone());
            extended_molg_atom_map.entry(*atom).and_modify(|atom_node| atom_node.color = 0);
        }
        return local_smax;
    }
    else {
        // not need to recurse further and print the signature-string
        for atom in mol_graph.node_indices() {
            let atom_inv = extended_molg_atom_map.get(&atom).unwrap().inv;
            let atom_color = extended_molg_atom_map.get(&atom).unwrap().color;
            let parent_len = DAG[*mol_g_dag_vertex_map.get(&atom).unwrap()].parents.len() as u32;
            // first update any atom without a color to be same as its invariant value
            if (atom_color == 0) && (parent_len >= 2) {
                extended_molg_atom_map.entry(atom).and_modify(|atom_node| atom_node.color = atom_inv);
            }
        }
        // start from root node of the DAG
        let root_node = DAG.node_indices().find(|vert| DAG.neighbors_directed(*vert, Incoming).count() == 0).unwrap();
        let local_smax = print_signature_string(root_node, &DAG, &mol_graph, &extended_molg_atom_map, &mut vec![]);
        if local_smax.len() > s_max.len() {
            return local_smax;
        }
        else {
            return s_max;    
        }
    }
}

fn print_signature_string(
    vertex: NodeIndex,
    DAG: &Graph::<DAGVert, &str>,
    mol_graph: &Graph<Atom, Bond, Undirected>,
    extended_molg_atom_map: &HashMap<NodeIndex, MolAtomNode>,
    mut edges: &mut Vec<(NodeIndex, NodeIndex)>
) -> String {
    let mut print_sign = String::new();
    print_sign.push('[');
    let atom_idx = DAG[vertex].atom_idx;
    let atom = &mol_graph[DAG[vertex].atom_idx];
    print_sign.push_str(&atom.get_element());
    let atom_color = extended_molg_atom_map.get(&atom_idx).unwrap().color;
    if atom_color != 0 {
        print_sign.push(',');
        print_sign.push_str(&atom_color.to_string());
    }
    print_sign.push(']');

    let mut child_vec = DAG.neighbors_directed(vertex, Outgoing).collect::<Vec<NodeIndex>>();
    if child_vec.len() == 0 { return print_sign; }
    else {
        // sort children in descending order of inv
        // let child_vec = dag_childs;
        child_vec.sort_by(|vert_a, vert_b| DAG[*vert_b].inv.cmp(&DAG[*vert_a].inv));
        
        print_sign.push('(');
        for child in child_vec {
            if let Some(_edge) = edges.iter().find(|egde| (egde.0 == vertex) && (egde.1 == child)) {}
            else {
                // if the edge is not already seen then add it to seen and generate signature-string for the child
                edges.push((vertex, child));
                print_sign.push_str(&print_signature_string(child, &DAG, &mol_graph, &extended_molg_atom_map, edges)); 
            }
        }
        print_sign.push(')');
        return print_sign;
    }
}

/*
    3. Generate Invariant for Atoms
     */
fn invariant_atom(
    mol_graph: &Graph<Atom, Bond, Undirected>,
    mut DAG: &mut Graph::<DAGVert, &str>,
    extended_molg_atom_map: &mut HashMap<NodeIndex, MolAtomNode>,
    dag_level_list: &HashMap<u32, Vec<NodeIndex>>,
    max_level: u32,
) {
    let mut count = 0;
    loop {
        // let start_inv_atoms = HashSet::<u32>::from_iter(mol_graph.node_indices()
        // .into_iter()
        // .map(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv)).len();

        // Better stopping condition ??
        let start_inv_atoms = mol_graph.node_indices()
        .into_iter()
        .filter(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv == 0).count();
        // println!("Begin Atom Invariants: {}", start_inv_atoms);

        /*
        3.1 Generate Invariants for DAG vertex
         */

        // first bottom-up
        invariant_dag_vert(&mut DAG, &extended_molg_atom_map, &dag_level_list, max_level, true);

        // println!("DAG after bottom invariants");
        // println!("{:?}", Dot::with_config(&*DAG, &[Config::EdgeNoLabel]));

        // then top-down
        invariant_dag_vert(&mut DAG, &extended_molg_atom_map, &dag_level_list, max_level, false);

        // println!("DAG after top invariants");
        // println!("{:?}", Dot::with_config(&*DAG, &[Config::EdgeNoLabel]));

        // Create a vector for each atom in molecule graph based on associated vertex in  
        let mut order_map_vert_atom: HashMap<NodeIndex, Vec<u32>> = HashMap::new();
        for atom in mol_graph.node_indices() {
            order_map_vert_atom.insert(atom, vec![0;(max_level+1).try_into().unwrap()]);
        }

        for vert in DAG.node_indices() {
            order_map_vert_atom.entry(DAG[vert].atom_idx).and_modify(|order_vec| order_vec[DAG[vert].level as usize] = DAG[vert].inv);
        }

        let mut order_to_atom: HashMap<String, Vec<NodeIndex>> = HashMap::new();

        // turn vectors into strings for sorting
        for atom in mol_graph.node_indices() {
            let order_str = order_map_vert_atom.get(&atom).unwrap().into_iter().map(|i| i.to_string()).collect::<String>();
            order_to_atom.entry(order_str).and_modify(|atom_list| atom_list.push(atom)).or_insert([atom].to_vec());
        }

        // lexico-sort the vectors-strings
        let mut atom_ordered_vec: Vec<_> = order_to_atom.keys().into_iter().collect();
        atom_ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);

        // assign the invariant of atom as the order of vectors-strings
        for (idx, order) in atom_ordered_vec.iter().enumerate() {
            for atom in order_to_atom.get(*order).unwrap() {
                extended_molg_atom_map.entry(*atom).and_modify(|atom_node| atom_node.inv = (idx as u32)+1);
            }
        }


        // let end_inv_atoms = HashSet::<u32>::from_iter(mol_graph.node_indices()
        // .into_iter()
        // .map(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv)).len();

        // Better stopping condition ??
        let end_inv_atoms = mol_graph.node_indices()
        .into_iter()
        .filter(|atom_idx| extended_molg_atom_map.get(&atom_idx).unwrap().inv == 0).count();


        // println!("End Atom Invariants: {}", end_inv_atoms);

        // compare the no. of invariants of all the atoms with the one's they started from
        if start_inv_atoms == end_inv_atoms {break;}

        // Naive way of stopping
        if count > mol_graph.node_count() {
            println!("breaking out because reached upper limit!");
            break;
        }
        count +=1;
    }
}

/*
    3. Generate Invariant for Vertices
     */
fn invariant_dag_vert(
    DAG: &mut Graph::<DAGVert, &str>,
    extended_molg_atom_map: &HashMap<NodeIndex, MolAtomNode>,
    dag_level_list: &HashMap<u32, Vec<NodeIndex>>,
    max_level: u32,
    bottom: bool) {
        // top-down or bottom-up calculation of invariants for each vertex in DAG
        let mut curr_lvl_range = if bottom {max_level} else {0};
        loop {
            // for each vertex generate a invariant-string based on assoc. atom color and atom invariant + directed neighbors
            let mut order_str_set: HashSet<String> = HashSet::new();
            for vert in dag_level_list.get(&curr_lvl_range).unwrap() {
                let atom_idx_for_vert = DAG[*vert].atom_idx;
                let atom_node = extended_molg_atom_map.get(&atom_idx_for_vert).unwrap();
                let (atom_color, atom_inv) = (atom_node.color, atom_node.inv);
                let mut vert_order = format!("{}{}", atom_color, atom_inv);
                let mut child_inv_set: Vec<u32> = Vec::new();

                if bottom {
                    for vert_neigh in DAG.neighbors_directed(*vert, Outgoing) {
                        child_inv_set.push(DAG[vert_neigh].inv);
                    }
                }
                else {
                    for vert_neigh in DAG.neighbors_directed(*vert, Incoming) {
                        child_inv_set.push(DAG[vert_neigh].inv);
                    }
                }

                child_inv_set.sort();
                child_inv_set.reverse();
                child_inv_set.iter().for_each(|val| vert_order.push_str(&format!("{}", *val)));

                let vec_string = format!("{:0>20}", vert_order);
                DAG[*vert].order = vec_string.clone();
                order_str_set.insert(vec_string);
            }

            // lexico-sort the invariant-strings in descending order
            let mut ordered_vec: Vec<String> = order_str_set.into_iter().collect();
            ordered_vec.string_sort_unstable(natural_cmp);
            ordered_vec.reverse();

            let mut order_idx:HashMap<String, u32> = HashMap::new();

            for (idx, order_str) in ordered_vec.iter().enumerate() {
                order_idx.insert(order_str.clone(), idx as u32);
            }
            
            // assign the invariant of vertex as the order of invariant-strings
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
