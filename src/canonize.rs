use crate::molecule::{AtomOrBond, CGraph, Molecule};
use lexical_sort::{lexical_cmp, lexical_only_alnum_cmp, natural_cmp, StringSort};
use petgraph::{
    graph::NodeIndex,
    Direction::{Incoming, Outgoing},
    Graph,
};
use std::collections::{HashMap, HashSet, VecDeque};

#[derive(Debug, Clone)]
struct DAGVert {
    atom_idx: NodeIndex,
    inv: u32,
    order: String,
    parents: Vec<NodeIndex>,
    level: u32,
}

impl DAGVert {
    pub fn new(atom_idx: NodeIndex, parents: Vec<NodeIndex>, level: u32) -> Self {
        DAGVert {
            atom_idx,
            inv: 0,
            parents,
            level,
            order: String::new(),
        }
    }
}

#[derive(Debug, Clone)]
struct MolAtomNode {
    color: u32,
    inv: u32,
    order: String,
    num_parents: u32,
}

impl MolAtomNode {
    pub fn new(color: u32, inv: u32, order: String, num_parents: u32) -> Self {
        MolAtomNode {
            color,
            inv,
            order,
            num_parents,
        }
    }
}

// Compute the assembly index of a molecule
pub fn canonize(molecule: &Molecule) -> String {
    let mgraph = molecule.graph();
    let mut mol_graph = CGraph::new_undirected();
    let mut vtx_map = vec![NodeIndex::default(); mgraph.node_count()];

    for bond_idx in mgraph.edge_indices() {
        let bond = mgraph.edge_weight(bond_idx).unwrap();
        let (start_atom_idx, end_atom_idx) = mgraph.edge_endpoints(bond_idx).unwrap();
        let start_atom = mgraph.node_weight(start_atom_idx).unwrap();
        let end_atom = mgraph.node_weight(end_atom_idx).unwrap();

        let new_bond_node_idx = mol_graph.add_node(AtomOrBond::Bond(*bond));

        if *vtx_map.get(start_atom_idx.index()).unwrap() == NodeIndex::default() {
            vtx_map[start_atom_idx.index()] = mol_graph.add_node(AtomOrBond::Atom(*start_atom))
        }

        mol_graph.add_edge(
            *vtx_map.get(start_atom_idx.index()).unwrap(),
            new_bond_node_idx,
            (),
        );

        if *vtx_map.get(end_atom_idx.index()).unwrap() == NodeIndex::default() {
            vtx_map[end_atom_idx.index()] = mol_graph.add_node(AtomOrBond::Atom(*end_atom))
        }

        mol_graph.add_edge(
            *vtx_map.get(end_atom_idx.index()).unwrap(),
            new_bond_node_idx,
            (),
        );
    }

    let mut max_string = String::new();
    for root in mol_graph.node_indices() {
        // for each node in the molecule graph create a signature
        /*
        1. create a dag from each start node
         */
        let mut dag = Graph::<DAGVert, &str>::new();
        let mut dag_vertex_map: HashMap<(NodeIndex, u32), NodeIndex> = HashMap::new();
        let mut mol_g_dag_vertex_map: Vec<Vec<NodeIndex>> = vec![vec![]; mol_graph.node_count()];
        let mut dag_level_list: Vec<Vec<NodeIndex>> = vec![vec![]; mol_graph.node_count()];
        let mut max_level: u32 = 0;

        {
            let mut seen_edges_cache: HashMap<(NodeIndex, NodeIndex), u32> = HashMap::new();
            let mut visited: VecDeque<(NodeIndex, u32)> = VecDeque::new();

            visited.push_back((root, 0));
            seen_edges_cache.insert((root, root), 0);

            let root_vertex_id = dag.add_node(DAGVert::new(root, [].to_vec(), 0));

            dag_vertex_map.insert((root, 0), root_vertex_id);
            dag_level_list[0].push(root_vertex_id);
            mol_g_dag_vertex_map[root.index()].push(root_vertex_id);

            loop {
                let (curr, level) = visited.pop_front().unwrap();

                for neigh in mol_graph.neighbors(curr) {
                    // let mut add_node_to_dag = false;

                    //check if curr -> neigh or neigh -> curr already exists
                    if let Some(seen_at_level) = seen_edges_cache.get(&(curr, neigh)) {
                        // edge already exists at a level above
                        if *seen_at_level < (level + 1) {
                            continue;
                        }
                    }

                    // if add_node_to_dag {

                    //check if a atom has already been processed during this current level's processing
                    if let Some(present_node_idx) = dag_vertex_map.get(&(neigh, (level + 1))) {
                        seen_edges_cache.insert((curr, neigh), level + 1);
                        seen_edges_cache.insert((neigh, curr), level + 1);
                        //get parent node's NodeIndex
                        if let Some(parent_node_idx) = dag_vertex_map.get(&(curr, level)) {
                            dag.add_edge(*parent_node_idx, *present_node_idx, "");
                            // add as parent in the DAGvert
                            dag[*present_node_idx].parents.push(*parent_node_idx);
                        }
                        //skip rest of the processing for the atom
                        continue;
                    }

                    // haven't seen the atom before so add it to dag
                    max_level = level + 1;
                    seen_edges_cache.insert((curr, neigh), level + 1);
                    seen_edges_cache.insert((neigh, curr), level + 1);
                    let child_node_idx = dag.add_node(DAGVert::new(neigh, [].to_vec(), level + 1));
                    dag_vertex_map.insert((neigh, level + 1), child_node_idx);

                    // Overriding the map!!! neigh can be seen before in previous layer
                    mol_g_dag_vertex_map[neigh.index()].push(child_node_idx);

                    // Insert into a level by level hashmap of dag nodes
                    dag_level_list[(level + 1) as usize].push(child_node_idx);

                    visited.push_back((neigh, level + 1));
                    //get parent node's NodeIndex
                    if let Some(parent_node_idx) = dag_vertex_map.get(&(curr, level)) {
                        dag.add_edge(*parent_node_idx, child_node_idx, "");
                        // add as parent in the DAGvert
                        dag[child_node_idx].parents.push(*parent_node_idx);
                    }
                    // }
                }

                if visited.is_empty() {
                    break;
                }
            }
        }

        /*
        2.1. Initialize the molecule graph with color = 0 and invariant no. for each atom from (atom_type,#parents in dag)
        2.2. Do lexicographical ordering of the (atom_type, #parents in dag)
         */
        let mut extended_molg_atom_map: Vec<MolAtomNode> =
            Vec::with_capacity(mol_graph.node_count());
        let mut order_str_set: HashSet<String> = HashSet::new();

        // Each atom does not have just one vertex in dag!!!
        for atom_node in mol_graph.node_indices() {
            // find unique parents for an atom's associated vertices in dag
            let atom_assoc_vert_list = &mol_g_dag_vertex_map[atom_node.index()];
            let mut parents = HashSet::new();
            for vert_id in atom_assoc_vert_list {
                for parent in &dag[*vert_id].parents {
                    parents.insert(parent);
                }
            }
            let parent_len = parents.len();
            let atom_str = mol_graph[atom_node].to_string();
            let atom_order_str = format!("{atom_str}{parent_len}");
            order_str_set.insert(atom_order_str.clone());
            extended_molg_atom_map.insert(
                atom_node.index(),
                MolAtomNode::new(0, 0, atom_order_str, parent_len as u32),
            );
        }

        // lexico-sort
        let mut ordered_vec: Vec<_> = order_str_set.into_iter().collect();
        ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);

        let mut order_idx: HashMap<String, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), (idx as u32) + 1);
        }

        // update the molecule graph invariant based on order idx of lexico-sort of (atom_type,#parents in dag)
        for atom_node in mol_graph.node_indices() {
            extended_molg_atom_map[atom_node.index()].inv = *order_idx
                .get(&extended_molg_atom_map[atom_node.index()].order)
                .unwrap();
        }

        // get the canonized string for current root atom
        let canon_string = canonize_signature(
            &mol_graph,
            &mut dag,
            &mut extended_molg_atom_map,
            &dag_level_list,
            max_level,
            1,
            "".to_string(),
        );

        // lexico-compare strings to save the max one.
        if lexical_cmp(&max_string, &canon_string).is_lt() {
            max_string = canon_string
        }
    }
    max_string
}

fn canonize_signature(
    mol_graph: &CGraph,
    dag: &mut Graph<DAGVert, &str>,
    extended_molg_atom_map: &mut [MolAtomNode],
    dag_level_list: &[Vec<NodeIndex>],
    max_level: u32,
    color_c: u32,
    s_max: String,
) -> String {
    // 1. get the invariants for each atom
    invariant_atom(
        mol_graph,
        dag,
        extended_molg_atom_map,
        dag_level_list,
        max_level,
    );

    // 2. generate orbits based on atom's invariant values
    let mut orbits: HashMap<u32, Vec<NodeIndex>> = HashMap::new();

    for atom in mol_graph.node_indices() {
        // let extended_atom = extended_molg_atom_map.get(&atom).unwrap();
        let extended_atom = &extended_molg_atom_map[atom.index()];
        let atom_inv = extended_atom.inv;
        let parent_len = extended_atom.num_parents;
        // only add atoms which have 2 or more parents in dag
        if parent_len >= 2 {
            orbits
                .entry(atom_inv)
                .and_modify(|atom_list| atom_list.push(atom))
                .or_insert([atom].to_vec());
        }
    }

    // 3. max length of any orbit
    let mut max_orbit_len = 0;
    orbits.values().for_each(|orbit| {
        if orbit.len() > max_orbit_len {
            max_orbit_len = orbit.len()
        }
    });

    if max_orbit_len >= 2 {
        // find the orbits with max len of atoms
        let max_orbits = orbits
            .keys()
            .filter(|orbit| orbits.get(orbit).unwrap().len() == max_orbit_len)
            .collect::<Vec<&u32>>();
        //  if multiple then use orbit with min value
        let min_orbit = (if max_orbits.len() > 1 {
            max_orbits.iter().min()
        } else {
            max_orbits.first()
        })
        .unwrap();

        let mut local_smax = s_max.clone();
        // recurse further for each of the atom in such a orbit and generate a canonized signature by diff. the atoms in same orbit
        for atom in orbits.get(min_orbit).unwrap() {
            extended_molg_atom_map[atom.index()].color = color_c;
            local_smax = canonize_signature(
                mol_graph,
                dag,
                extended_molg_atom_map,
                dag_level_list,
                max_level,
                color_c + 1,
                local_smax,
            );
            extended_molg_atom_map[atom.index()].color = 0;
        }
        local_smax
    } else {
        // no need to recurse further and print the signature-string
        for atom in mol_graph.node_indices() {
            let extended_atom = &extended_molg_atom_map[atom.index()];
            let atom_inv = extended_atom.inv;
            let atom_color = extended_atom.color;
            let parent_len = extended_atom.num_parents;
            // first update any atom without a color to be same as its invariant value
            if (atom_color == 0) && (parent_len >= 2) {
                extended_molg_atom_map[atom.index()].color = atom_inv;
            }
        }
        // start from root node of the dag
        let root_node = dag
            .node_indices()
            .find(|vert| dag.neighbors_directed(*vert, Incoming).count() == 0)
            .unwrap();
        let local_smax = print_signature_string(
            root_node,
            dag,
            mol_graph,
            extended_molg_atom_map,
            &mut vec![],
        );
        if local_smax.len() > s_max.len() {
            local_smax
        } else {
            s_max
        }
    }
}

fn print_signature_string(
    vertex: NodeIndex,
    dag: &Graph<DAGVert, &str>,
    mol_graph: &CGraph,
    extended_molg_atom_map: &[MolAtomNode],
    edges: &mut Vec<(NodeIndex, NodeIndex)>,
) -> String {
    let mut print_sign = String::new();
    print_sign.push('[');
    let atom_idx = dag[vertex].atom_idx;
    let atom = &mol_graph[dag[vertex].atom_idx];
    print_sign.push_str(&atom.to_string());
    let atom_color = extended_molg_atom_map[atom_idx.index()].color;
    if atom_color != 0 {
        print_sign.push(',');
        print_sign.push_str(&atom_color.to_string());
    }
    print_sign.push(']');

    let mut child_vec = dag
        .neighbors_directed(vertex, Outgoing)
        .collect::<Vec<NodeIndex>>();
    if child_vec.is_empty() {
        print_sign
    } else {
        // sort children in descending order of inv
        child_vec.sort_by(|vert_a, vert_b| dag[*vert_b].inv.cmp(&dag[*vert_a].inv));

        let mut sub_print_sign = String::new();

        for child in child_vec {
            if let Some(_edge) = edges
                .iter()
                .find(|egde| (egde.0 == vertex) && (egde.1 == child))
            {
            } else {
                // if the edge is not already seen then add it to seen and generate signature-string for the child
                edges.push((vertex, child));
                sub_print_sign.push_str(&print_signature_string(
                    child,
                    dag,
                    mol_graph,
                    extended_molg_atom_map,
                    edges,
                ));
            }
        }
        if !sub_print_sign.is_empty() {
            print_sign.push('(');
            print_sign.push_str(&sub_print_sign);
            print_sign.push(')');
        }
        print_sign
    }
}

/*
3. Generate Invariant for Atoms
 */
fn invariant_atom(
    mol_graph: &CGraph,
    dag: &mut Graph<DAGVert, &str>,
    extended_molg_atom_map: &mut [MolAtomNode],
    dag_level_list: &[Vec<NodeIndex>],
    max_level: u32,
) {
    let mut count = 0;
    let mut initial = true;
    loop {
        // Unique invariant values
        let start_inv_atoms = HashSet::<u32>::from_iter(
            mol_graph
                .node_indices()
                .map(|atom_idx| extended_molg_atom_map[atom_idx.index()].inv),
        )
        .len();

        /*
        3.1 Generate Invariants for dag vertex
         */

        // first bottom-up
        invariant_dag_vert(
            dag,
            extended_molg_atom_map,
            dag_level_list,
            max_level,
            true,
            initial,
        );

        initial = false;

        // then top-down
        invariant_dag_vert(
            dag,
            extended_molg_atom_map,
            dag_level_list,
            max_level,
            false,
            initial,
        );

        // Create a vector for each atom in molecule graph based on associated vertex in dag
        let mut order_map_vert_atom: Vec<Vec<u32>> =
            vec![vec![0; (max_level + 1).try_into().unwrap()]; mol_graph.node_count()];

        //for reverse sorting use: max_level - dag[vert].level as per paper
        for vert in dag.node_indices() {
            order_map_vert_atom[dag[vert].atom_idx.index()]
                [(max_level - dag[vert].level) as usize] = dag[vert].inv;
        }

        let mut order_to_atom: HashMap<String, Vec<NodeIndex>> = HashMap::new();

        // turn vectors into strings for sorting
        for atom in mol_graph.node_indices() {
            let order_str = order_map_vert_atom[atom.index()]
                .clone()
                .into_iter()
                .map(|i| i.to_string())
                .collect::<String>();
            order_to_atom
                .entry(order_str)
                .and_modify(|atom_list| atom_list.push(atom))
                .or_insert([atom].to_vec());
        }

        // lexico-sort the vectors-strings
        let mut atom_ordered_vec: Vec<_> = order_to_atom.keys().collect();
        atom_ordered_vec.string_sort_unstable(lexical_only_alnum_cmp);
        // atom_ordered_vec.string_sort_unstable(natural_cmp);
        // descend sort
        atom_ordered_vec.reverse();

        // assign the invariant of atom as the order of vectors-strings
        for (idx, order) in atom_ordered_vec.iter().enumerate() {
            for atom in order_to_atom.get(*order).unwrap() {
                // extended_molg_atom_map.entry(*atom).and_modify(|atom_node| atom_node.inv = (idx as u32)+1);
                extended_molg_atom_map[atom.index()].inv = (idx as u32) + 1;
            }
        }

        let end_inv_atoms = HashSet::<u32>::from_iter(
            mol_graph
                .node_indices()
                .map(|atom_idx| extended_molg_atom_map[atom_idx.index()].inv),
        )
        .len();

        // compare the no. of invariants of all the atoms with the one's they started from
        if start_inv_atoms == end_inv_atoms {
            break;
        }

        // Naive way of stopping
        if count > mol_graph.node_count() {
            println!("breaking out because reached upper limit!");
            break;
        }
        count += 1;
    }
}

/*
3. Generate Invariant for Vertices
 */
fn invariant_dag_vert(
    dag: &mut Graph<DAGVert, &str>,
    extended_molg_atom_map: &[MolAtomNode],
    dag_level_list: &[Vec<NodeIndex>],
    max_level: u32,
    bottom: bool,
    initial: bool,
) {
    // top-down or bottom-up calculation of invariants for each vertex in dag
    let mut curr_lvl_range = if bottom { max_level } else { 0 };
    loop {
        // for each vertex generate a invariant-string based on assoc. atom color and atom invariant + directed neighbors
        let mut order_str_set: HashSet<String> = HashSet::new();
        for vert in &dag_level_list[curr_lvl_range as usize] {
            let atom_idx_for_vert = dag[*vert].atom_idx;
            let atom_node = &extended_molg_atom_map[atom_idx_for_vert.index()];
            let (atom_color, atom_inv) = (atom_node.color, atom_node.inv);
            let vert_inv = dag[*vert].inv;
            let mut vert_order;
            let mut child_inv_set: Vec<u32> = Vec::new();
            vert_order = if initial {
                format!("{atom_color}{atom_inv}")
            } else {
                format!("{atom_color}{vert_inv}")
            };
            if bottom {
                // vert_order = format!("{}{}", atom_color, atom_inv);
                for vert_neigh in dag.neighbors_directed(*vert, Outgoing) {
                    child_inv_set.push(dag[vert_neigh].inv);
                }
            } else {
                // vert_order = format!("{}{}", atom_color, vert_inv);
                for vert_neigh in dag.neighbors_directed(*vert, Incoming) {
                    child_inv_set.push(dag[vert_neigh].inv);
                }
            }

            while child_inv_set.len() < 10 {
                child_inv_set.push(0);
            }

            child_inv_set.sort();
            child_inv_set.reverse();
            child_inv_set
                .iter()
                .for_each(|val| vert_order.push_str(&format!("{}", *val)));

            let vec_string = format!("{vert_order:0>20}");
            dag[*vert].order = vec_string.clone();
            order_str_set.insert(vec_string);
        }

        // lexico-sort the invariant-strings in descending order
        let mut ordered_vec: Vec<String> = order_str_set.into_iter().collect();
        ordered_vec.string_sort_unstable(natural_cmp);
        ordered_vec.reverse();

        let mut order_idx: HashMap<String, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), idx as u32);
        }

        // assign the invariant of vertex as the order of invariant-strings
        for vert in &dag_level_list[curr_lvl_range as usize] {
            dag[*vert].inv = (*order_idx.get(&dag[*vert].order).unwrap()) + 1;
        }

        if bottom {
            if curr_lvl_range == 0 {
                break;
            };
            curr_lvl_range -= 1;
        } else {
            if curr_lvl_range == max_level {
                break;
            };
            curr_lvl_range += 1;
        }
    }
}

mod tests {
    #![allow(unused_imports)]
    use super::*;
    use crate::loader;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn canonize_benzene() {
        let path = PathBuf::from(format!("./data/checks/benzene.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let canonical_repr = canonize(&molecule);

        println!("{}", canonical_repr);

        assert_eq!(
            canonical_repr,
            "[C]([2]([C]([1]([C]([2]([C,1])))))[1]([C]([2]([C]([1]([C,1]))))))"
        )
    }

    #[test]
    fn canonize_anthracene() {
        let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let canonical_repr = canonize(&molecule);

        println!("{}", canonical_repr);

        assert_eq!(canonical_repr, "[C]([2]([C]([1]([C]([2]([C]([1]([C]([2]([C]([1]([C]([2]([C,1])))))[1]([C,8]([2]([C]([1]([C,1])))))))))[1]([C,17]([2]([C]([1]([C,8])))))))))[1]([C]([2]([C]([1]([C,17]))))))")
    }

    // Dummy Molecule for testing
    /*
    fn canonize_dummy() {

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
    }
    */
}
