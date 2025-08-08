use crate::molecule::{AtomOrBond, CGraph, Molecule};
use bit_set::BitSet;
use petgraph::{
    graph::{EdgeIndex, NodeIndex},
    Direction::{Incoming, Outgoing},
    Graph,
};
use std::collections::{HashMap, HashSet, VecDeque};

// Struct for the vertex of the rooted-DAG, stores associated node's index
// invariant no. and auxiliary information needed while processing DAG.
#[derive(Debug, Clone)]
struct DAGVert {
    atom_idx: NodeIndex,
    inv: u32,
    order: Vec<u32>,
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
            order: vec![],
        }
    }
}

// Struct for the node of the molecule subgraph to store added information required
// by the Faulon et al. (2004) algorithm.
#[derive(Debug, Clone)]
struct MolAtomNode {
    color: u32,
    inv: u32,
    order: Vec<u8>,
    num_parents: u32,
}

impl MolAtomNode {
    pub fn new(color: u32, inv: u32, order: Vec<u8>, num_parents: u32) -> Self {
        MolAtomNode {
            color,
            inv,
            order,
            num_parents,
        }
    }
}

/// Our implementation of [Faulon et al. (2004)](https://doi.org/10.1021/ci0341823).
/// Returns a canonical byte array representation of a molecule's subgraph
/// such that two isomorphic subgraphs will have same representation.
pub fn canonize(molecule: &Molecule, subgraph: &BitSet) -> Option<Vec<u8>> {
    let mgraph = molecule.graph();
    let mut mol_graph = CGraph::new_undirected();
    let mut vtx_map = vec![NodeIndex::default(); mgraph.node_count()];

    // The Faulon et al. (2004) algorithm does not consider bond types while
    // generating the canonization representation. Convert bonds to nodes
    // and construct a new molecule subgraph
    for subgraph_bond_idx in subgraph {
        let bond_idx = EdgeIndex::new(subgraph_bond_idx);
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

    // Maintian lexicographically largest representation for the molecule subgraph
    let mut max_string = Vec::new();

    // Constuct a representation of a molecule subgraph for a starting node
    // and save the lexicographically largest representation by iterating over
    // all the nodes in the molecule subgraph
    for root in mol_graph.node_indices() {
        // Step 1: Create a rooted Directed Acyclic graph (DAG) of the
        // molecule subgraph with the current node as the root
        let mut dag = Graph::<DAGVert, &str>::new();
        // Maintain a molecule node to DAG vertex map
        let mut dag_vertex_map: HashMap<(NodeIndex, u32), NodeIndex> = HashMap::new();
        // Maintain a DAG vertex to molecule node map
        let mut mol_g_dag_vertex_map: Vec<Vec<NodeIndex>> = vec![vec![]; mol_graph.node_count()];
        // Maintain a level by level
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

            // A Breadth-First order traversal to process the molecule subgraph
            // level by level and using references instead of duplicating the
            // previously seen subtrees to construct a rooted-DAG
            loop {
                let (curr, level) = visited.pop_front().unwrap();

                for neigh in mol_graph.neighbors(curr) {
                    // skip further processing if the edge to the neighbor is seen in one of
                    // the previous levels
                    if let Some(seen_at_level) = seen_edges_cache.get(&(curr, neigh)) {
                        if *seen_at_level < (level + 1) {
                            continue;
                        }
                    }

                    // Process only the edge to the neighbor if the neighbor was previously processed at
                    // at the current level
                    if let Some(present_node_idx) = dag_vertex_map.get(&(neigh, (level + 1))) {
                        seen_edges_cache.insert((curr, neigh), level + 1);
                        seen_edges_cache.insert((neigh, curr), level + 1);

                        if let Some(parent_node_idx) = dag_vertex_map.get(&(curr, level)) {
                            dag.add_edge(*parent_node_idx, *present_node_idx, "");
                            dag[*present_node_idx].parents.push(*parent_node_idx);
                        }
                        continue;
                    }

                    // Process the both newly seen neighbor node and the edge to it
                    max_level = level + 1;
                    seen_edges_cache.insert((curr, neigh), level + 1);
                    seen_edges_cache.insert((neigh, curr), level + 1);
                    let child_node_idx = dag.add_node(DAGVert::new(neigh, [].to_vec(), level + 1));
                    dag_vertex_map.insert((neigh, level + 1), child_node_idx);

                    mol_g_dag_vertex_map[neigh.index()].push(child_node_idx);

                    dag_level_list[(level + 1) as usize].push(child_node_idx);

                    visited.push_back((neigh, level + 1));

                    if let Some(parent_node_idx) = dag_vertex_map.get(&(curr, level)) {
                        dag.add_edge(*parent_node_idx, child_node_idx, "");
                        dag[child_node_idx].parents.push(*parent_node_idx);
                    }
                }

                if visited.is_empty() {
                    break;
                }
            }
        }

        // Step 2:
        // First, initialize the molecule subgraph with color set to 0.
        // Next, set the invariant no. for each node to the order index
        // after lexicographical sorting based on the value (node_type, #parents in DAG)
        // associated with each node
        let mut extended_molg_atom_map: Vec<MolAtomNode> =
            Vec::with_capacity(mol_graph.node_count());
        let mut order_str_set: HashSet<Vec<u8>> = HashSet::new();

        // set the value (node_type, #parents in DAG) for each node
        for atom_node in mol_graph.node_indices() {
            let atom_assoc_vert_list = &mol_g_dag_vertex_map[atom_node.index()];
            let mut parents = BitSet::new();
            for vert_id in atom_assoc_vert_list {
                for parent in &dag[*vert_id].parents {
                    parents.insert(parent.index());
                }
            }
            let parent_len = parents.len();

            let mut atom_order_str = mol_graph[atom_node].to_string().into_bytes();
            atom_order_str.extend_from_slice(&parent_len.to_string().into_bytes());
            order_str_set.insert(atom_order_str.clone());
            extended_molg_atom_map.insert(
                atom_node.index(),
                MolAtomNode::new(0, 0, atom_order_str, parent_len as u32),
            );
        }

        // lexicographical sorting based on the value (node_type, #parents in DAG)
        let mut ordered_vec: Vec<_> = order_str_set.into_iter().collect();
        ordered_vec.sort();

        let mut order_idx: HashMap<Vec<u8>, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), (idx as u32) + 1);
        }

        // set the invariant no. for each node to the order index
        for atom_node in mol_graph.node_indices() {
            extended_molg_atom_map[atom_node.index()].inv = *order_idx
                .get(&extended_molg_atom_map[atom_node.index()].order)
                .unwrap();
        }

        // get the canonized representation for current root atom
        let canon_string: Vec<u8> = canonize_signature(
            &mol_graph,
            &mut dag,
            &mut extended_molg_atom_map,
            &dag_level_list,
            max_level,
            1,
            vec![],
        );

        // lexicographical compare the representations to save the larger one
        if max_string < canon_string {
            max_string = canon_string
        }
    }
    Some(max_string)
}

// Generate a canonical representation for the generated rooted-DAG.
fn canonize_signature(
    mol_graph: &CGraph,
    dag: &mut Graph<DAGVert, &str>,
    extended_molg_atom_map: &mut [MolAtomNode],
    dag_level_list: &[Vec<NodeIndex>],
    max_level: u32,
    color_c: u32,
    s_max: Vec<u8>,
) -> Vec<u8> {
    // Step 1: Calculate the invariants no. for each node in the molecule subgraph
    invariant_atom(
        mol_graph,
        dag,
        extended_molg_atom_map,
        dag_level_list,
        max_level,
    );

    // Step 2: Generate orbits based on nodes invariant no.
    // A single orbit is created for each invariant value. Assign a
    // node to an orbit if it has 2 or more parents in the DAG
    let mut orbits: HashMap<u32, Vec<NodeIndex>> = HashMap::new();

    for atom in mol_graph.node_indices() {
        let extended_atom = &extended_molg_atom_map[atom.index()];
        let atom_inv = extended_atom.inv;
        let parent_len = extended_atom.num_parents;

        if parent_len >= 2 {
            orbits
                .entry(atom_inv)
                .and_modify(|atom_list| atom_list.push(atom))
                .or_insert([atom].to_vec());
        }
    }

    let mut max_orbit_len = 0;
    orbits.values().for_each(|orbit| {
        if orbit.len() > max_orbit_len {
            max_orbit_len = orbit.len()
        }
    });

    // Find if any orbit exists with 2 or more nodes
    // then break the tie between these nodes by generating
    // canonized representation but with different node colors
    if max_orbit_len >= 2 {
        // First, find the orbits with max len of atoms
        let max_orbits = orbits
            .keys()
            .filter(|orbit| orbits.get(orbit).unwrap().len() == max_orbit_len)
            .collect::<Vec<&u32>>();
        // If multiple then use orbit with min value
        let min_orbit = (if max_orbits.len() > 1 {
            max_orbits.iter().min()
        } else {
            max_orbits.first()
        })
        .unwrap();

        let mut local_smax = s_max.clone();

        // recurse further for each of the atom in such a orbit and generate a canonized representation
        // by setting a different color for the atom. Use this new canonized representation if it is
        // larger than previously calculated representation.
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
        // Generate the signature repesentation. Use this new canonized representation
        // if it is larger than previously calculated representation.

        // first update any node without a color to be same as its invariant value
        for atom in mol_graph.node_indices() {
            let extended_atom = &extended_molg_atom_map[atom.index()];
            let atom_inv = extended_atom.inv;
            let atom_color = extended_atom.color;
            let parent_len = extended_atom.num_parents;

            if (atom_color == 0) && (parent_len >= 2) {
                extended_molg_atom_map[atom.index()].color = atom_inv;
            }
        }

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

// Constructs the signature representation for a vertex in DAG using the
// Depth-First order traversal on DAG. Called recursively.
fn print_signature_string(
    vertex: NodeIndex,
    dag: &Graph<DAGVert, &str>,
    mol_graph: &CGraph,
    extended_molg_atom_map: &[MolAtomNode],
    edges: &mut Vec<(NodeIndex, NodeIndex)>,
) -> Vec<u8> {
    let mut print_sign: Vec<u8> = vec![];
    print_sign.push(b'[');
    let atom_idx = dag[vertex].atom_idx;
    let atom = &mol_graph[dag[vertex].atom_idx];
    print_sign.extend_from_slice(atom.to_string().as_bytes());
    let atom_color = extended_molg_atom_map[atom_idx.index()].color;
    if atom_color != 0 {
        print_sign.push(b',');
        print_sign.extend_from_slice(&atom_color.to_string().into_bytes());
    }
    print_sign.push(b']');

    let mut child_vec = dag
        .neighbors_directed(vertex, Outgoing)
        .collect::<Vec<NodeIndex>>();
    if child_vec.is_empty() {
        print_sign
    } else {
        // sort children in descending order of invariant no.
        child_vec.sort_by(|vert_a, vert_b| dag[*vert_b].inv.cmp(&dag[*vert_a].inv));

        let mut sub_print_sign = vec![];

        for child in child_vec {
            if let Some(_edge) = edges
                .iter()
                .find(|egde| (egde.0 == vertex) && (egde.1 == child))
            {
            } else {
                // if the edge is not already seen then mark it seen and generate
                // signature representation for the child
                edges.push((vertex, child));
                sub_print_sign.extend_from_slice(&print_signature_string(
                    child,
                    dag,
                    mol_graph,
                    extended_molg_atom_map,
                    edges,
                ));
            }
        }
        if !sub_print_sign.is_empty() {
            print_sign.push(b'(');
            print_sign.extend_from_slice(&sub_print_sign);
            print_sign.push(b')');
        }
        print_sign
    }
}

// Calculate the invariant no. for the nodes in the molecule subgraph
// based on the invariant no. of vertices of the DAG. The process makes
// repeated passes to calculate invariant no. until it stabilizes ie.
// no. of unique invariant values don't change
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
        // Calculate unique invariant values at the start of the pass
        let start_inv_atoms = HashSet::<u32>::from_iter(
            mol_graph
                .node_indices()
                .map(|atom_idx| extended_molg_atom_map[atom_idx.index()].inv),
        )
        .len();

        // Step 1: Calculate the vertex invariants bottom up
        invariant_dag_vert(
            dag,
            extended_molg_atom_map,
            dag_level_list,
            max_level,
            true,
            initial,
        );

        initial = false;

        // Step 2: Calculate the vertex invariants top down
        invariant_dag_vert(
            dag,
            extended_molg_atom_map,
            dag_level_list,
            max_level,
            false,
            initial,
        );

        // Step 3: set the invariant no. for each node to the order index
        // after sorting based on the invariant no. of vertices associated
        // with each node

        // Create a vector to store invariant no. of the associated vertices
        let mut order_map_vert_atom: Vec<Vec<u32>> =
            vec![vec![0; (max_level + 1).try_into().unwrap()]; mol_graph.node_count()];

        for vert in dag.node_indices() {
            order_map_vert_atom[dag[vert].atom_idx.index()]
                [(max_level - dag[vert].level) as usize] = dag[vert].inv;
        }

        let mut order_to_atom: HashMap<Vec<u32>, Vec<NodeIndex>> = HashMap::new();

        for atom in mol_graph.node_indices() {
            let order_str = &order_map_vert_atom[atom.index()];

            order_to_atom
                .entry(order_str.to_vec())
                .and_modify(|atom_list| atom_list.push(atom))
                .or_insert([atom].to_vec());
        }

        // lexicographicaly sort the vectors in descending order
        let mut atom_ordered_vec: Vec<_> = order_to_atom.keys().collect();
        atom_ordered_vec.sort();
        atom_ordered_vec.reverse();

        // assign the invariant of atom as the order index of sorted vectors
        for (idx, order) in atom_ordered_vec.iter().enumerate() {
            for atom in order_to_atom.get(*order).unwrap() {
                extended_molg_atom_map[atom.index()].inv = (idx as u32) + 1;
            }
        }

        // Calculate unique invariant values at the end of the pass
        let end_inv_atoms = HashSet::<u32>::from_iter(
            mol_graph
                .node_indices()
                .map(|atom_idx| extended_molg_atom_map[atom_idx.index()].inv),
        )
        .len();

        // Stop the process if the invariant values stabilize
        if start_inv_atoms == end_inv_atoms {
            break;
        }

        // Hard stopping the process
        if count > mol_graph.node_count() {
            println!("breaking out because reached upper limit!");
            break;
        }
        count += 1;
    }
}

// Calculate the invariant no. for the vertices of DAG based on
// associated node's color and invariant no.
// The invariant calculation can proceed bottom-up or top-down
// specified by bottom flag
fn invariant_dag_vert(
    dag: &mut Graph<DAGVert, &str>,
    extended_molg_atom_map: &[MolAtomNode],
    dag_level_list: &[Vec<NodeIndex>],
    max_level: u32,
    bottom: bool,
    initial: bool,
) {
    let mut curr_lvl_range = if bottom { max_level } else { 0 };
    // for each vertex generate a ordering vector based on associated node's color,
    // invariant no. and vertex's directed neighbor's invariant values
    loop {
        let mut order_str_set: HashSet<Vec<u32>> = HashSet::new();
        for vert in &dag_level_list[curr_lvl_range as usize] {
            let atom_idx_for_vert = dag[*vert].atom_idx;
            let atom_node = &extended_molg_atom_map[atom_idx_for_vert.index()];
            let (atom_color, atom_inv) = (atom_node.color, atom_node.inv);
            let vert_inv = dag[*vert].inv;
            let mut vert_order: Vec<u32> = vec![];
            let mut child_inv_set: Vec<u32> = Vec::new();
            vert_order.push(atom_color);
            if initial {
                vert_order.push(atom_inv);
            } else {
                vert_order.push(vert_inv);
            }

            if bottom {
                for vert_neigh in dag.neighbors_directed(*vert, Outgoing) {
                    child_inv_set.push(dag[vert_neigh].inv);
                }
            } else {
                for vert_neigh in dag.neighbors_directed(*vert, Incoming) {
                    child_inv_set.push(dag[vert_neigh].inv);
                }
            }

            // sort the invariant values of the directed neighbors
            child_inv_set.sort();
            child_inv_set.reverse();
            vert_order.append(&mut child_inv_set);

            dag[*vert].order = vert_order.clone();
            order_str_set.insert(vert_order);
        }

        // lexicographicaly sort the vectors in descending order
        let mut ordered_vec: Vec<Vec<u32>> = order_str_set.into_iter().collect();
        ordered_vec.sort();
        ordered_vec.reverse();

        let mut order_idx: HashMap<Vec<u32>, u32> = HashMap::new();

        for (idx, order_str) in ordered_vec.iter().enumerate() {
            order_idx.insert(order_str.clone(), idx as u32);
        }

        // assign the invariant of the vertex as the order index of sorted vectors
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
        let canonical_repr = canonize(
            &molecule,
            &BitSet::from_iter(molecule.graph().edge_indices().map(|e| e.index())),
        )
        .unwrap();

        assert_eq!(
            String::from_utf8(canonical_repr).unwrap(),
            "[C]([2]([C]([1]([C]([2]([C,1])))))[1]([C]([2]([C]([1]([C,1]))))))"
        )
    }

    #[test]
    fn canonize_anthracene() {
        let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read the data file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
        let canonical_repr = canonize(
            &molecule,
            &BitSet::from_iter(molecule.graph().edge_indices().map(|e| e.index())),
        )
        .unwrap();

        assert_eq!(String::from_utf8(canonical_repr).unwrap(), "[C]([2]([C]([1]([C]([2]([C]([1]([C]([2]([C]([1]([C]([2]([C,1])))))[1]([C,8]([2]([C]([1]([C,1])))))))))[1]([C,17]([2]([C]([1]([C,8])))))))))[1]([C]([2]([C]([1]([C,17]))))))")
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
