use crate::molecule::{self, Atom, Bond, Element, Index, MGraph, Molecule};
use petgraph::graph::NodeIndex;

use std::{
    fs::{self},
    io::{self, Error},
    path::PathBuf,
};

pub fn parse(p: &PathBuf) -> io::Result<molecule::Molecule> {
    let contents = fs::read_to_string(p).expect("Should have been able to read the file");

    let mut graph: Option<MGraph> = None;
    let mut curr: Vec<String> = Vec::new();
    for line in contents.lines() {
        match line {
            "$$$$" | "M  END" => {
                if !curr.is_empty() {
                    graph = Some(parse_one_molecule(&curr));
                }
                curr.clear();
            }
            &_ => {
                curr.push(line.to_string());
            }
        }
    }

    if let Some(mol_graph) = graph {
        Ok(Molecule::from_graph(mol_graph))
    } else {
        Err(Error::new(
            io::ErrorKind::InvalidData,
            "Something broke while parsing",
        ))
    }
}

pub fn parse_one_molecule(mol_data: &[String]) -> MGraph {
    let mut mol_graph = MGraph::default();

    let (num_atoms, num_bonds) = parse_counts_line(&mol_data[3]);

    let atom_start_l: usize = (4 + num_atoms).try_into().unwrap();
    let bond_start_l: usize = (4 + num_atoms + num_bonds).try_into().unwrap();

    let mut atoms: Vec<&str> = Vec::with_capacity(num_atoms as usize);
    let mut atom_node_ids: Vec<NodeIndex<Index>> = Vec::new();
    // Atoms block parse
    for atom_line in mol_data[4..atom_start_l].iter() {
        let atom = parse_atom_line(atom_line);
        atoms.push(atom);
        match atom {
            // ignore Hydrogen for now
            "H" => atom_node_ids.push(NodeIndex::default()),
            &_ => {
                atom_node_ids.push(mol_graph.add_node(Atom::new(get_element(atom))));
            }
        }
    }

    // Bonds block parse
    for bond_line in mol_data[atom_start_l..bond_start_l].iter() {
        let (atom_one, atom_two, bond_type) = parse_bond_line(bond_line);
        let atom_one_idx = atom_one - 1;
        let atom_two_idx = atom_two - 1;
        // ignore Hydrogen for now
        match (atoms[atom_one_idx as usize] == "H") | (atoms[atom_two_idx as usize] == "H") {
            true => {}
            false => {
                mol_graph.add_edge(
                    atom_node_ids[atom_one_idx as usize],
                    atom_node_ids[atom_two_idx as usize],
                    get_bond(bond_type),
                );
            }
        }
    }
    mol_graph
}

fn parse_counts_line(counts_line: &str) -> (u32, u32) {
    (
        counts_line[0..3].trim().parse().unwrap(),
        counts_line[3..6].trim().parse().unwrap(),
    )
}

fn parse_atom_line(atom_line: &str) -> &str {
    atom_line[31..34].trim()
}

fn parse_bond_line(bond_line: &str) -> (u32, u32, u32) {
    (
        bond_line[0..3].trim().parse().unwrap(),
        bond_line[3..6].trim().parse().unwrap(),
        bond_line[6..9].trim().parse().unwrap(),
    )
}

fn get_element(atom: &str) -> Element {
    match atom {
        "O" => Element::Oxygen,
        "N" => Element::Nitrogen,
        "C" => Element::Carbon,
        &_ => Element::Hydrogen,
    }
}

fn get_bond(bond_type: u32) -> Bond {
    match bond_type {
        1 => Bond::Single,
        2 => Bond::Double,
        _ => Bond::Single,
    }
}
