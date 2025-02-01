use crate::molecule::{self, Atom, Bond, Element, MGraph, Molecule};
use petgraph::graph::NodeIndex;

use std::{
    fs::{self},
    io::{self, Error},
    num::ParseIntError,
    path::PathBuf,
};

pub fn parse(p: &PathBuf) -> io::Result<molecule::Molecule> {
    let mut curr: Vec<String> = Vec::new();
    let mut graph = None;

    let contents = fs::read_to_string(p)?;

    for line in contents.lines() {
        if let "$$$$" | "M  END" = line {
            if curr.is_empty() {
                curr.clear();
                continue;
            }
            graph = parse_one_molecule(&curr).ok();
        } else {
            curr.push(line.to_string())
        }
    }
    if let Some(g) = graph {
        Ok(Molecule::from_graph(g))
    } else {
        Err(Error::new(
            io::ErrorKind::InvalidData,
            "Something broke while parsing",
        ))
    }
}

pub fn parse_one_molecule(mol_data: &[String]) -> Result<MGraph, ParseIntError> {
    let mut mol_graph = MGraph::default();

    let (num_atoms, num_bonds) = parse_counts_line(&mol_data[3])?;

    let atom_start = 4 + num_atoms;
    let bond_start = 4 + num_atoms + num_bonds;

    let mut atoms = Vec::with_capacity(num_atoms);
    let mut atom_node_ids = Vec::new();
    // Atoms block parse
    for atom_line in mol_data[4..atom_start].iter() {
        let atom = parse_atom_line(atom_line);
        atoms.push(atom);
        if let "H" = atom {
            // ignore Hydrogen for now
            atom_node_ids.push(NodeIndex::default())
        } else {
            atom_node_ids.push(mol_graph.add_node(Atom::new(get_element(atom))));
        }
    }

    // Bonds block parse
    for bond_line in mol_data[atom_start..bond_start].iter() {
        let (atom_one, atom_two, bond_type) = parse_bond_line(bond_line)?;
        let atom_one_idx = atom_one - 1;
        let atom_two_idx = atom_two - 1;
        // ignore Hydrogen for now
        if (atoms[atom_one_idx] == "H") && (atoms[atom_two_idx] == "H") {
            continue;
        }
        mol_graph.add_edge(
            atom_node_ids[atom_one_idx],
            atom_node_ids[atom_two_idx],
            get_bond(bond_type),
        );
    }
    Ok(mol_graph)
}

fn parse_counts_line(counts_line: &str) -> Result<(usize, usize), ParseIntError> {
    Ok((
        counts_line[0..3].trim().parse()?,
        counts_line[3..6].trim().parse()?,
    ))
}

fn parse_atom_line(atom_line: &str) -> &str {
    atom_line[31..34].trim()
}

fn parse_bond_line(bond_line: &str) -> Result<(usize, usize, usize), ParseIntError> {
    Ok((
        bond_line[0..3].trim().parse()?,
        bond_line[3..6].trim().parse()?,
        bond_line[6..9].trim().parse()?,
    ))
}

fn get_element(atom: &str) -> Element {
    match atom {
        "O" => Element::Oxygen,
        "N" => Element::Nitrogen,
        "C" => Element::Carbon,
        &_ => Element::Hydrogen,
    }
}

fn get_bond(bond_type: usize) -> Bond {
    match bond_type {
        1 => Bond::Single,
        2 => Bond::Double,
        _ => Bond::Single,
    }
}
