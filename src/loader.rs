use crate::molecule::{Atom, Bond, MGraph, Molecule};
use clap::error::Result;
use std::fmt::Display;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParserError {
    AtomCountNotInt(usize),
    BondCountNotInt(usize),
    FileVersionIsNotV2000(usize),
    BadElementSymbol(usize, String),
    BadBondNumber(usize),
    BondTypeNotInt(usize),
    BondTypeOutOfBounds(usize),
    ThisShouldNotHappen,
    NotEnoughLines,
}

pub fn parse_sdfile_str(_input: &str) -> Result<Molecule, ParserError> {
    todo!("SDfile parser unimplemented!")
}

pub fn parse_molfile_str(input: &str) -> Result<Molecule, ParserError> {
    let mut lines = input.lines().enumerate().skip(3); // Skip the header block, 3 lines
    let (ix, counts_line) = lines.next().ok_or_else(|| ParserError::NotEnoughLines)?;
    let (n_atoms, n_bonds) = parse_counts_line(ix, counts_line)?;

    let mut graph = MGraph::new_undirected();
    let mut atom_indices = Vec::new();

    lines
        .by_ref()
        .take(n_atoms)
        .try_fold(&mut graph, |g, (i, l)| {
            parse_atom_line(i, l).map(|atom| {
                atom_indices.push(g.add_node(atom));
                g
            })
        })?;

    lines
        .by_ref()
        .take(n_bonds)
        .try_fold(&mut graph, |g, (i, l)| {
            parse_bond_line(i, l).map(|(first, second, bond)| {
                g.add_edge(atom_indices[first - 1], atom_indices[second - 1], bond);
                g
            })
        })?;

    Ok(Molecule::from_graph(graph))
}

fn parse_counts_line(line_ix: usize, counts_line: &str) -> Result<(usize, usize), ParserError> {
    let n_atoms = counts_line[0..3]
        .trim()
        .parse()
        .map_err(|_| ParserError::AtomCountNotInt(line_ix))?;
    let n_bonds = counts_line[3..6]
        .trim()
        .parse()
        .map_err(|_| ParserError::BondCountNotInt(line_ix))?;
    let version_number = counts_line[33..39].trim().to_uppercase();
    if version_number != "V2000" {
        Err(ParserError::FileVersionIsNotV2000(line_ix))
    } else {
        Ok((n_atoms, n_bonds))
    }
}

fn parse_atom_line(line_ix: usize, atom_line: &str) -> Result<Atom, ParserError> {
    let elem_str = atom_line[31..34].trim();
    let element = elem_str
        .parse()
        .map_err(|_| ParserError::BadElementSymbol(line_ix, elem_str.to_owned()))?;
    let capacity = atom_line[44..47].trim().parse::<u32>().unwrap_or(0);
    Ok(Atom::new(element, capacity))
}

fn parse_bond_line(line_ix: usize, bond_line: &str) -> Result<(usize, usize, Bond), ParserError> {
    let first_atom = bond_line[0..3]
        .trim()
        .parse()
        .map_err(|_| ParserError::BadBondNumber(line_ix))?;
    let second_atom = bond_line[3..6]
        .trim()
        .parse()
        .map_err(|_| ParserError::BadBondNumber(line_ix))?;
    let bond = bond_line[6..9]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParserError::BondTypeNotInt(line_ix))?
        .try_into()
        .map_err(|_| ParserError::BondTypeOutOfBounds(line_ix))?;
    Ok((first_atom, second_atom, bond))
}

impl Display for ParserError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AtomCountNotInt(line) => {
                write!(f, "Line {line}: Atom count is not an integer")
            }
            Self::BondCountNotInt(line) => {
                write!(f, "Line {line}: Bond count is not an integer")
            }
            Self::FileVersionIsNotV2000(line) => {
                write!(f, "Line {line}: File version is not V2000")
            }
            Self::BondTypeNotInt(line) => {
                write!(f, "Line {line}: Bond type is not an integer")
            }
            Self::BondTypeOutOfBounds(line) => {
                write!(f, "Line {line}: Bond type is not 1, 2, or 3")
            }
            Self::BadElementSymbol(line, sym) => {
                write!(f, "Line {line}: Bad element symbol {sym}")
            }
            Self::BadBondNumber(line) => {
                write!(f, "Line {line}: Bad bond number")
            }
            Self::NotEnoughLines => {
                write!(f, "File does not have enough lines")
            }
            Self::ThisShouldNotHappen => {
                write!(f, "This should not happen, report it as a bug")
            }
        }
    }
}
