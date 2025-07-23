//! Parse molecules in the `.mol` file format.
//!
//! # Example
//! ```
//! # use std::{fs, path::PathBuf};
//! use assembly_theory::loader::parse_molfile_str;
//!
//! # fn main() -> Result<(), std::io::Error> {
//! let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
//! let molfile = fs::read_to_string(path)?;
//! let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
//! # Ok(())
//! # }
//! ```

use std::{error::Error, fmt::Display};

use clap::error::Result;

use crate::molecule::{Atom, Bond, Element::Hydrogen, MGraph, Molecule};

/// Thrown by [`parse_molfile_str`] when errors occur.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParserError {
    /// In the counts line, atom count is not an integer value.
    AtomCountNotInt(usize),
    /// In the counts line, bond count is not an integer value.
    BondCountNotInt(usize),
    /// In the counts line, `.mol` file version is not `V2000`.
    FileVersionIsNotV2000(usize),
    /// In an atom line, element symbol is not one of those recognized by
    /// [`Atom`].
    BadElementSymbol(usize, String),
    /// In a bond line, bond number is not an integer value.
    BondNumberNotInt(usize),
    /// In a bond line, bond type is not an integer value.
    BondTypeNotInt(usize),
    /// In a bond line, bond type is not one of those recognized by [`Bond`].
    BadBondType(usize),
    /// The `.mol` file has insufficient lines to reconstruct the molecule.
    NotEnoughLines,
    /// An unknown error that should be reported to the crate maintainers.
    ThisShouldNotHappen,
}

impl Error for ParserError {}

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
            Self::BadElementSymbol(line, sym) => {
                write!(f, "Line {line}: Bad element symbol '{sym}'")
            }
            Self::BondNumberNotInt(line) => {
                write!(f, "Line {line}: Bond number is not an integer")
            }
            Self::BondTypeNotInt(line) => {
                write!(f, "Line {line}: Bond type is not an integer")
            }
            Self::BadBondType(line) => {
                write!(f, "Line {line}: Bond type is not 1, 2, or 3")
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

/// Parse the contents of a `.mol` file as a [`Molecule`].
///
/// If the `.mol` file contents are malformed, a [`ParserError`] is thrown.
///
/// # Example
/// ```
/// # use std::{fs, path::PathBuf};
/// use assembly_theory::loader::parse_molfile_str;
///
/// # fn main() -> Result<(), std::io::Error> {
/// let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
/// let molfile = fs::read_to_string(path)?;
/// let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
/// # Ok(())
/// # }
/// ```
pub fn parse_molfile_str(input: &str) -> Result<Molecule, ParserError> {
    let mut lines = input.lines().enumerate().skip(3); // Skip header block
    let (ix, counts_line) = lines.next().ok_or(ParserError::NotEnoughLines)?;
    let (n_atoms, n_bonds) = parse_counts_line(ix, counts_line)?;

    let mut graph = MGraph::new_undirected();
    let mut atom_indices = Vec::with_capacity(n_atoms); // original atom index -> Option<NodeIndex>

    // Atom parsing with hydrogen exclusion
    lines.by_ref().take(n_atoms).try_for_each(|(i, line)| {
        let atom = parse_atom_line(i, line)?;
        if atom.element() == Hydrogen {
            atom_indices.push(None); // skip H
        } else {
            let idx = graph.add_node(atom);
            atom_indices.push(Some(idx));
        }
        Ok(())
    })?;

    // Bond parsing with skipped H handling
    lines.by_ref().take(n_bonds).try_for_each(|(i, line)| {
        let (first, second, bond) = parse_bond_line(i, line)?;
        let a = atom_indices.get(first - 1).copied().flatten();
        let b = atom_indices.get(second - 1).copied().flatten();
        if let (Some(ai), Some(bi)) = (a, b) {
            graph.add_edge(ai, bi, bond);
        }
        Ok(())
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
        .map_err(|_| ParserError::BondNumberNotInt(line_ix))?;
    let second_atom = bond_line[3..6]
        .trim()
        .parse()
        .map_err(|_| ParserError::BondNumberNotInt(line_ix))?;
    let bond = bond_line[6..9]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParserError::BondTypeNotInt(line_ix))?
        .try_into()
        .map_err(|_| ParserError::BadBondType(line_ix))?;
    Ok((first_atom, second_atom, bond))
}
