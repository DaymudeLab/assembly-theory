//! ORCA (*O*pen *R*eproducible *C*omputation of *A*ssembly Indices) is a free, open-source, and
//! high-performance library that computes the assembly index of a molecule. It is built with
//! modern CPU architectures and parallelism in mind.
//!
//! The crown jewel of the ORCA crate is the [`crate::assembly::index`] function, which computes
//! the assembly index of a given molecule.
//!
//! ORCA comes with batteries included. You can quickly get started with a  parser for the `.mol`
//! file spec ([`crate::loader::parse_molfile_str`]) and an associated graph-theoretic
//! representation of a molecule ([`crate::molecule::Molecule`]), as described in assembly index
//! literature.
//!
//! ORCA can also be used as a python library or as a standalone command-line tool
//!
//! # Example
//! ```
//! # use std::fs;
//! # use std::path::PathBuf;
//! use assembly_theory::*;
//! # fn main() -> Result<(), std::io::Error> {
//! # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
//! // Read a molecule data file
//! let molfile = fs::read_to_string(path)?;
//! let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
//!
//! // Compute assembly index of benzene
//! assert_eq!(assembly::index(&benzene), 3);
//! # Ok(())
//! # }
//! ```

// TODO: Cite ORCA JOSS paper when it's out.

// Molecule definition, joining operation
pub mod molecule;

// Data IO
pub mod loader;

// The hard bit: compute assembly index
pub mod assembly;

// Utility functions
mod utils;

// Python library
#[cfg(feature = "python")]
pub mod python;
