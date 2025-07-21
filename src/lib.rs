//! `assembly_theory` is an open-source, high-performance library for computing
//! *assembly indices* of molecular structures (see, e.g.,
//! [Sharma et al., 2023](https://doi.org/10.1038/s41586-023-06600-9);
//! [Walker et al., 2024](https://doi.org/10.1098/rsif.2024.0367)).
//!
//! This crate is specific to the Rust library; see the
//! [GitHub repository](https://github.com/DaymudeLab/assembly-theory) for ways
//! to use `assembly_theory` as a standalone executable or as a Python library.
//!
//! # Example
//!
//! Install the crate as usual:
//! ```shell
//! cargo add assembly-theory
//! ```
//!
//! Load a molecule from a `.mol` file and calculate its assembly index:
//! ```
//! # use std::{fs, path::PathBuf};
//! use assembly_theory::{
//!     assembly::index,
//!     loader::parse_molfile_str
//! };
//!
//! # fn main() -> Result<(), std::io::Error> {
//! // Load a molecule from a `.mol` file.
//! let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
//! let molfile = fs::read_to_string(path)?;
//! let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");
//!
//! // Compute the molecule's assembly index using an efficient algorithm.
//! assert_eq!(index(&anthracene), 6);
//! # Ok(())
//! # }
//! ```
//!
//! See [`assembly`] for more usage examples.

// TODO: Cite ORCA JOSS paper when it's out.

// Molecule definition, joining operation
pub mod molecule;

// Data IO
pub mod loader;

// The hard bit: compute assembly index
pub mod assembly;

// Bounding strategies for the search phase.
pub mod bounds;

// Algorithms for enumerating connected subgraphs of a molecular graph.
pub mod enumerate;

// Molecule graph canonization algorithms.
pub mod canonize;

// Graph kernelization algorithms.
pub mod kernels;

// Utility functions
mod utils;

// Python library
#[cfg(feature = "python")]
pub mod python;

mod vf3;
