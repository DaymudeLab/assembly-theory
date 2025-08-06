//! `assembly_theory` is an open-source, high-performance library for computing
//! *assembly indices* of molecular structures (see, e.g.,
//! [Sharma et al., 2023](https://doi.org/10.1038/s41586-023-06600-9);
//! [Walker et al., 2024](https://doi.org/10.1098/rsif.2024.0367)).
//!
//! This crate is specific to the Rust library; see the
//! [GitHub repository](https://github.com/DaymudeLab/assembly-theory) for ways
//! to use `assembly_theory` as a standalone executable or as a Python package.
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

// Graph-theoretic utility functions.
mod utils;

// Graph representations of molecules and associated parsing.
pub mod loader;
pub mod molecule;

// Assembly index calculation and supporting functions.
pub mod assembly;
pub mod bounds;
pub mod canonize;
pub mod enumerate;
pub mod kernels;
pub mod memoize;
mod nauty;
mod vf3;

// Python wrapper.
#[cfg(feature = "python")]
pub mod python;
