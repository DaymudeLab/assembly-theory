//! Expose public `assembly_theory` functionality to a Python package using
//! [`pyo3`](https://docs.rs/pyo3/latest/pyo3/).
//!
//! The package is available [on PyPI](https://pypi.org/project/assembly-theory/);
//! see that README for installation and usage instructions. To build the
//! Python package directly from this crate's source code, see the instructions
//! in the [GitHub README](https://github.com/DaymudeLab/assembly-theory).
//!
//! Note that all Rust functions in this module have the form `_fn_name`, which
//! correspond to the actual Rust function `fn_name` elsewhere in the crate and
//! are exposed to the Python package as `fn_name`.
//!
//! # Python Example
//!
//! ```custom,{class=language-python}
//! import assembly_theory as at
//! 
//! # Load a mol block from file.
//! with open('data/checks/anthracene.mol') as f:
//!     mol_block = f.read()
//! 
//! # Calculate the molecule's assembly index.
//! at.index(mol_block)  # 6
//! ```

use std::{collections::HashSet, str::FromStr};

use pyo3::{
    exceptions::{PyOSError, PyValueError},
    prelude::*,
    PyErr,
};

use crate::{
    assembly::{index, index_search, ParallelMode},
    bounds::Bound as OurBound,
    canonize::CanonizeMode,
    enumerate::EnumerateMode,
    kernels::KernelMode,
    loader::{parse_molfile_str, ParserError},
};

/// Implement a Python version of [`crate::loader::ParserError`].
impl From<ParserError> for PyErr {
    fn from(err: ParserError) -> PyErr {
        PyOSError::new_err(err.to_string())
    }
}

// TODO: Is there a clean way of avoiding the duplication of all our various
// algorithm variant enums?

/// Mirrors the `enumerate::EnumerateMode` enum.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
enum PyEnumerateMode {
    Extend,
    GrowErode,
}

/// Mirrors the `canonize::CanonizeMode` enum.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
enum PyCanonizeMode {
    Nauty,
    Faulon,
    TreeNauty,
    TreeFaulon,
}

/// Mirrors the `assembly::ParallelMode` enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyParallelMode {
    None,
    DepthOne,
    Always,
}

/// Mirrors the `kernels::KernelMode` enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyKernelMode {
    None,
    Once,
    DepthOne,
    Always,
}

/// Mirrors the `bounds::Bound` enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyBound {
    Log,
    Int,
    VecSimple,
    VecSmallFrags,
    CoverSort,
    CoverNoSort,
    CliqueBudget,
}

/// Converts bound options in `&str` format to `PyEnumerateMode`.
impl FromStr for PyEnumerateMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "extend" => Ok(PyEnumerateMode::Extend),
            "grow-erode" => Ok(PyEnumerateMode::GrowErode),
            _ => Err(PyValueError::new_err(format!("Invalid enumerate: {s}"))),
        }
    }
}

/// Converts bound options in `&str` format to `PyCanonizeMode`.
impl FromStr for PyCanonizeMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "nauty" => Ok(PyCanonizeMode::Nauty),
            "faulon" => Ok(PyCanonizeMode::Faulon),
            "tree-nauty" => Ok(PyCanonizeMode::TreeNauty),
            "tree-faulon" => Ok(PyCanonizeMode::TreeFaulon),
            _ => Err(PyValueError::new_err(format!("Invalid canonize: {s}"))),
        }
    }
}

/// Converts bound options in `&str` format to `PyParallelMode`.
impl FromStr for PyParallelMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" => Ok(PyParallelMode::None),
            "depth-one" => Ok(PyParallelMode::DepthOne),
            "always" => Ok(PyParallelMode::Always),
            _ => Err(PyValueError::new_err(format!("Invalid parallel: {s}"))),
        }
    }
}

/// Converts bound options in `&str` format to `PyKernelMode`.
impl FromStr for PyKernelMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" => Ok(PyKernelMode::None),
            "once" => Ok(PyKernelMode::Once),
            "depth-one" => Ok(PyKernelMode::DepthOne),
            "always" => Ok(PyKernelMode::Always),
            _ => Err(PyValueError::new_err(format!("Invalid kernel: {s}"))),
        }
    }
}

/// Converts bound options in `&str` format to `PyBound`.
impl FromStr for PyBound {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "log" => Ok(PyBound::Log),
            "int" => Ok(PyBound::Int),
            "vec-simple" => Ok(PyBound::VecSimple),
            "vec-small-frags" => Ok(PyBound::VecSmallFrags),
            "cover-sort" => Ok(PyBound::CoverSort),
            "cover-no-sort" => Ok(PyBound::CoverNoSort),
            "clique-budget" => Ok(PyBound::CliqueBudget),
            _ => Err(PyValueError::new_err(format!("Invalid bound: {s}"))),
        }
    }
}

/// Converts a `HashSet<String>` of bound strings from Python into a
/// `Vec<PyBound>`, raising an error if any bound string is invalid.
fn process_bound_strs(bound_strs: HashSet<String>) -> PyResult<Vec<PyBound>> {
    bound_strs
        .iter()
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
}

/// Converts a slice of `PyBound`s into a vector of `bounds::Bound`s.
fn make_boundlist(pybounds: &[PyBound]) -> Vec<OurBound> {
    let mut boundlist = pybounds
        .iter()
        .flat_map(|b| match b {
            PyBound::Log => vec![OurBound::Log],
            PyBound::Int => vec![OurBound::Int],
            PyBound::VecSimple => vec![OurBound::VecSimple],
            PyBound::VecSmallFrags => vec![OurBound::VecSmallFrags],
            PyBound::CoverSort => vec![OurBound::CoverSort],
            PyBound::CoverNoSort => vec![OurBound::CoverNoSort],
            PyBound::CliqueBudget => vec![OurBound::CliqueBudget],
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

/// Get a pretty-printable string of this molecule's graph representation.
///
/// Python version of [`crate::molecule::Molecule::info`].
///
/// # Python Parameters
/// - `mol_block`: The contents of a `.mol` file as a `str`.
///
/// # Python Returns
/// - A pretty-printable `str` detailing the molecule's atoms and bonds.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
/// 
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
/// 
/// # Print the molecule's graph structure.
/// print(at.mol_info(mol_block))
///
/// # graph {
/// #     0 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     1 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     2 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     ...
/// #     0 -- 1 [ label = "Double" ]
/// #     1 -- 2 [ label = "Single" ]
/// #     2 -- 5 [ label = "Double" ]
/// #     ...
/// # }
/// ```
#[pyfunction(name = "mol_info")]
pub fn _mol_info(mol_block: String) -> PyResult<String> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Return molecule info.
    Ok(mol.info())
}

/// Computes a molecule's assembly index using an efficient default strategy.
///
/// Python version of [`index`].
///
/// # Python Parameters
/// - `mol_block`: The contents of a `.mol` file as a `str`.
///
/// # Python Returns
/// - The molecule's `int` assembly index.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
/// 
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
/// 
/// # Calculate the molecule's assembly index.
/// at.index(mol_block)  # 6
/// ```
#[pyfunction(name = "index")]
pub fn _index(mol_block: String) -> PyResult<u32> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Calculate the assembly index.
    Ok(index(&mol))
}

/// Computes a molecule's assembly index and related information using a
/// top-down recursive search, parameterized by the specified options.
///
/// Python version of [`index_search`].
///
/// # Python Parameters
///
/// - `mol_block`: The contents of a `.mol` file as a `str`.
/// - `enumerate_str`: An enumeration mode from [`"extend"`, `"grow-erode"`].
/// See [`EnumerateMode`] for details.
/// - `canonize_str`: A canonization mode from [`"nauty"`, `"faulon"`,
/// `"tree-nauty"`, `"tree-faulon"`]. See [`CanonizeMode`] for details.
/// - `parallel_str`: A parallelization mode from [`"none"`, `"depth-one"`,
/// `"always"`]. See [`ParallelMode`] for details.
/// - `kernel_str`: A kernelization mode from [`"none"`, `"once"`,
/// `"depth-one"`, `"always"`]. See [`KernelMode`] for details.
/// - `bound_strs`: A `set` of bounds containing zero or more of [`"log"`,
/// `"int"`, `"vec-simple"`, `"vec-small-frags"`, `"cover-sort"`,
/// `"cover-no-sort"`, `"clique-budget"`]. See [`crate::bounds::Bound`] for
/// details.
/// - `memoize`: `True` iff memoization should be used in search.
///
/// # Python Returns
///
/// A 3-tuple containing:
/// - The molecule's `int` assembly index.
/// - The molecule's `int` number of non-overlapping isomorphic subgraph pairs.
/// - The `int` number of assembly states searched.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
/// 
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
/// 
/// # Calculate the molecule's assembly index using the specified options.
/// (index, num_matches, states_searched) = at.index_search(
///     mol_block,
///     "grow-erode",
///     "nauty",
///     "none",
///     "none",
///     set(["int", "vec-simple", "vec-small-frags"]),
///     False)
///
/// print(f"Assembly Index: {index}")  # 6
/// print(f"Non-Overlapping Isomorphic Subgraph Pairs: {num_matches}")  # 466
/// print(f"Assembly States Searched: {states_searched}")  # 2562
/// ```
#[pyfunction(name = "index_search")]
pub fn _index_search(
    mol_block: String,
    enumerate_str: String,
    canonize_str: String,
    parallel_str: String,
    kernel_str: String,
    bound_strs: HashSet<String>,
    memoize: bool,
) -> PyResult<(u32, u32, usize)> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Parse the various modes and bound options.
    let enumerate_mode = match PyEnumerateMode::from_str(&enumerate_str) {
        Ok(PyEnumerateMode::Extend) => EnumerateMode::Extend,
        Ok(PyEnumerateMode::GrowErode) => EnumerateMode::GrowErode,
        _ => {
            panic!("Unrecognized enumerate mode {enumerate_str}.")
        }
    };
    let canonize_mode = match PyCanonizeMode::from_str(&canonize_str) {
        Ok(PyCanonizeMode::Nauty) => CanonizeMode::Nauty,
        Ok(PyCanonizeMode::Faulon) => CanonizeMode::Faulon,
        Ok(PyCanonizeMode::TreeNauty) => CanonizeMode::TreeNauty,
        Ok(PyCanonizeMode::TreeFaulon) => CanonizeMode::TreeFaulon,
        _ => {
            panic!("Unrecognized canonize mode {canonize_str}.")
        }
    };
    let parallel_mode = match PyParallelMode::from_str(&parallel_str) {
        Ok(PyParallelMode::None) => ParallelMode::None,
        Ok(PyParallelMode::DepthOne) => ParallelMode::DepthOne,
        Ok(PyParallelMode::Always) => ParallelMode::Always,
        _ => {
            panic!("Unrecognized parallel mode {parallel_str}.")
        }
    };
    let kernel_mode = match PyKernelMode::from_str(&kernel_str) {
        Ok(PyKernelMode::None) => KernelMode::None,
        Ok(PyKernelMode::Once) => KernelMode::Once,
        Ok(PyKernelMode::DepthOne) => KernelMode::DepthOne,
        Ok(PyKernelMode::Always) => KernelMode::Always,
        _ => {
            panic!("Unrecognized parallel mode {parallel_str}.")
        }
    };
    let pybounds = process_bound_strs(bound_strs)?;
    let boundlist = make_boundlist(&pybounds);

    // Compute assembly index.
    Ok(index_search(
        &mol,
        enumerate_mode,
        canonize_mode,
        parallel_mode,
        kernel_mode,
        &boundlist,
        memoize,
    ))
}

/// A Python wrapper for the assembly_theory Rust crate.
// Registers the listed functions as a Python module named 'assembly_theory';
// the above line is used as a docstring.
#[pymodule(name = "assembly_theory")]
fn _assembly_theory(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_mol_info, m)?)?;
    m.add_function(wrap_pyfunction!(_index, m)?)?;
    m.add_function(wrap_pyfunction!(_index_search, m)?)?;
    Ok(())
}
