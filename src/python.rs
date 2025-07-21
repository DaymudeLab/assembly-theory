//! Expose public `assembly_theory` functionality to a Python library using
//! [`pyo3`](https://docs.rs/pyo3/latest/pyo3/).
//!
//! The Python library is available on PyPI as
//! [`assembly-theory`](https://pypi.org/project/assembly-theory/); see that
//! README for installation and usage instructions. To build the Python library
//! directly from this crate's source code, see the instructions in the
//! [GitHub README](https://github.com/DaymudeLab/assembly-theory).

use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
};

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
            "growerode" => Ok(PyEnumerateMode::GrowErode),
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
            "treenauty" => Ok(PyCanonizeMode::TreeNauty),
            "treefaulon" => Ok(PyCanonizeMode::TreeFaulon),
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
            "depthone" => Ok(PyParallelMode::DepthOne),
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
            "depthone" => Ok(PyKernelMode::DepthOne),
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
            "vecsimple" => Ok(PyBound::VecSimple),
            "vecsmallfrags" => Ok(PyBound::VecSmallFrags),
            "coversort" => Ok(PyBound::CoverSort),
            "covernosort" => Ok(PyBound::CoverNoSort),
            "cliquebudget" => Ok(PyBound::CliqueBudget),
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

/// Python version of [`molecule::Molecule::info`].
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
///
/// # Returns
/// - A `String` containing molecular information.
#[pyfunction]
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

/// Python version of [`assembly::index`].
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
///
/// # Returns
/// - The molecule's assembly index as a `u32`.
#[pyfunction]
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

/// Python version of [`assembly::index_search`] returning only the assembly
/// index.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `enumerate_str`: The enumeration mode as a string.
/// - `canonize_str`: The canonization mode as a string.
/// - `parallel_str`: The parallelization mode as a string.
/// - `kernel_str`: The kernelization mode as a string.
/// - `bound_strs`: A set of bounds as strings (from Python).
/// - `memoize`: True iff memoization should be used in search.
///
/// # Returns
/// - The molecule's assembly index as a `u32`.
#[pyfunction]
pub fn _index_search(
    mol_block: String,
    enumerate_str: String,
    canonize_str: String,
    parallel_str: String,
    kernel_str: String,
    bound_strs: HashSet<String>,
    memoize: bool,
) -> PyResult<u32> {
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
    let (index, _, _) = index_search(
        &mol,
        enumerate_mode,
        canonize_mode,
        parallel_mode,
        kernel_mode,
        &boundlist,
        memoize,
    );

    Ok(index)
}

/// Python version of [`assembly::index_search`] returning the assembly index
/// and related information.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `enumerate_str`: The enumeration mode as a string.
/// - `canonize_str`: The canonization mode as a string.
/// - `parallel_str`: The parallelization mode as a string.
/// - `kernel_str`: The kernelization mode as a string.
/// - `bound_strs`: A set of bounds as strings (from Python).
/// - `memoize`: True iff memoization should be used in search.
///
/// # Returns
/// - A `HashMap<String, usize>` containing:
///   - `"index"`: The molecule's assembly index.
///   - `"num_matches"`: The molecule's number of non-overlapping isomorphic
///   subgraph pairs.
///   - `"states_searched"`: The number of assembly states searchede.
#[pyfunction]
pub fn _index_search_verbose(
    mol_block: String,
    enumerate_str: String,
    canonize_str: String,
    parallel_str: String,
    kernel_str: String,
    bound_strs: HashSet<String>,
    memoize: bool,
) -> PyResult<HashMap<String, usize>> {
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
    let (index, num_matches, states_searched) = index_search(
        &mol,
        enumerate_mode,
        canonize_mode,
        parallel_mode,
        kernel_mode,
        &boundlist,
        memoize,
    );

    // Package results and return.
    let mut data = HashMap::new();
    data.insert("index".to_string(), index as usize);
    data.insert("num_matches".to_string(), num_matches as usize);
    data.insert("states_searched".to_string(), states_searched);

    Ok(data)
}

/// A python wrapper for the assembly_theory Rust crate.
// Registers the listed functions as a Python module accessible as _pyat; the
// above line is used as a docstring.
#[pymodule]
fn _pyat(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_mol_info, m)?)?;
    m.add_function(wrap_pyfunction!(_index, m)?)?;
    m.add_function(wrap_pyfunction!(_index_search, m)?)?;
    m.add_function(wrap_pyfunction!(_index_search_verbose, m)?)?;
    Ok(())
}
