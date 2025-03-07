use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::loader::parse_molfile_str;
use crate::assembly::{
    index_search
};
use crate::assembly::Bound as AssemblyBound;

// This needs to be combined with the Bounds Enum in main but I'm not sure the 
// best way to do that. Could move it to utils
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyBounds {
    Log,
    IntChain,
    VecChain,
}

/// Implements conversion from `&str` to `PyBounds`
impl FromStr for PyBounds {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "log" => Ok(PyBounds::Log),
            "intchain" => Ok(PyBounds::IntChain),
            "vecchain" => Ok(PyBounds::VecChain),
            _ => Err(PyValueError::new_err(format!("Invalid bound: {}", s))),
        }
    }
}

/// Converts a slice of `PyBounds` to a vector of `AssemblyBound`
fn make_boundlist(u: &[PyBounds]) -> Vec<AssemblyBound> {
    let mut boundlist = u
        .iter()
        .flat_map(|b| match b {
            PyBounds::Log => vec![AssemblyBound::Log],
            PyBounds::IntChain => vec![AssemblyBound::IntChain],
            PyBounds::VecChain => vec![
                AssemblyBound::VecChainSimple,
                AssemblyBound::VecChainSmallFrags,
            ],
        })
        .collect::<Vec<_>>();

    boundlist.dedup(); // Ensure no duplicate bounds
    boundlist
}

/// Processes a `HashSet<String>` from Python and converts it into a `Vec<PyBounds>`
/// Raises an error if any string is invalid.
fn process_bound_set(bound_set: HashSet<String>) -> PyResult<Vec<PyBounds>> {
    bound_set.iter().map(|s| s.parse()).collect::<Result<_, _>>() // Try parsing each string
}

/// Computes the assembly index using specified bounds.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `bound_set`: A set of bounds as strings (from Python).
///
/// # Returns
/// - The computed molecular index as a `u32`.
#[pyfunction]
pub fn _compute_index(mol_block: String, bound_set: HashSet<String>) -> PyResult<u32> {
    let mol_result = parse_molfile_str(&mol_block);
    let py_bounds = process_bound_set(bound_set)?;

    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    let (index, _duplicates, _space) = index_search(&mol, &make_boundlist(&py_bounds));

    Ok(index)
}

/// Computes the assembly index with additional details.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `bound_set`: A set of bounds as strings (from Python).
///
/// # Returns
/// - A `HashMap<String, u32>` containing:
///   - `"index"`: The computed molecular index.
///   - `"duplicates"`: Duplicate count.
///   - `"space"`: Space calculation.
#[pyfunction]
pub fn _compute_verbose_index(mol_block: String, bound_set: HashSet<String>) -> PyResult<HashMap<String, u32>> {
    let mol_result = parse_molfile_str(&mol_block);
    let py_bounds = process_bound_set(bound_set)?;

    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert error to PyErr
    };

    let (ix, duplicates, space) = index_search(&mol, &make_boundlist(&py_bounds));

    let mut data = HashMap::new();
    data.insert("index".to_string(), ix);
    data.insert("duplicates".to_string(), duplicates);
    data.insert("space".to_string(), space);

    Ok(data)
}

/// Retrieves molecular information from a given mol block.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
///
/// # Returns
/// - A `String` containing molecular information.
#[pyfunction]
pub fn _molecule_info(mol_block: String) -> PyResult<String> {
    let mol_result = parse_molfile_str(&mol_block);

    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert error to PyErr
    };

    Ok(mol.info()) // Retrieve molecular info
}

// Registers the Rust functions as a Python module.
//
// This function must match the `lib.name` setting in `Cargo.toml`,
// otherwise, Python will not be able to import the module.
#[pymodule]
fn _pyorca(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_compute_index, m)?)?;
    m.add_function(wrap_pyfunction!(_compute_verbose_index, m)?)?;
    m.add_function(wrap_pyfunction!(_molecule_info, m)?)?;
    Ok(())
}
