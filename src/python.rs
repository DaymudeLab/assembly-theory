//! Expose functionality to python library using pyo3.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::assembly::{index_search, serial_index_search};
use crate::bounds::Bound as PyBound;
use crate::loader::parse_molfile_str;

/// Mirrors the `BoundOption` enum in main.rs.
// TODO: Is there a clean way of combining these so we don't have to maintain
// two identical lists? Move to utils?
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyBoundOption {
    Log,
    Int,
    VecSimple,
    VecSmallFrags,
    CoverSort,
    CoverNoSort,
    CliqueBudget,
}

/// Converts bound options in `&str` format to `PyBoundOption`.
impl FromStr for PyBoundOption {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "log" => Ok(PyBoundOption::Log),
            "int" => Ok(PyBoundOption::Int),
            "vecsimple" => Ok(PyBoundOption::VecSimple),
            "vecsmallfrags" => Ok(PyBoundOption::VecSmallFrags),
            "coversort" => Ok(PyBoundOption::CoverSort),
            "covernosort" => Ok(PyBoundOption::CoverNoSort),
            "cliquebudget" => Ok(PyBoundOption::CliqueBudget),
            _ => Err(PyValueError::new_err(format!("Invalid bound: {s}"))),
        }
    }
}

/// Converts a `HashSet<String>` of bound strings from Python into a
/// `Vec<PyBoundOption>`, raising an error if any bound string is invalid.
fn process_bound_strs(bound_strs: HashSet<String>)
    -> PyResult<Vec<PyBoundOption>> {
    bound_strs
        .iter()
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
}

/// Converts a slice of `PyBoundOption`s into a vector of `bounds::Bound`s.
fn make_boundlist(pybounds: &[PyBoundOption]) -> Vec<PyBound> {
    let mut boundlist = pybounds
        .iter()
        .flat_map(|b| match b {
            PyBoundOption::Log => vec![PyBound::Log],
            PyBoundOption::Int => vec![PyBound::Int],
            PyBoundOption::VecSimple => vec![PyBound::VecSimple],
            PyBoundOption::VecSmallFrags => vec![PyBound::VecSmallFrags],
            _ => {
                println!("WARNING: Ignoring bound not implemented yet");
                vec![]
            },
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

/// Computes the molecular assembly index using specified bounds.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `bound_strs`: A set of bounds as strings (from Python).
///
/// # Returns
/// - The computed molecular index as a `u32`.
#[pyfunction]
pub fn _molecular_assembly(
    mol_block: String,
    bound_strs: HashSet<String>,
    serial: bool,
) -> PyResult<u32> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Parse bound options and compute assembly index.
    let pybounds = process_bound_strs(bound_strs)?;
    let boundlist = make_boundlist(&pybounds);
    let (index, _, _) = if serial {
        serial_index_search(&mol, &boundlist)
    } else {
        index_search(&mol, &boundlist)
    };

    Ok(index)
}

/// Computes the molecular assembly index with additional details.
///
/// # Parameters
/// - `mol_block`: The contents of a .mol file as a string.
/// - `bound_strs`: A set of bounds as strings (from Python).
///
/// # Returns
/// - A `HashMap<String, u32>` containing:
///   - `"index"`: The computed molecular index.
///   - `"duplicates"`: Duplicate count.
///   - `"space"`: Space calculation.
#[pyfunction]
pub fn _molecular_assembly_verbose(
    mol_block: String,
    bound_strs: HashSet<String>,
    serial: bool,
) -> PyResult<HashMap<String, usize>> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Parse bound options and compute assembly index.
    let pybounds = process_bound_strs(bound_strs)?;
    let boundlist = make_boundlist(&pybounds);
    let (index, duplicates, space) = if serial {
        serial_index_search(&mol, &boundlist)
    } else {
        index_search(&mol, &boundlist)
    };

    // Package results and return.
    let mut data = HashMap::new();
    data.insert("index".to_string(), index as usize);
    data.insert("duplicates".to_string(), duplicates as usize);
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
    // Parse the .mol file contents as a molecule::Molecule.
    let mol_result = parse_molfile_str(&mol_block);
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into()), // Convert the error to PyErr
    };

    // Return molecule info.
    Ok(mol.info())
}

// Registers the Rust functions as a Python module.
//
// This function must match the `lib.name` setting in `Cargo.toml`,
// otherwise, Python will not be able to import the module.
#[pymodule]
fn _pyat(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_molecular_assembly, m)?)?;
    m.add_function(wrap_pyfunction!(_molecular_assembly_verbose, m)?)?;
    m.add_function(wrap_pyfunction!(_molecule_info, m)?)?;
    Ok(())
}
