use pyo3::prelude::*;
use std::collections::HashMap;
use pyo3::exceptions::PyValueError;

use crate::loader::parse_molfile_str;
use crate::assembly::{index_search, log_bound, vec_bound_simple, vec_bound_small_frags, addition_bound};
use crate::assembly::Bound as AssemblyBound;
//use crate::utils::Bounds;


use std::collections::HashSet;
use std::str::FromStr;

#[pyclass]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyBounds {
    Log,
    Addition,
    Vector,
}


impl FromStr for PyBounds {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "log" => Ok(PyBounds::Log),
            "addition" => Ok(PyBounds::Addition),
            "vector" => Ok(PyBounds::Vector),
            _ => Err(PyValueError::new_err(format!("Invalid bound: {}", s)).into()),
        }
    }
}

fn make_boundlist(u: &[PyBounds]) -> Vec<AssemblyBound> {
    let mut boundlist = u
        .iter()
        .flat_map(|b| match b {
            PyBounds::Log => vec![AssemblyBound::Log(log_bound)],
            PyBounds::Addition => vec![AssemblyBound::Addition(addition_bound)],
            PyBounds::Vector => vec![
                AssemblyBound::Vector(vec_bound_simple),
                AssemblyBound::Vector(vec_bound_small_frags),
            ],
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

fn process_bound_set(bound_set: HashSet<String>) -> PyResult<Vec<PyBounds>> {

    let bounds: Vec<PyBounds> = bound_set
        .iter()
        .map(|s| s.parse()) // Try parsing
        .collect::<Result<_, _>>()?; 
    Ok(bounds)
}


#[pyfunction]
pub fn compute_index(mol_block: String, bound_set: HashSet<String>) -> PyResult<u32> {

    let mol_result = parse_molfile_str(&mol_block);
    let py_bounds = process_bound_set(bound_set)?;

    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into())
    };


    let (index, _duplicates, _space) = index_search(&mol, &make_boundlist(&py_bounds));

    Ok(index)
}

#[pyfunction]
pub fn compute_verbose_index(mol_block: String, bound_set: HashSet<String>) -> PyResult<HashMap<String, u32>> {

    let mol_result = parse_molfile_str(&mol_block);
    let py_bounds = process_bound_set(bound_set)?;
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into())
    };

    let (ix, duplicates, space) = index_search(&mol, &make_boundlist(&py_bounds));

    let mut data = HashMap::new();
    data.insert("index".to_string(), ix);
    data.insert("duplicates".to_string(), duplicates);
    data.insert("space".to_string(), space);
    
    Ok(data)
}

#[pyfunction]
pub fn molecule_info(mol_block: String) -> PyResult<String> {

    let mol_result = parse_molfile_str(&mol_block);
    
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into())
    };

    return Ok(mol.info())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn orca(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_index, m)?)?;
    m.add_function(wrap_pyfunction!(compute_verbose_index, m)?)?;
    m.add_function(wrap_pyfunction!(molecule_info, m)?)?;
    Ok(())
}
