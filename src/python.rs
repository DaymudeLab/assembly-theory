use pyo3::prelude::*;
use std::collections::HashMap;

use crate::loader::parse_molfile_str;
use crate::assembly::{index, index_search};



#[pyfunction]
pub fn compute_index(mol_block: String) -> PyResult<u32> {

    let mol_result = parse_molfile_str(&mol_block);
    
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into())
    };

    let ix = index(&mol);

    Ok(ix)
}

#[pyfunction]
pub fn compute_verbose_index(mol_block: String) -> PyResult<HashMap<String, u32>> {

    let mol_result = parse_molfile_str(&mol_block);
    
    let mol = match mol_result {
        Ok(mol) => mol,
        Err(e) => return Err(e.into())
    };

    let (ix, duplicates, space) = index_search(&mol, &[]);

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
