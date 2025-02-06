

use pyo3::prelude::*;
use crate::loader::parse_mol_block;
use crate::assembly::index;
//use pyo3::types::PyString;

#[pyfunction]
fn assembly_index(mol_block: String) -> PyResult<u32> {

    let mol = parse_mol_block(mol_block)?;
    if mol.is_malformed() {
        panic!("Bad input! Molecule has self-loops or doubled edges")
    }
    let ix = index(&mol);

    Ok(ix)
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
#[pyo3(name = "orca")]
fn _orca(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(assembly_index, m)?)?;

    Ok(())
}
