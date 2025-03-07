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
