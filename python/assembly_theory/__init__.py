from .assembly_theory import *

# Currently, maturin and pyo3 are working together to directly expose functions
# in src/python.rs marked with #[pyfunction] to this assembly_theory Python
# package. If ever we need to write Python-specific wrapper code (e.g., some
# preprocessing with RDKit), we can write those extensions here.

# Make assembly_theory docstrings and functions directly accessible using:
# '>>> import assembly_theory'
# instead of requiring:
# '>>> from assembly_theory import assembly_theory'
__doc__ = assembly_theory.__doc__
if hasattr(assembly_theory, "__all__"):
    __all__ = assembly_theory.__all__
