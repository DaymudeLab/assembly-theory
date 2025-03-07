# ORCA Python Library
ORCA is a Python library for computing the assembly index of molecules using Rust for high-performance calculations. It integrates with `RDKit` and is built using `maturin`.

To install first create a virtual environment 

`python -m venv orca-env`

Activate the environment:

Windows:
    `orca-env\Scripts\activate`


macOS \& Unix:
    `source tutorial-env/bin/activate`

Install maturin: 

`pip install maturin`

Build the library:

`maturin develop`

# Running Tests
To run the test suite, install pytest and execute the tests from the top-level orca directory:

```
pip install pytest

pytest python/tests
```

# Example usage

ORCA can compute the assembly index of molecules using RDKit's Mol class. Below is a basic example:

```python
import pyorca
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
pyorca.compute_ma(anthracene) # 6
```


# Details

`pyorca` provides three primary functions:

* `compute_ma(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False) -> int`
Computes the assembly index of a given molecule.
* `compute_ma_verbose(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False) -> dict`
Returns additional details, including the number of duplicated isomorphic subgraphs (duplicates) and the size of the search space (space).
* `get_molecule_info(mol: Chem.Mol) -> str`
Provides a string representation of the molecule’s atom and bond structure, useful for debugging.

## Search Strategy Options
Both compute_ma and compute_ma_verbose support optional parameters for controlling the search strategy in the branch-and-bound algorithm:

* `bounds: set[str]` – Specifies the heuristic bounds used to optimize the search.
* * Options: `{"log"}`, `{"intchain"}`, or `{"log", "intchain"}`.
* * Defaults to the best-performing option when not provided.
* `no_bounds: bool` – If `True`, disables all bounds, forcing a complete search of all pathways.

The effect of these on the search space can be see using `compute_ma_verbose` for example:

```python
from rdkit import Chem
import pyorca

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")

pyorca.compute_ma_verbose(mol, bounds={"log"})
# {'index': 6, 'duplicates': 418, 'space': 40507}
pyorca.compute_ma_verbose(mol, bounds={"intchain"})
# {'index': 6, 'duplicates': 418, 'space': 3484}
pyorca.compute_ma_verbose(mol, bounds={"intchain", "log"})
# {'index': 6, 'duplicates': 418, 'space': 3081}
pyorca.compute_ma_verbose(mol, no_bounds=True)
# {'index': 6, 'duplicates': 418, 'space': 129409}

# Invalid combination
pyorca.compute_ma_verbose(mol, no_bounds=True, bounds={"log"})
# ValueError("bounds specified but `no_bounds` is True.")
```

Note that due to the multiprocessing features of `orca` the `space` outputs may be slightly different. More details can be found in [Seet et al 2024](https://arxiv.org/abs/2410.09100), and *TODO* JOSS Link.

# Cross-platform support
Rust, `maturin`, and `cargo` facilitate robust cross platform support for ORCA and `pyorca` requires `RDKit` as a dependency. Accordingly, `pyorca` is only available on those platforms that have `RDKit` support through PyPI. These include `windows-x64`, `macos-x86`, `macos-aarch64`, `ubuntu-x86`, and `ubuntu-aarch64`. If you're using a different platform and already have `RDKit` installed (for example through `conda`), then `pyorca` may work but we offer no promises.   