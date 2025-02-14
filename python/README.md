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
import orca
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
orca.compute_ma(anthracene) # 7
```


# Details

ORCA provides three primary functions:

* `compute_ma(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False) -> int`
Computes the assembly index of a given molecule.
* `compute_ma_verbose(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False) -> dict`
Returns additional details, including the number of duplicated isomorphic subgraphs (duplicates) and the size of the search space (space).
* `get_molecule_info(mol: Chem.Mol) -> str`
Provides a string representation of the molecule’s atom and bond structure, useful for debugging.

## Search Strategy Options
Both compute_ma and compute_ma_verbose support optional parameters for controlling the search strategy in the branch-and-bound algorithm:

* `bounds: set[str]` – Specifies the heuristic bounds used to optimize the search.
* * Options: `{"log"}`, `{"addition"}`, or `{"log", "addition"}`.
* * Defaults to the best-performing option when not provided.
* `no_bounds: bool` – If `True`, disables all bounds, forcing a complete search of all pathways.

The effect of these on the search space can be see using `compute_ma_verbose` for example:

```python
from rdkit import Chem
import orca
anthra_smi = "c1ccc2cc3ccccc3cc2c1"
mol = Chem.MolFromSmiles(anthra_smi)
orca.compute_ma_verbose(mol, bounds={"log"})
# {'index': 6, 'duplicates': 418, 'space': 39015}
orca.compute_ma_verbose(mol, bounds={"addition"})
# {'index': 6, 'duplicates': 418, 'space': 2784}
orca.compute_ma_verbose(mol, bounds={"addition", "log"})
# {'index': 6, 'duplicates': 418, 'space': 2784}
orca.compute_ma_verbose(mol, no_bounds=True)
# {'index': 6, 'duplicates': 418, 'space': 129409}

# Invalid combination
orca.compute_ma_verbose(mol, no_bounds=True, bounds={"log"})
# ValueError("bounds specified but `no_bounds` is True.")
```

More details can be found in Ref JOSS and Ref Seet. TODO (add links)