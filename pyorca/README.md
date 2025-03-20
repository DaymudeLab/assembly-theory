# `assembly-theory` Python Library

This is a Python library for computing the assembly index of molecules with high-performance Rust-based calculations. It integrates with `RDKit` and is built using `maturin`.

## Installation  

First, create and activate a virtual environment:  

### Windows:
```shell
python -m venv at_env
at_env\Scripts\activate
```

### macOS & Unix:
```shell
python -m venv at_env
source at_env/bin/activate
```

Next, install `maturin` and build the library:
```shell
pip install maturin
maturin develop
```

## Running Tests

To run the test suite, install `pytest` and execute the tests from the top-level `assembly-theory` directory:

```shell
pip install pytest
pytest test
```

## Example Usage

`assembly-theory` computes the assembly index of molecules using RDKit's `Mol` class. Here's a basic example:

```python
import assembly_theory as at
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
at.molecular_assembly(anthracene)  # 6
```

## Core Functions  

`assembly-theory` provides three main functions:

- **`molecular_assembly(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> int`**  
  Computes the assembly index of a given molecule.
  - `timeout` (in seconds) sets a limit on computation time, raising a `TimeoutError` if exceeded.  
  - `serial=True` forces a serial execution mode, mainly useful for debugging.


- **`molecular_assembly_verbose(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> dict`**  
  Returns additional details, including the number of duplicated isomorphic subgraphs (`duplicates`) and the size of the search space (`space`).  
  - `timeout` (in seconds) sets a limit on computation time, raising a `TimeoutError` if exceeded.  
  - `serial=True` forces a serial execution mode, mainly useful for debugging.

- **`molecule_info(mol: Chem.Mol) -> str`**  
  Returns a string representation of the molecule’s atom and bond structure for debugging.

## Search Strategy Options  

Both `molecular_assembly` and `molecular_assembly_verbose` support optional parameters for controlling the search strategy in the branch-and-bound algorithm:

- **`bounds: set[str]`** – Specifies heuristic bounds used to optimize the search.  
  - Options: `{"log"}`, `{"intchain"}`, or `{"log", "intchain"}`.  
  - Defaults to the best-performing option when not specified.

- **`no_bounds: bool`** – If `True`, disables all bounds, forcing an exhaustive search of all pathways.

The effect of these options can be observed using `molecular_assembly_verbose`:

```python
from rdkit import Chem
import assembly_theory as at

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")

at.molecular_assembly_verbose(anthracene, bounds={"log"})
# {'index': 6, 'duplicates': 418, 'space': 40507}

at.molecular_assembly_verbose(anthracene, bounds={"intchain"})
# {'index': 6, 'duplicates': 418, 'space': 3484}

at.molecular_assembly_verbose(anthracene, bounds={"intchain", "log"})
# {'index': 6, 'duplicates': 418, 'space': 3081}

at.molecular_assembly_verbose(anthracene, no_bounds=True)
# {'index': 6, 'duplicates': 418, 'space': 129409}
```

Due to multiprocessing, `space` outputs may vary slightly. More details can be found in [Seet et al. 2024](https://arxiv.org/abs/2410.09100) (*TODO: JOSS Link*).

## Cross-Platform Support  

`assembly-theory` leverages Rust, `maturin`, and `cargo` for robust cross-platform support. However, since `assembly-theory` depends on `RDKit`, it is only available on platforms where `RDKit` is supported via PyPI, including:  

- `windows-x64`  
- `macos-x86`  
- `macos-aarch64`  
- `ubuntu-x86`  
- `ubuntu-aarch64`  

If you are using a different platform and have `RDKit` installed (e.g., via `conda`), `assembly-theory` **may** work, but we do not guarantee compatibility.
