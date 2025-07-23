# assembly-theory

`assembly-theory` is an open-source, high-performance library for computing *assembly indices* of molecular structures (see, e.g., [Sharma et al., 2023](https://doi.org/10.1038/s41586-023-06600-9); [Walker et al., 2024](https://doi.org/10.1098/rsif.2024.0367)).
This Python package wraps the functionality of the [`assembly-theory` Rust crate](https://crates.io/crates/assembly-theory) for easier interoperability with existing cheminformatic libraries (e.g., [RDKit](https://pypi.org/project/rdkit-pypi/)) and computational pipelines.

If you're looking for the Rust crate or standalone executable versions of `assembly-theory`, or if you're trying to build this Python package from source, see our [GitHub repository](https://github.com/DaymudeLab/assembly-theory).


## Installation

Install `assembly-theory` using a Python virtual environment manager of your choosing:

```shell
pip install assembly-theory   # Using pip.
pipx install assembly-theory  # Using pipx.
uv add assembly-theory        # Using uv.
```


## Usage

`assembly-theory` exposes several functions, all of which expecta "mol block" (i.e., the contents of a `.mol` file as a string) as input.
For example:

```python
import assembly_theory as at  # Note the '_' instead of '-'!

# Load a mol block from file.
with open('anthracene.mol') as f:
    mol_block = f.read()

# Calculate the molecule's assembly index.
at.index(mol_block)  # 6
```

Combine `assembly-theory` with [RDKit](https://pypi.org/project/rdkit-pypi/) (installed separately) if you need to manipulate molecular representations or incorporate assembly index calculations in a broader cheminformatics pipeline.

```python
import assembly_theory as at
from rdkit import Chem

# Get a mol block from a molecule's SMILES representation.
anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
anthracene = Chem.MolToMolBlock(anthracene)

# Calculate the molecule's assembly index.
at.index(anthracene)  # 6
```

See the [`assembly_theory::python` documentation](https://docs.rs/assembly-theory/latest/assembly_theory/python) for a complete list of functions exposed to this package along with usage examples.


## Citation

If you use this package in your own scientific work, please consider citing us:

```bibtex
Coming soon!
```


## License

`assembly-theory` is licensed under the [Apache License, Version 2.0](https://choosealicense.com/licenses/apache-2.0/) or the [MIT License](https://choosealicense.com/licenses/mit/), at your option.
