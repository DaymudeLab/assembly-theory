# Open, Reproducible Calculation of Assembly Indices

This repository contains the Rust source code for a faster
calculation of molecular assembly indices. This is a collaboration in the
Biodesign Center for Biocomputing, Security and Society at Arizona State
University involving (in alphabetical order) Cole Mathis, Devansh Vimal,
Devendra Parkar, Joshua Daymude, Garrett Parzych, Olivia Smith, and Sean
Bergen.

## Functionality and Usage Examples

`assembly-theory` can be used to compute assembly indices as a standalone executable, as a library imported by other Rust code, or via a Python interface.
Here, we provide usage examples of each; in the next section, we demonstrate testing and benchmarking functionality.


### Building and Running the Executable

Rust provides the `cargo` build system and package manager for dependency management, compilation, packaging, and versioning.
To build the standalone executable, run:

```shell
cargo build --release
```

This creates an optimized, portable, standalone executable named `target/release/assembly-theory`.
It takes as input a path to a `.mol` file and returns that molecule's integer assembly index:

```shell
> ./target/release/assembly-theory data/checks/anthracene.mol
6
```

Running with the `--verbose` flag provides additional information, including the input molecule's *number of disjoint, isomorphic subgraph pairs* (i.e., the number of times any molecular substructure is repeated inside the molecule) and the size of the top-down algorithm's *search space* (i.e., its total number of recursive calls).

```shell
> ./target/release/assembly-theory data/checks/anthracene.mol --verbose
Assembly Index: 6
Duplicate subgraph pairs: 406
Search Space: 3143
```

By default, `assembly-theory` uses its fastest algorithm for assembly index calculation (currently `assembly-theory`-allbounds, see the previous section).
To use a specific bound or disable bounds altogether, set the `--bounds` or `--no-bounds` flags:

```shell
# naive, no bounds
./target/release/assembly-theory <molpath> --no-bounds

# logbound, only logarithmic bound (Jirasek et al., 2024)
./target/release/assembly-theory <molpath> --bounds log

# intbound, only integer addition chain bound (Seet et al., 2024)
./target/release/assembly-theory <molpath> --bounds int-chain

# allbounds, both integer and vector addition chain bounds
./target/release/assembly-theory <molpath> --bounds int-chain vec-chain
```

Finally, the `--molecule-info` flag prints the molecule's graph representation as a vertex and edge list, the `--help` flag prints a guide to this command line interface, and the `--version` flag prints the current `assembly-theory` version.


### Installing and using the Python library

The python library uses `maturin` as a build tool. This needs to be run in a virtual environment. Use the following commands to build and install the library:
```shell
pip install maturin
maturin develop
```

This library computes the assembly index of molecules using RDKit's `Mol` class. Here's a basic example:

```python
import assembly_theory as at
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
at.molecular_assembly(anthracene)  # 6
```

### Core Functions  

The python library provides three main functions:

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

## Contributing
See [`HACKING`](HACKING.md)

## Known Issues

The current implementation follows the approach of Seet et. al. 2024 and
enumerates all duplicable subgraphs of the input molecule. The size of the
enumeration is stored in a usize. If there are enough duplicate subgraph pairs
in the molecule, then it is possible for the usize to overflow, resulting in a
panic. This behavior was observed in previous versions which used a u32 on
molecules in the coconut_220 benchmark, only using the naive implementation.
This problem is unavoidable as chemical space is vast, and naive enumeration of
it is a bad idea. In principle a molecular graph could be constructed such that
its duplicable subgraph enumeration would overflow arbitrary memory. Such an
error is unlikely to occur on reasonable compute time-scales. The discovery and
discussion of this issue is documented in
[Issue #49](https://github.com/DaymudeLab/assembly-theory/issues/49).

## License

Licensed under either of the [Apache License, Version 2.0](https://choosealicense.com/licenses/apache-2.0/) or the [MIT License](https://choosealicense.com/licenses/mit/), at your option.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in this repository/package (as defined by the Apache-2.0 License) shall be dual-licensed as above, without any additional terms or conditions.
