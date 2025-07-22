# Open, Reproducible Calculation of Assembly Indices

`assembly-theory` is an open-source, high-performance library for computing *assembly indices* of molecular structures (see, e.g., [Sharma et al., 2023](https://doi.org/10.1038/s41586-023-06600-9); [Walker et al., 2024](https://doi.org/10.1098/rsif.2024.0367)).
It is implemented in Rust and is available as a [Rust crate](https://crates.io/crates/assembly-theory), [Python package](https://pypi.org/project/assembly-theory/), and standalone executable. 


## Getting Started

If you want to use the Rust crate in a Rust project, refer to the [docs.rs](https://docs.rs/assembly-theory) documentation for installation and usage examples.
If you want to use the Python library (e.g., to take advantage of RDKit-compatible molecule loaders), refer to the documentation on [PyPI](https://pypi.org/project/assembly-theory/).

Otherwise, clone/download this repository if you want to:

- Build and run the standalone executable
- Build and run tests and benchmarks on our reference datasets
- Build the Python package locally


### Requirements

Currently, this project only supports Unix-like systems (macOS and Linux).
Windows is supported through [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

You need Rust, installed either [using `rustup`](https://www.rust-lang.org/tools/install) or via your system package manager of choice.
This provides the `cargo` build system and dependency manager for compilation, testing, benchmarking, documentation, and packaging.

You will also need the [clang](https://clang.llvm.org) toolkit for C/C++.


### Using the Standalone Executable

Build an optimized (release) version of the standalone executable with:

```shell
cargo build --release
```

Simply pass this executable a `.mol` file path to compute that molecule's assembly index:

```shell
> ./target/release/assembly-theory data/checks/anthracene.mol
6
```

A full list of options for returning more information or customizing the assembly index calculation procedure can be found by running:

```shell
./target/release/assembly-theory --help
```


### Tests and Benchmarks

`assembly-theory` comes with a variety of unit, integration, and documentation example tests ensuring the correct calculation of assembly indices.
To run all tests, use:

```shell
cargo test
```

We actively encourage the use of `assembly-theory` as a framework within which new algorithmic improvements can be implemented and tested.
To measure the performance of a potential improvement, we've implemented benchmarks using the [`criterion`](https://crates.io/crates/criterion) crate.
These benchmarks run assembly index calculation against reference datasets of molecules, timing only the calculation part (skipping molecule parsing, etc.).
To run all benchmarks, use:

```shell
cargo bench
```

See the [`criterion` command line options](https://bheisler.github.io/criterion.rs/book/user_guide/command_line_options.html) for details on how to run only specific benchmarks or save baselines for comparison.


### Building the Python Package Locally

We use the [`pyo3`](https://crates.io/crates/pyo3) crate to package functionality from our Rust library as a Python package called `assembly_theory`.
To build this package locally, first install [`maturin`](https://pypi.org/project/maturin/) using your virtual environment manager of choice:

```shell
pip install maturin      # using pip
pipx install maturin     # using pipx
uv tool install maturin  # using uv
```

Then, within the corresponding virtual environment, build and install this Rust library as a Python package:

```shell
maturin develop
```

Once installed, this Python package can be combined with standard cheminformatic packages like [`RDKit`](https://www.rdkit.org/docs/index.html#)) to flexibly manipulate molecular representations and compute their assembly indices.

```python
import assembly_theory as at
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
at.index(anthracene)  # 6
```

See the documentation on [PyPI](https://pypi.org/project/assembly-theory/) for a complete list of functions exposed to the Python package.


## Known Issues

- The current implementation tallies the number of states searched during recursive assembly index calculation.
If the number of states searched exceeds the (very large) limit of a `usize`, the code panics.
This is unlikely to occur when various kernelization, memoization, and bounding strategies are enabled to prune the search space, but is theoretically always possible given a large enough molecule and sufficiently long search time.
See [#49](https://github.com/DaymudeLab/assembly-theory/issues/49) for details.


## Contributing

This project is under active development!
Any code on the `main` branch is considered usable, but not necessarily stable or feature-complete.
See our [releases](https://github.com/DaymudeLab/assembly-theory/releases) for more reliable snapshots of the project.

Have a suggestion for new features or a bug you need fixed?
Open a [new issue](https://github.com/DaymudeLab/assembly-theory/issues/new).

Want to contribute your own code?

- Familiarize yourself with the [Rust API Guidelines](https://github.com/DaymudeLab/assembly-theory/compare) and overall architecture of `assembly-theory`.
- Development team members should work in individual feature branches.
External contributors work in repository forks.
- Before opening a pull request onto `main`, make sure you have (1) rebased onto `main` and (2) run `cargo fmt` and `cargo clippy`.
- Open a [new pull request](https://github.com/DaymudeLab/assembly-theory/compare), provide a descriptive list of your changes (along with referencing any issues your PR resolves) and assign one of @AgentElement, @jdaymude, or @colemathis as a reviewer. 
Your PR will not be reviewed until it passes all GitHub Actions (compilation, formatting, tests, etc.).


## Governance

`assembly-theory` is maintained by Devansh Vimal (@AgentElement), Josh Daymude (@jdaymude), and Cole Mathis (@colemathis), with support from other members of the [Biodesign Center for Biocomputing, Security and Society](https://biodesign.asu.edu/biocomputing-security-and-society/) at Arizona State University including Garrett Parzych (@Garrett-Pz), Olivia M. Smith (@omsmith161), Devendra Parkar (@devrz45), and Sean Bergen(@ARandomCl0wn).

The maintainers govern the project using the committee model: high-level decisions about the project's direction require maintainer consensus, major code changes require majority approval, hotfixes and patches require just one maintainer approval, new maintainers can be added by unanimous decision of the existing maintainers, and existing maintainers can step down with advance notice.


## License

`assembly-theory` is licensed under either of the [Apache License, Version 2.0](https://choosealicense.com/licenses/apache-2.0/) or the [MIT License](https://choosealicense.com/licenses/mit/), at your option.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in this repository/package (as defined by the Apache-2.0 License) shall be dual-licensed as above, without any additional terms or conditions.
