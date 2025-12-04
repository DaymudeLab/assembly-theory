---
title: 'Open, Reproducible Calculation of Assembly Indices'
tags:
  - assembly theory
  - biochemistry
  - astrobiology
  - Rust
authors:
  - name: Devansh Vimal
    orcid: 0009-0006-2794-8995
    affiliation: 1
  - name: Garrett Parzych
    orcid: 0009-0008-4789-9603
    affiliation: "1, 2"
  - name: Olivia M. Smith
    orcid: 0009-0004-2299-3522
    affiliation: "1, 3"
  - name: Devendra Parkar
    orcid: 0009-0009-0133-8875
    affiliation: "1, 2"
  - name: Sean Bergen
    orcid: 0009-0004-3570-5120
    affiliation: "1, 2"
  - name: Joshua J. Daymude
    orcid: 0000-0001-7294-5626
    affiliation: "1, 2"
  - name: Cole Mathis
    orcid: 0000-0001-8424-9169
    corresponding: true
    affiliation: "1, 3"
affiliations:
  - name: Biodesign Center for Biocomputing, Security and Society, Arizona State University, United States
    index: 1
  - name: School of Computing and Augmented Intelligence, Arizona State University, United States
    index: 2
  - name: School of Complex Adaptive Systems, Arizona State University, United States
    index: 3
date: 6 August 2025
bibliography: paper.bib
---

# Summary

We present `assembly-theory`, an open-source, high-performance library for computing *assembly indices* of covalently bonded molecular structures.
This is a key complexity measure of *assembly theory*, a recent theoretical framework quantifying selection across diverse systems, most importantly chemistry.
`assembly-theory` is designed for researchers and practitioners alike, providing (i) extensible, high-performance Rust implementations of assembly index calculation algorithms, (ii) comprehensive tests and benchmarks against which current and future algorithmic improvements can be evaluated, and (iii) Python bindings to support integration with existing computational pipelines.



# Background

*Assembly theory* (AT) is a recently developed body of theoretical and empirical work characterizing selection in diverse physical systems [@Sharma2023-assemblytheory; @Walker2024-experimentallymeasured].
In AT, objects are entities that are finite, distinguishable, decomposable, and persistent in time.
AT characterizes objects by their *assembly index*, the minimum number of recursive subconstructions required to construct the object starting from a given set of building blocks [@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation].
To date, AT has most commonly been applied to molecular chemistry, where bonds are the basic building blocks and the quantity of interest is the *molecular assembly index* (MA); see \autoref{fig:assemblyindex} for an example.
MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy [@Jirasek2024-investigatingquantifying], enabling a novel approach to life detection [@Marshall2021-identifyingmolecules].
It has also been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into evolutionary history [@Liu2021-exploringmapping; @Kahana2024-constructingmolecular].



# Statement of Need

Computing MA efficiently remains a challenge.
In general, exact MA calculation is NP-hard [@Kempes2024-assemblytheory].
Previous software to compute MA have been approximate, closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
The original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines [@Marshall2021-identifyingmolecules].
Machine learning methods only provide approximate MA values [@Gebhard2022-inferringmolecular].
The more recent `assembly_go` implementation computes MA exactly but is written in Go, yielding worse performance and posing an accessibility barrier for most scientists who are unfamiliar with the language [@Jirasek2024-investigatingquantifying].
Finally, the latest `assemblycpp` C++ implementation achieves significant performance milestones through a branch-and-bound approach, but until very recently was not publicly available for comparison or verification by the community [@Seet2024-rapidcomputation].

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=100% }

With `assembly-theory`, we provide an open-source, fully documented, extensible, and high-performance library for assembly index calculation, implementing not just one algorithm, but a framework and source of ground truth within which current and future algorithmic approaches can be validated and compared.
The main implementation is written in Rust, which we chose for its cross-platform support, memory-safety, performant runtime, convenient parallelism, and integrated testing and documentation [@Perkel2020-whyscientists].
We also provide Python bindings, enabling ease of use for scientific practitioners.



# Functionality and Usage Examples

`assembly-theory` is available as a standalone executable, a Rust crate, and a Python package.
Here, we provide usage examples of each; in the next section, we describe testing and benchmarking functionality.


## Standalone Executable

Rust provides the `cargo` build system and dependency manager for compilation, testing, benchmarking, documentation, and packaging.
From source, build the executable with:

```shell
> cargo build --release
```

Simply pass this executable a `.mol` file to compute that molecule's assembly index:

```shell
> ./target/release/assembly-theory data/checks/anthracene.mol  # 6
```

A full list of options for customizing the assembly index calculation procedure is obtained with:

```shell
> ./target/release/assembly-theory --help
```


## Rust Crate

The publicly released `assembly-theory` Rust crate is hosted on [crates.io](https://crates.io/crates/assembly-theory) and can be included in a broader Rust project with:

```shell
> cargo add assembly-theory
```

Complete documentation of the crate is available on [docs.rs](https://docs.rs/assembly-theory); a simple usage example is:

```rust
use assembly_theory::{
    assembly::index,
    loader::parse_molfile_str
};

// Load a molecule from a `.mol` file.
let path = PathBuf::from(format!("./data/checks/anthracene.mol"));
let molfile = fs::read_to_string(path)?;
let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");

// Compute the molecule's assembly index.
assert_eq!(index(&anthracene), 6);
```


## Python Package

We use the [`pyo3`](https://crates.io/crates/pyo3) crate and [`maturin`](https://github.com/PyO3/maturin) to create a Python package that mirrors the Rust crate's public functionality.
This Python package can be installed from [PyPI](https://pypi.org/project/assembly-theory/) in the usual way:

```shell
> pip install assembly-theory
```

This package is designed for compatibility with `RDKit`, the standard Python library for cheminformatics [@2024-rdkitopensource].
Molecules can be loaded and manipulated using the `rdkit.Chem.Mol` class and then passed to our functions for assembly index calculation:

```python
import assembly_theory as at
from rdkit import Chem

# Get a mol block from a molecule's SMILES representation.
anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
anthracene = Chem.MolToMolBlock(anthracene)

# Calculate the molecule's assembly index.
at.index(anthracene)  # 6
```



# Tests and Benchmarks

`assembly-theory` includes test and benchmark suites for software validation and performance evaluation, respectively.
Both use curated reference datasets representing different classes of molecules, chosen for their structural diversity and approachable runtime on commodity hardware.
These reference data are sampled from:

- GDB-13, a database of enumerated chemical structures containing Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur, and Chlorine that are constrained only by valence rules and quantum mechanics [@Blum2009-970million].
- GDB-17, an extension of GDB-13 which includes additional nuclei such as the halogens Flourine and Iodine [@Ruddigkeit2012-enumeration166].
- KEGG COMPOUND, a database of small molecules, biopolymers, and other biologically relevant substances [@Kanehisa2000-keggkyoto; @Kanehisa2019-understandingorigin; @Kanehisa2023-keggtaxonomybased].
- COCONUT, a database of natural products (secondary metabolites) offering a rich source of evolved chemical complexity [@Sorokina2021-coconutonline; @Chandrasekhar2025-coconut20].

The `assembly-theory` test suite (run with `cargo test`) contains unit tests validating internal functionality and integration tests verifying the calculation of correct assembly indices for all molecules in our reference datasets.
Each reference dataset contains an `ma-index.csv` file with ground truth assembly indices calculated using `assemblycpp` [@Seet2024-rapidcomputation].

Our benchmark suite (run with `cargo bench`) evaluates the performance of each granular phase of assembly index calculation over entire reference datasets.
We leverage the [`criterion`](https://bheisler.github.io/criterion.rs/criterion/) Rust crate to automatically collect detailed timing statistics and create performance reports.



# Availability and Governance

`assembly-theory` is available as a source code repository on [GitHub](https://github.com/DaymudeLab/assembly-theory), as a Rust crate on [crates.io](https://crates.io/crates/assembly-theory), and as a Python package on [PyPI](https://pypi.org/project/assembly-theory/).
Following the standard practice for Rust projects, `assembly-theory` is dual-licensed under the MIT and Apache-2.0 licenses.
Contributing guidelines and project governance are described in our `README`.



# Author Contributions

DV was the primary software developer (architecture, command line interface, molecule representations, unit tests, parallelism, performance engineering).
GP, DV, and CM formalized the core algorithm design.
GP and SB implemented the algorithm's bounding strategies.
DP and DV implemented the `.mol` file parser.
CM and JJD implemented the Python interface.
OMS curated all reference datasets and assembly index ground truths with input from CM.
JJD created the integration tests and benchmarks.
JJD and CM wrote the paper.



# Acknowledgements

GP and JJD are supported in part by NSF award CCF-2312537.
DV, OMS, and CM acknowledge support from the ASU Biodesign Institute.



# References
