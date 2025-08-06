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
date: 6 March 2025
bibliography: paper.bib
---

# Summary

We present `assembly-theory`, a Rust package for computing *assembly indices* of covalently bonded molecular structures.
This is a key complexity measure of *assembly theory*, a recent theoretical framework quantifying selection across diverse systems, most importantly chemistry [@Walker2024-experimentallymeasured; @Sharma2023-assemblytheory].
`assembly-theory` is designed for researchers and practitioners alike, providing (i) extensible, high-performance implementations of assembly index calculation algorithms, (ii) comprehensive benchmarks against which current and future algorithmic improvements can be tested, and (iii) Python bindings and `RDKit`-compatible data loaders to support integration with existing computational pipelines.



# Background

*Assembly theory* (AT) is a recently developed body of theoretical and empirical work focused on characterizing selection in diverse physical systems, most importantly in chemistry [@Sharma2023-assemblytheory; @Walker2024-experimentallymeasured].
Objects are defined in AT as entities that are finite, distinguishable, decomposable, and persistent in time.
AT characterizes objects based on their *assembly index*, the minimum number of recursive subconstructions required to construct the object starting from a given set of building blocks [@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation].
The most commonly studied application domain of AT to date is molecular chemistry, where bonds act as the basic building blocks and the quantity of interest is the *molecular assembly index* (MA); see \autoref{fig:assemblyindex} for an example.
It has previously been shown that MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy [@Jirasek2024-investigatingquantifying], enabling a novel approach to life detection based on AT [@Marshall2021-identifyingmolecules].
Beyond life detection, AT and MA have been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into evolutionary history by inferring relationships directly from metabolomic data [@Liu2021-exploringmapping; @Kahana2024-constructingmolecular].

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=100% }



# Statement of Need

Despite AT's promising applications, computing MA efficiently remains a challenge.
In general, exact MA calculation is an NP-hard problem [@Kempes2024-assemblytheory]; i.e., the necessary computing resources are likely to grow exponentially with a molecule's number of bonds.
Previous software to compute MA have been approximate, closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
For example, the original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines [@Marshall2021-identifyingmolecules].
Machine learning methods only provide approximate MA values [@Gebhard2022-inferringmolecular].
The more recent `assembly_go` implementation computes MA exactly but is written in Go, yielding worse performance and posing an accessibility barrier for most scientists who are unfamiliar with the language [@Jirasek2024-investigatingquantifying].
Finally, the latest implementation achieves significant performance milestones through a branch-and-bound approach [@Seet2024-rapidcomputation].
It is again written in C++ but is not publicly available, prohibiting its use and verification by the community.

With `assembly-theory`, we provide a high-performance Rust package for MA calculation while also providing Python bindings for key functionality, offering the best efficiency without sacrificing ease of use for scientific practitioners.
We chose Rust for its advantages of cross-platform support, memory-safety, performant runtime, convenient parallelism, and integrated testing and documentation [@Perkel2020-whyscientists].
By including test and benchmark suites, we also lay a foundation for fair, reproducible comparisons of future algorithmic improvements.



# Design and Current Algorithms

`assembly-theory` is not a single algorithmic implementation of assembly index calculations; rather, it is a framework and source of ground truth within which a diversity of algorithmic approaches can be validated and compared.
We purposely designed `assembly-theory` with a modular algorithm interface and data structures that can be easily extended to handle new algorithmic developments introduced as AT matures.

Following prior work [@Marshall2021-identifyingmolecules; @Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation], `assembly-theory` currently implements a top-down approach with two phases that execute in sequence: (1) a search space enumeration and (2) a parallel branch-and-bound search.
The enumeration phase finds all pairs of isomorphic, edge-disjoint subgraphs of the given molecule.
Isomorphic subgraphs are binned into equivalence classes by their `nauty` canonical representations [@McKay2014-practicalgraph] and pairs of subgraphs within the same class are yielded if they are edge-disjoint.
In the search phase, the given molecule is recursively fragmented by removing duplicate subgraphs enumerated in the first phase.
The MA of smaller fragments are used to determine the MA of their parents.
Optionally, a bounding strategy can be used to improve search efficiency.
We briefly summarize the implemented branch-and-bound ("`bb`") variants below, but emphasize that `assembly-theory` is not limited to this top-down, recursive approach.

- `bb-naive` fully enumerates all non-duplicate assembly pathways to find the shortest one.
- `bb-logbound` improves over the naive method by eliminating any assembly pathways whose current length plus $\log_2b$ exceeds the length of the shortest assembly pathway found so far, where $b$ is the number of remaining bonds [@Jirasek2024-investigatingquantifying].
- `bb-intbound` uses a stronger lower bound on the number of remaining assembly steps provided by an integer addition chain [@Seet2024-rapidcomputation].
- `bb-allbounds` simultaneously applies the previous integer addition chain bound and a novel bound provided by a vector addition chain.



# Functionality and Usage Examples

`assembly-theory` can be used to compute assembly indices as a standalone executable, as a library imported by other Rust code, or via a Python interface.
Here, we provide usage examples of each; in the next section, we demonstrate testing and benchmarking functionality.


## Standalone Executable

Rust provides the `cargo` build system and package manager for dependency management, compilation, packaging, and versioning and the [crates.io](https://crates.io/crates/assembly-theory) registry for package distribution.
To install the standalone executable, run:

```shell
> cargo install assembly-theory
```

This executable takes as input a path to a `.mol` file and returns that molecule's assembly index:

```shell
> assembly-theory data/checks/anthracene.mol
6
```

Running with the `--verbose` flag provides additional information, including the input molecule's *number of disjoint, isomorphic subgraph pairs* (i.e., the number of times any molecular substructure is repeated inside the molecule) and the size of the top-down algorithm's *search space* (i.e., its total number of recursive calls).

```shell
> assembly-theory data/checks/anthracene.mol --verbose
Assembly Index: 6
Duplicate subgraph pairs: 406
Search Space: 2306
```

By default, `assembly-theory` parallelizes its recursive search over as many threads as the OS allows.
To disable parallelism, use the `--serial` flag.

Also by default, `assembly-theory` uses its fastest algorithm for assembly index calculation (currently `bb-allbounds`, see the previous section).
To use a specific bound or disable bounds altogether, set the `--bounds` or `--no-bounds` flags:

```shell
# bb-naive, no bounds
> assembly-theory data/checks/anthracene.mol --no-bounds

# bb-logbound, only logarithmic bound (Jirasek et al., 2024)
> assembly-theory data/checks/anthracene.mol --bounds log

# bb-intbound, only integer addition chain bound (Seet et al., 2024)
> assembly-theory data/checks/anthracene.mol --bounds int-chain

# bb-allbounds, both integer and vector addition chain bounds
> assembly-theory data/checks/anthracene.mol --bounds int-chain vec-chain
```

Finally, the `--molecule-info` flag prints the molecule's graph representation as a vertex and edge list, the `--help` flag prints a guide to this command line interface, and the `--version` flag prints the current `assembly-theory` version.


## Rust Library

To include `assembly-theory` in a broader Rust project, run:

```shell
> cargo add assembly-theory
```

Complete documentation of the library is available on [docs.rs](https://docs.rs/crate/assembly-theory/0.2.0); a simple usage example is:

```rust
use assembly_theory::*;

// Read a .mol file as an assembly_theory::molecule::Molecule.
let molfile = fs::read_to_string(path)?;
let anthracene = loader::parse_molfile_str(&molfile).expect("Parsing failed");

// Calculate the molecule's assembly index.
assert_eq!(assembly::index(&anthracene), 6);
```


## Python Interface

We use [`maturin`](https://github.com/PyO3/maturin) to repackage the `assembly-theory` Rust binaries as the `assembly_theory` package for Python.
Instructions for this build process can be found in our `README`; otherwise, the Python package can be obtained from PyPI in the usual way:

```shell
> pip install assembly_theory
```

Our Python interface is built for compatibility with `RDKit`, the standard Python library for cheminformatics [@2024-rdkitopensource].
Molecules can be loaded and manipulated using the `rdkit.Chem.Mol` class and then passed to our functions for assembly index calculation:

```python
>>> import assembly_theory as at
>>> from rdkit import Chem
>>> anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
>>> at.molecular_assembly(anthracene)
6
```

In detail, `assembly_theory` exposes three main functions:

1. **`molecular_assembly`**`(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> int`
2. **`molecular_assembly_verbose`**`(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> dict` 
3. **`molecule_info`**`(mol: Chem.Mol) -> str`

These correspond to (1) running the Rust `assembly-theory` executable to obtain only an assembly index, (2) running with the `--verbose` flag to also obtain the number of disjoint isomorphic subgraph pairs (`duplicates`) and search space size (`space`), and (3) running with the `--molecule-info` flag to obtain molecule information, respectively.
The `timeout` parameter is specific to the Python interface: when set to a non-`None` integer value, a `TimeoutError` is raised if assembly index calculation exceeds `timeout` seconds.


# Tests and Benchmarks

`assembly-theory` includes test and benchmark suites for software validation and performance evaluation, respectively.
Both use curated reference datasets representing different classes of molecules, arranged roughly in order of increasing molecular size and complexity:

- `gdb13_1201`: 1,201 small, organic molecular structures sampled from GDB-13, a database of enumerated chemical structures containing Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur, and Chlorine that are constrained only by valence rules and quantum mechanics but may not be chemically stable or synthesizable [@Blum2009-970million].
Our sample includes all 201 molecules with 4&ndash;5 heavy atoms and 200 randomly sampled molecules for each number of heavy atoms from 6&ndash;10.
These molecules' MA range from 2&ndash;9.
- `gdb17_200`: 200 organic molecular structures sampled from the GDB-17 database, which includes additional nuclei beyond GDB-13 such as the halogens Flourine and Iodine [@Ruddigkeit2012-enumeration166].
Compared to GDB-13, these molecules are typically larger and represent more structural diversity.
Our sample includes 50 randomly sampled molecules for each number of heavy atoms from 14&ndash;17.
These molecules' MA range from 7&ndash;16.
- `checks`: 15 named molecules (e.g., anthracene, aspirin, caffeine, morphine) from KEGG COMPOUND [@Kanehisa2000-keggkyoto; @Kanehisa2019-understandingorigin; @Kanehisa2023-keggtaxonomybased] primarily used for rapid testing.
These molecules' number of heavy atoms range from 5&ndash;21 and have MA from 3&ndash;14.
- `coconut_55`: 55 natural products sampled from the COCONUT database [@Sorokina2021-coconutonline], accessed in late 2024, prior to COCONUT 2.0 [@Chandrasekhar2025-coconut20].
Natural products (or secondary metabolites) are a rich source of evolved chemical complexity, often exhibiting drug-like properties.
Subsets of this database were used to benchmark recent algorithmic progress in [@Seet2024-rapidcomputation].
Our sample includes five randomly sampled molecules for each number of heavy atoms from 15&ndash;25.
These molecules' MA range from 7&ndash;16.

We curated these reference datasets for their structural diversity and approachable runtime on commodity hardware.
Larger, more demanding datasets can be added as needed.

The `assembly-theory` test suite contains unit tests validating internal functionality and integration tests verifying the calculation of correct assembly indices for all molecules in our reference datasets.
Each reference dataset contains an `ma-index.csv` file with ground truth assembly indices calculated using the closed-source [@Seet2024-rapidcomputation] algorithm, privately provided to us by the authors for this use only.

Our benchmark suite evaluates `assembly-theory` performance by running repeated assembly index calculations over individual molecules or entire reference datasets.
We leverage the [`criterion`](https://bheisler.github.io/criterion.rs/criterion/) package for Rust to automatically collect detailed timing statistics, charts, and estimates of performance improvements and regressions.
As an example, \autoref{tab:benchtimes} shows `assembly-theory` performance across our four reference datasets against that of `assembly_go` [@Jirasek2024-investigatingquantifying].
Depending on the dataset and choice of `assembly-theory` algorithm, `assembly-theory` outperforms `assembly_go` by one or two orders of magnitude.
The 6.5&ndash;120.0x speedup of `bb-logbound` over `assembly_go` most clearly represents the efficiency of Rust over Go, since both use the same branch-and-bound approach with a logarithmic bound.
Algorithmic improvements such as the `bb-allbounds` combination of an integer addition chain bound [@Seet2024-rapidcomputation] and our novel vector addition chain bound yield more dramatic speedups for larger molecules, like those up to 410.0x for `coconut_55`.
This internal comparison showcases `assembly-theory` as a framework capable of comparing multiple algorithmic approaches on equal footing, free of differences in underlying datasets or language-specific efficiency issues.

: \label{tab:benchtimes} Mean benchmark execution times for `assembly_go` [@Jirasek2024-investigatingquantifying] vs. `assembly-theory` across reference datasets.
The benchmark times the MA calculation of all molecules in a given dataset in sequence, excluding the time required to parse and load `.mol` files into internal molecular graph representations.
`assembly_go` uses its default parameters while `assembly-theory` tests each of its algorithm variants independently.
Each benchmark was run on a Linux machine with a 5.7 GHz Ryzen 9 7950X CPU (16 cores) and 64 GB of memory.
Means are reported over 20 samples per software&ndash;dataset pair, except those marked with an $\ast$ which have prohibitively long runtimes and thus ran only once.

|              | `assembly_go` | `bb-naive`   | `bb-logbound` | `bb-intbound`   | `bb-allbounds`   |
| ---------- | ----------: | --------: | -----------: | -----------: | -----------: | 
| `gdb13_1201` |       0.968 s |      0.147 s |       0.149 s |         **0.140 s** |          0.142 s |
| `gdb17_200`  |      46.189 s |      4.586 s |       3.446 s |         2.988 s |          **2.946 s** |
| `checks`     |     212.000 s |     12.635 s |       1.767 s |         1.357 s |          **1.297 s** |
| `coconut_55` | 1.48 h$^\ast$ |    175.687 s |      57.123 s |        13.188 s |         **12.996 s** |

If finer-grained timing insights are needed, `assembly-theory` can also benchmark assembly index calculations for each individual molecule in a reference dataset.
For example, \autoref{fig:timescatter} shows the calculation time of each molecule in `gdb17_200` for the four branch-and-bound algorithms.
This is useful for teasing out which molecules are "hard" and characterizing where algorithmic improvements make the largest impact.

![*Per-Molecule Benchmark Times*. The mean assembly index calculation time across 20 samples for each molecule (dot) in `gdb17_200` as a function of the molecule's number of duplicate isomorphic subgraphs, a measure roughly correlated with the molecule's size and complexity. The same four `assembly-theory` branch-and-bound algorithms from \autoref{tab:benchtimes} are shown here.\label{fig:timescatter}](figures/jossplot.pdf){ width=75% }



# Availability and Governance

`assembly-theory` source code and documentation are openly available on [GitHub](https://github.com/DaymudeLab/assembly-theory).
Following the standard practice for Rust packages, `assembly-theory` is dual-licensed under the MIT and Apache-2.0 licenses.
External feedback and code contributions are handled through the usual Issues and Pull Request interfaces; guidelines for contributions are listed in `HACKING.md`.
The project's *maintainers* (initially Vimal, Daymude, and Mathis) will govern the project using the committee model: high-level decisions about the project's direction require maintainer consensus, major code changes require majority approval, hotfixes and patches require at least one approval, new maintainers may be added by unanimous decision of the existing maintainers, and existing maintainers may step down with advance notice.



# Author Contributions

GP, DV, and CM formalized the branch-and-bound algorithm design.
GP and SB formalized the integer and vector addition chain bounds.
DV was the primary software developer (architecture, command line interface, molecule representations, unit tests, parallelism, performance engineering).
GP implemented all bound calculations.
DP and DV implemented the `.mol` file parser and dataset-based benchmarks.
CM implemented the Python interface.
OMS curated all reference datasets and assembly index ground truths with input from CM.
SB and JJD wrote the `assembly_go` benchmarks.
JJD conducted and analyzed the benchmarks shown in \autoref{tab:benchtimes}.
DP, SB, GP, and JJD produced \autoref{fig:timescatter}.
JJD and CM wrote the paper.



# Acknowledgements

GP and JJD are supported in part by NSF award CCF-2312537.
DV, OMS, and CM acknowledge support from the ASU Biodesign Institute.



# References
