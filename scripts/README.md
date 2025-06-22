# `assembly-theory` Scripts

As with many open-source projects, `scripts/` is the junk drawer of supporting code external to the `assembly-theory` library itself.


## Installation

We use [`uv`](https://docs.astral.sh/uv/) to manage Python environments. [Install it](https://docs.astral.sh/uv/getting-started/installation/) and then run the following to get all dependencies:

```shell
# Make sure you're in the scripts/ directory!
uv sync
```


## Scripts Catalog

### `dataset_curation.ipynb`

This Python notebook generates and documents the molecule datasets in `data/` curated from existing databases.
To view and interact with the notebook, run:

```shell
uv run jupyter notebook dataset_curation.ipynb
```

Or, if you just want to run the code:

```shell
uv run jupyter execute dataset_curation.ipynb
```

> [!NOTE]
> The notebook may take a while to run, since it downloads full datasets from Zenodo and then extracts the desired subsets. 
> Jupyter does not provide a way of redirecting the internal kernel's stdout to the terminal, so the various progress update messages are invisible when run in this way.


### `generate-ma-index.sh`

While the above Python notebook curates the reference dataset `.mol` files, this script generates their ground truth assembly indices in `data/<dataset>/ma-index.csv`.
This file is needed before a dataset can be used for testing.

From the `scripts/` directory, run this script to be presented with interactive menus for generation:

```shell
./generate-ma-index.sh
```

The first menu asks you to choose the dataset to generate ground truth for.
Enter a number to choose the dataset:

```shell
1) ../data/checks       3) ../data/gdb13_1201
2) ../data/coconut_55   4) ../data/gdb17_200
Generate ma-index.csv for: 
```

The next menu asks you which program to use for assembly index calculation.
Again, enter a number to choose:

```shell
1) assembly_go (Jirasek et al., 2024)
2) assembly_cpp (Seet et al., 2024)
3) assembly-theory
Calculate assembly indices using: 
```

[`assembly_go`](https://github.com/croningp/assembly_go) is existing, open-source software for calculating assembly indices.
We do not package its source code or its executable with our library, but it can be obtained [on GitHub](https://github.com/croningp/assembly_go) if non-self-referential ground truth is desired.
[`assembly_cpp`] is the current state-of-the-art algorithm by Seet et al. (2024) and was provided to us by its authors on the condition that it remains private and is used only for this ground-truth generation.
Otherwise, a release build of `assembly-theory` is created and used.
