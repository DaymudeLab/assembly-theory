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
