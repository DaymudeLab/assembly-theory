# `assembly-theory` JOSS Manuscript

This `paper/` directory contains the source files for our [Journal of Open Source Science](https://joss.theoj.org/) manuscript on the `assembly-theory` software package titled *Open, Reproducible Calculation of Assembly Indices*.
This is a collaboration in the Biodesign Center for Biocomputing, Security and Society at Arizona State University involving Devansh Vimal, Garrett Parzych, Olivia M. Smith, Devendra Parkar, Holly Bergen, Joshua J. Daymude, and Cole Mathis.

> [!NOTE]
> You should only be able to see this directory on the `joss-paper` branch.
> This branch is never meant to be merged into `main`.


## What's Here?

The source files for our JOSS manuscript are `paper.md`, `paper.bib`, and `figures/`.
Each time commits are pushed to this `joss-paper` branch, a draft PDF is compiled using [GitHub Actions](https://github.com/DaymudeLab/assembly-theory/actions/workflows/paper.yml).
Instructions for compiling a paper PDF locally can be found in the [JOSS docs](https://joss.readthedocs.io/en/latest/paper.html#locally).

If you additionally wish to reproduce our benchmarks and associated figures (instructions below), you will need Rust, Go, Python, and `cmake`.
After cloning this repository, run the following to also pull the submodules:

```shell
git submodule init
git submodule update
```

We use [`uv`](https://docs.astral.sh/uv/) to manage Python environments.
[Install it](https://docs.astral.sh/uv/getting-started/installation/) and then run the following to get all dependencies:

```shell
uv sync
```


## Benchmarking Method and Instructions for Reproduction

This paper includes benchmarks of our `assembly-theory` assembly index calculations against those of [`assembly_go`](https://github.com/croningp/assembly_go), an earlier implementation written in Go ([Jirasek et al., 2024](https://doi.org/10.1021/acscentsci.4c00120)), and [`assemblycpp-v5`](https://github.com/croningp/assemblycpp-v5), a state-of-the-art implementation written in C++ ([Seet et al., 2024](https://arxiv.org/abs/2410.09100)).
For the purposes of reproducibilty, this repository includes the versions of `assembly_go` and `assemblycpp-v5` that we benchmarked as submodules.
The benchmarks measure only the runtime required to compute assembly pathways and indices, but exclude all setup and teardown (e.g., loading `.mol` files into internal molecule/graph representations).
The molecule datasets used for benchmarking are described in `paper.md` and can be found in the top-level `data/` directory.

> [!WARNING]
> Some of these benchmarks take a long time to run, especially when averaging over many samples on large reference datasets.


### Benchmarking `assembly-theory`

Set up the `assembly-theory` benchmark by copying the paper-specific Rust benchmark file into the top-level `benches/` directory:

```shell
cp paper/scripts/benchmark.rs benches/
```

Then run the benchmark with

```shell
cargo bench
```


### Benchmarking `assembly_go`

Set up the `assembly_go` benchmark by copying the Go benchmark file into the appropriate submodule and then going to the corresponding directory:

```shell
cp paper/scripts/main_test.go paper/assembly_go/cmd/app/
cd paper/assembly_go/cmd/app
```

Then run the benchmark with

```shell
go test -bench=. -cpu=<cpus> -count=<iters> -timeout=0 > datasets_bench.tsv
```

where `<cpus>` is replaced by the number of CPUs you want to let `assembly_go` parallelize over and `<iters>` is replaced by the number of iterations you want to run the benchmark and average the times over (see the [go testing flags](https://pkg.go.dev/cmd/go#hdr-Testing_flags) for details).
For our paper, we used `-cpu=16` and `-count=20`.

> [!NOTE]
> The benchmark for `assembly_go` on `coconut_55` is very slow, so we only ran that version of the benchmark once (i.e., `-count=1`).


### Benchmarking `assemblycpp-v5`

TODO


### Getting Benchmark Results

From this `paper/` directory, run the following to get the benchmark statistics.

```
uv run scripts/bench_stats.py
```

This script reports the mean benchmark time and 95% confidence interval of the mean for each algorithm&ndash;dataset pair.
