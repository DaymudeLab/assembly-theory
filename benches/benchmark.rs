use std::{
    ffi::OsStr,
    fs,
    path::Path,
    sync::{atomic::AtomicUsize, Arc},
    time::{Duration, Instant},
};

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use assembly_theory::{
    assembly::{recurse_index_search, ParallelMode},
    bounds::Bound,
    canonize::CanonizeMode,
    loader::parse_molfile_str,
    matches::Matches,
    memoize::{Cache, MemoizeMode},
    molecule::Molecule,
    state::State,
};

/// Parse all .mol files in `dataset` as [`Molecule`]s.
fn load_dataset_molecules(dataset: &str) -> Vec<Molecule> {
    let paths = fs::read_dir(Path::new("data").join(dataset)).unwrap();
    let mut mol_list: Vec<Molecule> = Vec::new();
    for path in paths {
        let name = path.unwrap().path();
        if name.extension().and_then(OsStr::to_str) == Some("mol") {
            mol_list.push(
                parse_molfile_str(
                    &fs::read_to_string(name.clone())
                        .expect(&format!("Could not read file {name:?}")),
                )
                .expect(&format!("Failed to parse {name:?}")),
            );
        }
    }
    mol_list
}

/// Benchmark the first step of [`index_search`] which computes all structural
/// information related to "matches" (i.e., pairs of edge-disjoint, isomorphic,
/// connected, non-induced subgraphs).
///
/// This benchmark preloads all dataset .mol files as [`Molecule`]s and then
/// times only the [`Matches::new`] function for each [`CanonizeMode`].
pub fn bench_matches(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_matches");

    // Define datasets and canonization modes.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let canonize_modes = [
        (CanonizeMode::Nauty, "nauty"),
        (CanonizeMode::TreeNauty, "tree-nauty"),
    ];

    // Run a benchmark for each dataset and canonization mode.
    for dataset in &datasets {
        let mol_list = load_dataset_molecules(dataset);
        for (canonize_mode, name) in &canonize_modes {
            bench_group.bench_with_input(
                BenchmarkId::new(*dataset, &name),
                &canonize_mode,
                |b, &canonize_mode| {
                    b.iter(|| {
                        for mol in &mol_list {
                            Matches::new(mol, *canonize_mode);
                        }
                    });
                },
            );
        }
    }

    bench_group.finish();
}

/// Benchmark the search step of [`index_search`] using different [`Bound`]s.
///
/// This benchmark precomputes matches information using the fastest options
/// and times only the search step for different combinations of [`Bound`]s.
/// This benchmark otherwise uses the default search options.
pub fn bench_bounds(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_bounds");

    // Define datasets and bound lists.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let bound_lists = [
        (vec![], "no-bounds"),
        (vec![Bound::Log], "log"),
        (vec![Bound::Int], "int"),
        (
            vec![Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
            "int-vec",
        ),
        (vec![Bound::Int, Bound::MatchableEdges], "int-matchable"),
    ];

    // Run the benchmark for each dataset and bound list.
    for dataset in &datasets {
        let mol_list = load_dataset_molecules(dataset);
        for (bounds, name) in &bound_lists {
            bench_group.bench_with_input(
                BenchmarkId::new(*dataset, &name),
                &bounds,
                |b, &bounds| {
                    b.iter_custom(|iters| {
                        let mut total_time = Duration::new(0, 0);
                        for mol in &mol_list {
                            // Precompute the molecule's matches and setup.
                            let matches = Matches::new(mol, CanonizeMode::TreeNauty);
                            let state = State::new(mol);
                            let edge_count = mol.graph().edge_count();

                            // Benchmark the search phase.
                            for _ in 0..iters {
                                let mut cache =
                                    Cache::new(MemoizeMode::CanonIndex, CanonizeMode::TreeNauty);
                                let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
                                let states_searched = Arc::new(AtomicUsize::from(0));
                                let start = Instant::now();
                                recurse_index_search(
                                    mol,
                                    &matches,
                                    &state,
                                    best_index,
                                    states_searched,
                                    bounds,
                                    &mut cache,
                                    ParallelMode::DepthOne,
                                );
                                total_time += start.elapsed();
                            }
                        }
                        total_time
                    });
                },
            );
        }
    }

    bench_group.finish();
}

/// Benchmark the search step of [`index_search`] using different
/// [`MemoizeMode`]s.
///
/// This benchmark precomputes matches information using the fastest options
/// and times only the search step for different [`MemoizeMode`]s. This
/// benchmark otherwise uses the default search options.
pub fn bench_memoize(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_memoize");

    // Define datasets and bound lists.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let memoize_modes = [
        (MemoizeMode::None, CanonizeMode::Nauty, "no-memoize"),
        (MemoizeMode::CanonIndex, CanonizeMode::Nauty, "nauty-index"),
        (
            MemoizeMode::CanonIndex,
            CanonizeMode::TreeNauty,
            "tree-nauty-index",
        ),
    ];

    // Run the benchmark for each dataset and bound list.
    for dataset in &datasets {
        let mol_list = load_dataset_molecules(dataset);
        for (memoize_mode, canonize_mode, name) in &memoize_modes {
            bench_group.bench_with_input(
                BenchmarkId::new(*dataset, &name),
                &(memoize_mode, canonize_mode),
                |b, (&memoize_mode, &canonize_mode)| {
                    b.iter_custom(|iters| {
                        let mut total_time = Duration::new(0, 0);
                        for mol in &mol_list {
                            // Precompute the molecule's matches and setup.
                            let matches = Matches::new(mol, CanonizeMode::TreeNauty);
                            let state = State::new(mol);
                            let edge_count = mol.graph().edge_count();

                            // Benchmark the search phase.
                            for _ in 0..iters {
                                let mut cache = Cache::new(memoize_mode, canonize_mode);
                                let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
                                let states_searched = Arc::new(AtomicUsize::from(0));
                                let start = Instant::now();
                                recurse_index_search(
                                    mol,
                                    &matches,
                                    &state,
                                    best_index,
                                    states_searched,
                                    &[Bound::Int, Bound::MatchableEdges],
                                    &mut cache,
                                    ParallelMode::DepthOne,
                                );
                                total_time += start.elapsed();
                            }
                        }
                        total_time
                    });
                },
            );
        }
    }

    bench_group.finish();
}

criterion_group! {
    name = benchmark;
    config = Criterion::default().sample_size(20);
    targets = bench_matches, bench_bounds, bench_memoize
}
criterion_main!(benchmark);
