use std::{
    collections::HashMap,
    ffi::OsStr,
    fs,
    path::Path,
    sync::{atomic::AtomicUsize, Arc},
    time::{Duration, Instant},
};

use bit_set::BitSet;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use assembly_theory::{
    assembly::{matches, recurse_index_search, ParallelMode},
    bounds::Bound,
    canonize::{canonize, CanonizeMode},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    loader::parse_molfile_str,
    memoize::{Cache, MemoizeMode},
    molecule::Molecule,
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

/// Benchmark the first step of [`index_search`] which enumerates all connected
/// non-induced subgraphs with at most |E|/2 edges.
///
/// This benchmark preloads all dataset .mol files as [`Molecule`]s and then
/// times only the [`enumerate_subgraphs`] function for each [`EnumerateMode`].
pub fn bench_enumerate(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_enumerate");

    // Define datasets and enumeration modes. EnumerateMode::ExtendIsomorphic
    // is not included here because it combines enumeration and canonization.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let enumerate_modes = [
        (EnumerateMode::Extend, "extend"),
        (EnumerateMode::GrowErode, "grow-erode"),
    ];

    // Run a benchmark for each dataset and enumeration mode.
    for dataset in &datasets {
        let mol_list = load_dataset_molecules(dataset);
        for (enumerate_mode, name) in &enumerate_modes {
            bench_group.bench_with_input(
                BenchmarkId::new(*dataset, &name),
                &enumerate_mode,
                |b, &enumerate_mode| {
                    b.iter(|| {
                        for mol in &mol_list {
                            enumerate_subgraphs(mol, *enumerate_mode);
                        }
                    });
                },
            );
        }
    }

    bench_group.finish();
}

/// Benchmark the second step of [`index_search`] which bins connected,
/// non-induced subgraphs into isomorphism classes.
///
/// This benchmark preloads all dataset .mol files as [`Molecule`]s and uses
/// the fastest option for [`enumerate_subgraphs`] to get their subgraphs. It
/// times only the creation of isomorphism classes for each [`CanonizeMode`].
pub fn bench_canonize(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_canonize");

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
                    b.iter_custom(|iters| {
                        let mut total_time = Duration::new(0, 0);
                        for mol in &mol_list {
                            // Precompute subgraph enumeration.
                            let subgraphs = enumerate_subgraphs(mol, EnumerateMode::GrowErode);

                            // Benchmark the isomorphism class creation.
                            for _ in 0..iters {
                                let start = Instant::now();
                                let mut isomorphism_classes = HashMap::<_, Vec<BitSet>>::new();
                                subgraphs.iter().for_each(|subgraph| {
                                    isomorphism_classes
                                        .entry(canonize(mol, subgraph, *canonize_mode))
                                        .and_modify(|bucket| bucket.push(subgraph.clone()))
                                        .or_insert(vec![subgraph.clone()]);
                                });
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

/// Benchmark the search step of [`index_search`] using different [`Bound`]s.
///
/// This benchmark precomputes the enumeration and isomorphism steps using the
/// fastest options and times only the search step for different combinations
/// of [`Bound`]s. This benchmark otherwise uses the default search options.
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
                            let matches = matches(
                                mol,
                                EnumerateMode::GrowErode,
                                CanonizeMode::TreeNauty,
                                ParallelMode::DepthOne,
                            );
                            let mut init = BitSet::new();
                            init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
                            let edge_count = mol.graph().edge_count();

                            // Benchmark the search phase.
                            for _ in 0..iters {
                                let mut cache =
                                    Cache::new(MemoizeMode::CanonIndex, CanonizeMode::TreeNauty);
                                let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
                                let start = Instant::now();
                                recurse_index_search(
                                    mol,
                                    &matches,
                                    Vec::new(),
                                    &[init.clone()],
                                    edge_count - 1,
                                    best_index,
                                    edge_count,
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
/// This benchmark precomputes the enumeration and isomorphism steps using the
/// fastest options and times only the search step for different
/// [`MemoizeMode`]s. This benchmark otherwise uses the default search options.
pub fn bench_memoize(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_memoize");

    // Define datasets and bound lists.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let memoize_modes = [
        (MemoizeMode::None, CanonizeMode::Nauty, "no-memoize"),
        (MemoizeMode::FragsIndex, CanonizeMode::Nauty, "frags-index"),
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
                            let matches = matches(
                                mol,
                                EnumerateMode::GrowErode,
                                CanonizeMode::TreeNauty,
                                ParallelMode::DepthOne,
                            );
                            let mut init = BitSet::new();
                            init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
                            let edge_count = mol.graph().edge_count();

                            // Benchmark the search phase.
                            for _ in 0..iters {
                                let mut cache = Cache::new(memoize_mode, canonize_mode);
                                let best_index = Arc::new(AtomicUsize::from(edge_count - 1));
                                let start = Instant::now();
                                recurse_index_search(
                                    mol,
                                    &matches,
                                    Vec::new(),
                                    &[init.clone()],
                                    edge_count - 1,
                                    best_index,
                                    edge_count,
                                    &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
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
    targets = bench_enumerate, bench_canonize, bench_bounds, bench_memoize
}
criterion_main!(benchmark);
