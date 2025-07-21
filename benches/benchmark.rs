use std::{
    collections::HashMap,
    ffi::OsStr,
    fs,
    path::Path,
    time::{Duration, Instant},
};

use bit_set::BitSet;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use assembly_theory::{
    assembly::{index_search, ParallelMode},
    bounds::Bound,
    canonize::{canonize, CanonizeMode},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    kernels::KernelMode,
    loader::parse_molfile_str,
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
        (EnumerateMode::GrowErode, "growerode"),
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
    let canonize_modes = [(CanonizeMode::Nauty, "nauty")];

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
                        for _ in 0..iters {
                            for mol in &mol_list {
                                let subgraphs = enumerate_subgraphs(mol, EnumerateMode::GrowErode);
                                let start = Instant::now();
                                let mut isomorphism_classes = HashMap::<_, Vec<BitSet>>::new();
                                for subgraph in &subgraphs {
                                    isomorphism_classes
                                        .entry(canonize(mol, &subgraph, *canonize_mode))
                                        .and_modify(|bucket| bucket.push(subgraph.clone()))
                                        .or_insert(vec![subgraph.clone()]);
                                }
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
/// of [`Bound`]s. This benchmark disables other search optimizations (e.g.,
/// kernelization or memoization) to focus on the effects of bounds only.
pub fn bench_bounds(c: &mut Criterion) {}

pub fn bench_index_search(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bench_index_search");

    // Define datasets and algorithm modes.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let bound_lists = [
        (vec![], "naive"),
        (vec![Bound::Log], "logbound"),
        (vec![Bound::Int], "intbound"),
        (
            vec![Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
            "allbounds",
        ),
    ];

    // Run the benchmark for each dataset and algorithm mode.
    for dataset in &datasets {
        let mol_list = load_dataset_molecules(dataset);
        for (bounds, name) in &bound_lists {
            bench_group.bench_with_input(BenchmarkId::new(*dataset, &name), bounds, |b, bounds| {
                b.iter(|| {
                    for mol in &mol_list {
                        index_search(
                            &mol,
                            EnumerateMode::GrowErode,
                            CanonizeMode::Nauty,
                            ParallelMode::Always,
                            KernelMode::None,
                            &bounds,
                            false,
                        );
                    }
                });
            });
        }
    }

    bench_group.finish();
}

criterion_group! {
    name = benchmark;
    config = Criterion::default().sample_size(20);
    targets = bench_enumerate, bench_canonize, bench_index_search
}
criterion_main!(benchmark);
