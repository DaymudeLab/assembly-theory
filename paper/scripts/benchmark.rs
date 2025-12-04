use std::{ffi::OsStr, fs, path::Path};

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use assembly_theory::{
    assembly::{index_search, ParallelMode},
    bounds::Bound,
    canonize::CanonizeMode,
    kernels::KernelMode,
    loader::parse_molfile_str,
    memoize::MemoizeMode,
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

/// Benchmark the entire [`index_search`] function for different [`Bound`]s.
pub fn joss_bench(c: &mut Criterion) {
    // Define a new criterion benchmark group for the JOSS manuscript.
    let mut bench_group = c.benchmark_group("joss");

    // Define datasets and bound lists.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let bound_lists = [
        (vec![], "no-bounds"),
        (vec![Bound::Log], "log"),
        (vec![Bound::Int], "int"),
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
                    b.iter(|| {
                        for mol in &mol_list {
                            index_search(
                                &mol,
                                CanonizeMode::TreeNauty,
                                ParallelMode::DepthOne,
                                MemoizeMode::CanonIndex,
                                KernelMode::None,
                                &bounds,
                            );
                        }
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
    targets = joss_bench
}
criterion_main!(benchmark);
