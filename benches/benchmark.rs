use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::ffi::OsStr;
use std::fs;
use std::iter::zip;
use std::path::Path;

use assembly_theory::{
    assembly::index_search,
    bounds::Bound,
    loader,
    molecule::Molecule,
};

pub fn reference_datasets(c: &mut Criterion) {
    // Define a new criterion benchmark group of dataset benchmarks.
    let mut group = c.benchmark_group("reference_datasets");

    // Define datasets, bounds, and labels.
    let datasets = ["gdb13_1201", "gdb17_200", "checks", "coconut_55"];
    let bounds = [
        vec![],
        vec![Bound::Log],
        vec![Bound::Int],
        vec![Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
    ];
    let bound_strs = ["naive", "logbound", "intbound", "allbounds"];

    // Loop over all datasets of interest.
    for dataset in datasets.iter() {
        // Load all molecules from the given dataset.
        let paths = fs::read_dir(Path::new("data").join(dataset)).unwrap();
        let mut mol_list: Vec<Molecule> = Vec::new();
        for path in paths {
            let name = path.unwrap().path();
            if name.extension().and_then(OsStr::to_str) != Some("mol") {
                continue;
            }
            mol_list.push(
                loader::parse_molfile_str(
                    &fs::read_to_string(name.clone())
                        .expect(&format!("Could not read file {name:?}")),
                )
                .expect(&format!("Failed to parse {name:?}")),
            );
        }

        // For each of the bounds options, run the benchmark over all molecules
        // in this dataset.
        for (bound, bound_str) in zip(&bounds, &bound_strs) {
            group.bench_with_input(BenchmarkId::new(*dataset, &bound_str), bound, |b, bound| {
                b.iter(|| {
                    for mol in &mol_list {
                        index_search(&mol, &bound);
                    }
                });
            });
        }
    }

    group.finish();
}

criterion_group! {
    name = benchmark;
    config = Criterion::default().sample_size(20);
    targets = reference_datasets
}
criterion_main!(benchmark);
