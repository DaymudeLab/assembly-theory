use bit_set::BitSet;
use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use std::{ffi::OsStr, time::Instant};
use std::fs;
use std::path::Path;

use assembly_theory::{
    assembly::{clique_index_search_bench},
    loader,
    molecule::Molecule,
};

pub fn reference_datasets(c: &mut Criterion) {
    // Define a new criterion benchmark group of dataset benchmarks.
    let mut group = c.benchmark_group("reference_datasets");

    // Define datasets, bounds, and labels.
    let datasets = ["gdb13_1201", "gdb17_200", "coconut_55"];

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
        
        group.bench_function(*dataset, move |b| {
            b.iter_custom(|iters| {
                    let mut total = Duration::new(0, 0);
                    for _ in 0..iters {
                        for mol in mol_list.iter() {
                            let matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
                            let start = Instant::now();
                            clique_index_search_bench(mol, matches, assembly_theory::assembly::Kernel::Once);
                            total += start.elapsed()
                        }
                    }

                    total
                },
            );
        });
    }

    group.finish();
}

criterion_group! {
    name = benchmark;
    config = Criterion::default().sample_size(10);
    targets = reference_datasets
}
criterion_main!(benchmark);
