use std::{ffi::OsStr, fs, path::Path};

use criterion::{criterion_group, criterion_main, Criterion};

use assembly_theory::{assembly::index, loader::parse_molfile_str, molecule::Molecule};

/// Parse all .mol files in the given dataset as [`Molecule`]s.
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

/// Benchmark the [`index`] function for different reference datasets.
pub fn joss_bench(c: &mut Criterion) {
    // Define a new criterion benchmark group for the JOSS manuscript.
    let mut bench_group = c.benchmark_group("joss");

    // Run the benchmark for each reference dataset.
    for dataset in ["gdb13_1201", "gdb17_200", "checks", "coconut_55"] {
        let mol_list = load_dataset_molecules(dataset);
        bench_group.bench_function(dataset, |b| {
            b.iter(|| {
                for mol in &mol_list {
                    index(&mol);
                }
            });
        });
    }

    bench_group.finish();
}

criterion_group! {
    name = benchmark;
    config = Criterion::default().sample_size(20);
    targets = joss_bench
}
criterion_main!(benchmark);
