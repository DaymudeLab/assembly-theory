use std::fs;
use std::path::PathBuf;

use assembly_theory::assembly::{index, naive_index_search};

use assembly_theory::loader;
use criterion::{criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(c: &mut Criterion) {
    for str in ["aspartic", "benzene", "aspirin", "morphine"] {
        let path = PathBuf::from(format!("./data/checks/{str}.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(str, |b| b.iter(|| index(&molecule)));
    }

    for str in ["aspartic", "benzene", "aspirin"] {
        let path = PathBuf::from(format!("./data/checks/{str}.mol"));
        let molfile = fs::read_to_string(path).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(&format!("naive-{str}"), |b| {
            b.iter(|| naive_index_search(&molecule))
        });
    }
}

pub fn gdb13_benchmark(c: &mut Criterion) {
    let paths = fs::read_dir("data/gdb13_1201").unwrap();

    for path in paths {
        let name = path.unwrap().path();
        let molfile = fs::read_to_string(name.clone()).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(name.to_str().unwrap(), |b| b.iter(|| index(&molecule)));
    }
}

criterion_group!(benches, criterion_benchmark, gdb13_benchmark);
criterion_main!(benches);
