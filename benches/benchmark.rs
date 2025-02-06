use std::fs;
use std::path::PathBuf;

use orca::assembly::{index, naive_index};

use criterion::{criterion_group, criterion_main, Criterion};
use orca::loader;

pub fn criterion_benchmark(c: &mut Criterion) {
    for str in ["aspartic", "benzene", "cubane", "aspirin", "morphine"] {
        let path = PathBuf::from(format!("./data/checks/{str}.sdf"));
        let molfile = fs::read_to_string(path).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(str, |b| b.iter(|| index(&molecule)));
    }

    for str in ["aspartic", "benzene", "aspirin"] {
        let path = PathBuf::from(format!("./data/checks/{str}.sdf"));
        let molfile = fs::read_to_string(path).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(&format!("naive-{str}"), |b| {
            b.iter(|| naive_index(&molecule))
        });
    }
}

pub fn gdb13_benchmark(c: &mut Criterion) {
    let paths = fs::read_dir("data/gdb13").unwrap();

    for path in paths {
        let name = path.unwrap().path();
        let molfile = fs::read_to_string(name.clone()).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        c.bench_function(name.to_str().unwrap(), |b| b.iter(|| index(&molecule)));
    }
}

criterion_group!(benches, criterion_benchmark, gdb13_benchmark);
criterion_main!(benches);
