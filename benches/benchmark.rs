use std::path::PathBuf;

use orca::assembly::{index, naive_index};

use criterion::{criterion_group, criterion_main, Criterion};
use orca::loader;

pub fn criterion_benchmark(c: &mut Criterion) {
    for str in ["aspartic", "benzene", "cubane", "aspirin", "morphine"] {
        let molecule = loader::parse(&PathBuf::from(format!("data/checks/{str}.sdf"))).unwrap();
        c.bench_function(str, |b| b.iter(|| index(&molecule)));
    }

    for str in ["aspartic", "benzene", "cubane", "aspirin"] {
        let molecule = loader::parse(&PathBuf::from(format!("data/checks/{str}.sdf"))).unwrap();
        c.bench_function(&format!("naive-{str}"), |b| {
            b.iter(|| naive_index(&molecule))
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
