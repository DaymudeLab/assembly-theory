use std::path::PathBuf;

use orca::assembly::index;

use criterion::{criterion_group, criterion_main, Criterion};
use orca::loader;

pub fn criterion_benchmark(c: &mut Criterion) {
    for str in ["aspartic", ""] {
        let molecule = loader::parse(&PathBuf::from(format!("data/{str}.sdf"))).unwrap();
        c.bench_function("aspartic", |b| b.iter(|| index(&molecule)));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
