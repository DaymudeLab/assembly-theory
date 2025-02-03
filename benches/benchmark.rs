use std::path::PathBuf;

use fastassembly::assembly::index;

use criterion::{criterion_group, criterion_main, Criterion};
use fastassembly::loader;

pub fn criterion_benchmark(c: &mut Criterion) {
    let molecule = loader::parse(&PathBuf::from("data/aspartic.sdf")).unwrap();
    c.bench_function("aspartic", |b| b.iter(|| index(&molecule)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
