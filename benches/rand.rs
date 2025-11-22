use criterion::{Criterion, criterion_group, criterion_main};
use rand::{Rng, SeedableRng, rngs::StdRng};

use nucs::{AmbiAmino, AmbiNuc, Amino, Nuc};

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = StdRng::from_os_rng();
    c.bench_function("random Nuc", |b| b.iter(|| rng.random::<Nuc>()));
    c.bench_function("random AmbiNuc", |b| b.iter(|| rng.random::<AmbiNuc>()));
    c.bench_function("random Amino", |b| b.iter(|| rng.random::<Amino>()));
    c.bench_function("random AmbiAmino", |b| b.iter(|| rng.random::<AmbiAmino>()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
