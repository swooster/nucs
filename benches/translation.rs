use criterion::{BatchSize, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{Rng, SeedableRng, rngs::StdRng};

use nucs::{AmbiNuc, DnaSlice, NCBI1, NCBI1_RC, Nuc};

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = StdRng::from_os_rng();

    let lens = [16, 256, 1024];

    let mut g = c.benchmark_group("Nuc translation");
    for len in lens {
        g.throughput(Throughput::Elements(len));
        g.bench_function(BenchmarkId::from_parameter(len), |b| {
            b.iter_batched_ref(
                || (0..len).map(|_| rng.random()).collect::<Vec<Nuc>>(),
                |dna| dna.translate_to_vec(NCBI1),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();

    let mut g = c.benchmark_group("AmbiNuc translation");
    for len in lens {
        g.throughput(Throughput::Elements(len));
        g.bench_function(BenchmarkId::from_parameter(len), |b| {
            b.iter_batched_ref(
                || (0..len).map(|_| rng.random()).collect::<Vec<AmbiNuc>>(),
                |dna| dna.translate_to_vec(NCBI1),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();

    let mut g = c.benchmark_group("Nuc RC translation");
    for len in lens {
        g.throughput(Throughput::Elements(len));
        g.bench_function(BenchmarkId::from_parameter(len), |b| {
            b.iter_batched_ref(
                || (0..len).map(|_| rng.random()).collect::<Vec<Nuc>>(),
                |dna| dna.rev_translate_to_vec(NCBI1_RC),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();

    let mut g = c.benchmark_group("AmbiNuc RC translation");
    for len in lens {
        g.throughput(Throughput::Elements(len));
        g.bench_function(BenchmarkId::from_parameter(len), |b| {
            b.iter_batched_ref(
                || (0..len).map(|_| rng.random()).collect::<Vec<AmbiNuc>>(),
                |dna| dna.rev_translate_to_vec(NCBI1_RC),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
