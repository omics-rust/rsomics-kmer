use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use rsomics_kmer::{RollingKmers, canonical, reverse_complement};

fn fixture(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| if i % 4_093 == 0 { b'N' } else { b"ACGT"[i % 4] })
        .collect()
}

fn rolling(c: &mut Criterion) {
    let sequence = fixture(1_048_576);
    let mut group = c.benchmark_group("rolling");
    group.throughput(Throughput::Bytes(sequence.len() as u64));
    group.bench_function("k31_1m", |b| {
        b.iter(|| {
            RollingKmers::new(black_box(sequence.as_slice()), 31)
                .flatten()
                .fold(0, |acc, kmer| acc ^ kmer)
        });
    });
    group.finish();
}

fn strand_operations(c: &mut Criterion) {
    let sequence = fixture(65_536);
    let kmers: Vec<_> = RollingKmers::new(&sequence, 31).flatten().collect();
    let mut group = c.benchmark_group("strand_operations");
    group.throughput(Throughput::Elements(kmers.len() as u64));
    group.bench_function("reverse_complement_k31", |b| {
        b.iter(|| {
            black_box(&kmers)
                .iter()
                .fold(0, |acc, &kmer| acc ^ reverse_complement(kmer, 31))
        });
    });
    group.bench_function("canonical_k31", |b| {
        b.iter(|| {
            black_box(&kmers)
                .iter()
                .fold(0, |acc, &kmer| acc ^ canonical(kmer, 31))
        });
    });
    group.finish();
}

criterion_group!(benches, rolling, strand_operations);
criterion_main!(benches);
