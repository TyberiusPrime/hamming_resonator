use bstr::BString;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use hamming_resonate::HammingResonator;

/// Generate N random-ish DNA sequences of the given length.
/// Uses a simple LCG for determinism without pulling in rand.
fn gen_seqs(n: usize, len: usize, seed: u64) -> Vec<BString> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut state = seed;
    (0..n)
        .map(|_| {
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    BASES[((state >> 33) & 3) as usize]
                })
                .collect();
            BString::from(seq)
        })
        .collect()
}

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");
    let seqs_1m = gen_seqs(1_000_000, 16, 42);
    group.throughput(Throughput::Elements(1_000_000));

    for max_dist in [1u32, 2, 3] {
        group.bench_with_input(
            BenchmarkId::new("1M_seqs_len16", max_dist),
            &max_dist,
            |b, &d| {
                b.iter(|| {
                    HammingResonator::with_max_dist(seqs_1m.clone(), d).unwrap()
                });
            },
        );
    }
    group.finish();
}

fn bench_query(c: &mut Criterion) {
    let seqs_1m = gen_seqs(1_000_000, 16, 42);
    let queries = gen_seqs(10_000, 16, 99);

    let mut group = c.benchmark_group("query");
    group.throughput(Throughput::Elements(10_000));

    for max_dist in [1u32, 2, 3] {
        let resonator = HammingResonator::with_max_dist(seqs_1m.clone(), max_dist).unwrap();

        group.bench_with_input(
            BenchmarkId::new("serial_10k", max_dist),
            &max_dist,
            |b, _| {
                b.iter(|| {
                    for q in &queries {
                        let _ = resonator.query(q.as_ref()).unwrap();
                    }
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("batch_parallel_10k", max_dist),
            &max_dist,
            |b, _| {
                b.iter(|| resonator.query_batch(&queries));
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_build, bench_query);
criterion_main!(benches);
