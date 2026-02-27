/// Benchmark comparing two Hamming distance implementations for DNA sequences.
///
/// Method A — **byte_cmp** (current):
///   Encode one byte per base (A=0, C=1, G=2, T=4).
///   Distance = count of positions where bytes differ.
///
/// Method B — **nibble_fast** (proposed):
///   Encode two bases per byte using 4-bit one-hot nibbles
///   (A=0x1, C=0x2, G=0x4, T=0x8; even index → high nibble, odd → low nibble).
///   Any two different bases XOR to exactly 2 bits in a nibble, so
///   `hamming_bitwise_fast(a, b) / 2` gives the correct base-level distance.
///   This halves the memory footprint and lets the bitwise path use wider SIMD words.
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use hamming_bitwise_fast::hamming_bitwise_fast;

// ---------------------------------------------------------------------------
// Encoders (self-contained; don't depend on pub(crate) internals)
// ---------------------------------------------------------------------------

/// Encode DNA bytes (ACGT, any case) to one byte per base: A=0, C=1, G=2, T=4.
fn encode_byte_per_base(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b | 0x20 {
            b'a' => 0u8,
            b'c' => 1u8,
            b'g' => 2u8,
            b't' => 4u8,
            _ => panic!("invalid base {b}"),
        })
        .collect()
}

/// Encode DNA bytes to nibble-packed form: 2 bases per byte, one-hot nibbles.
/// Even-indexed base → high nibble (bits 7:4), odd-indexed → low nibble (bits 3:0).
/// One-hot nibble values: A=0x1, C=0x2, G=0x4, T=0x8.
/// A sequence of length L produces ceil(L/2) bytes.
fn encode_nibble_packed(seq: &[u8]) -> Vec<u8> {
    let packed_len = seq.len().div_ceil(2);
    let mut out = vec![0u8; packed_len];
    for (i, &b) in seq.iter().enumerate() {
        let nibble: u8 = match b | 0x20 {
            b'a' => 0x1,
            b'c' => 0x2,
            b'g' => 0x4,
            b't' => 0x8,
            _ => panic!("invalid base {b}"),
        };
        if i % 2 == 0 {
            out[i / 2] |= nibble << 4;
        } else {
            out[i / 2] |= nibble;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Distance functions
// ---------------------------------------------------------------------------

/// Method A: byte comparison loop (current implementation).
#[inline]
fn hamming_byte_cmp(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b.iter()).map(|(x, y)| (x != y) as u32).sum()
}

/// Method B: nibble-packed + hamming_bitwise_fast / 2.
/// Correct because each one-hot nibble pair XORs to exactly 2 bits when bases differ,
/// so bitwise Hamming = 2 * base-level Hamming.
#[inline]
fn hamming_nibble_fast(a: &[u8], b: &[u8]) -> u32 {
    hamming_bitwise_fast(a, b) / 2
}

// ---------------------------------------------------------------------------
// Test data generation
// ---------------------------------------------------------------------------

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn lcg_next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

fn gen_dna(n: usize, len: usize, seed: u64) -> Vec<Vec<u8>> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            (0..len)
                .map(|_| BASES[((lcg_next(&mut state) >> 33) & 3) as usize])
                .collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Benchmarks
// ---------------------------------------------------------------------------

fn bench_hamming_methods(c: &mut Criterion) {
    let n_pairs = 10_000usize;

    for seq_len in [16usize, 64, 256] {
        let raw = gen_dna(n_pairs * 2, seq_len, 12345);

        // Pre-encode into both formats so encoding cost is excluded.
        let byte_encoded: Vec<Vec<u8>> = raw.iter().map(|s| encode_byte_per_base(s)).collect();
        let nibble_encoded: Vec<Vec<u8>> = raw.iter().map(|s| encode_nibble_packed(s)).collect();

        // Sanity check: both methods agree on every pair.
        for i in (0..n_pairs * 2).step_by(2) {
            let expected = hamming_byte_cmp(&byte_encoded[i], &byte_encoded[i + 1]);
            let got = hamming_nibble_fast(&nibble_encoded[i], &nibble_encoded[i + 1]);
            assert_eq!(expected, got, "mismatch at pair {i}");
        }

        let mut group = c.benchmark_group(format!("hamming_len{seq_len}"));
        group.throughput(Throughput::Elements(n_pairs as u64));

        group.bench_with_input(
            BenchmarkId::new("byte_cmp", seq_len),
            &seq_len,
            |b, _| {
                b.iter(|| {
                    let mut total = 0u32;
                    for i in (0..n_pairs * 2).step_by(2) {
                        total = total
                            .wrapping_add(hamming_byte_cmp(&byte_encoded[i], &byte_encoded[i + 1]));
                    }
                    total
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("nibble_fast", seq_len),
            &seq_len,
            |b, _| {
                b.iter(|| {
                    let mut total = 0u32;
                    for i in (0..n_pairs * 2).step_by(2) {
                        total = total.wrapping_add(hamming_nibble_fast(
                            &nibble_encoded[i],
                            &nibble_encoded[i + 1],
                        ));
                    }
                    total
                });
            },
        );

        group.finish();
    }
}

criterion_group!(benches, bench_hamming_methods);
criterion_main!(benches);
