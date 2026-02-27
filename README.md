# hamming-resonate

Fast nearest-neighbour search over fixed-length DNA sequences using Hamming distance.

Build an index from a set of reference sequences, then query it to find all
references within a given number of mismatches.  Two flavours are provided:

- **`HammingResonator`** — returns every reference within `max_dist`.
- **`HammingResonatorWeighted`** — returns the single highest-scoring reference
  within `max_dist`, with lowest-index tie-breaking.

Sequences are ASCII strings over `{A, C, G, T}` (case-insensitive).  All
sequences in one index must be the same length.

## Quick start

```rust
use bstr::{BString, ByteSlice};
use hamming_resonate::{HammingResonator, ResonateError};

fn main() -> Result<(), ResonateError> {
    let refs: Vec<BString> = ["AAAA", "AAAC", "TTTT", "CCCC"]
        .iter().map(|&s| BString::from(s)).collect();

    let index = HammingResonator::with_max_dist(refs, 1)?;

    let mut hits: Vec<_> = index.query("AAAG".as_bytes().as_bstr())?
        .into_iter().map(|h| h.to_str().unwrap()).collect();
    hits.sort();
    // "AAAA" (1 mismatch at pos 3) and "AAAC" (1 mismatch at pos 3) are within distance 1.
    assert_eq!(hits, ["AAAA", "AAAC"]);
    Ok(())
}
```

### Weighted variant

```rust
use bstr::{BString, ByteSlice};
use hamming_resonate::{HammingResonatorWeighted, ResonateError};

fn main() -> Result<(), ResonateError> {
    let refs: Vec<(BString, f64)> = vec![
        (BString::from("AAAA"), 1.0),
        (BString::from("AAAC"), 9.0),
    ];

    let index = HammingResonatorWeighted::with_max_dist(refs, 1)?;

    // Both "AAAA" (d=0) and "AAAC" (d=1) are within distance 1 of "AAAA",
    // but "AAAC" has a higher score so it wins.
    let best = index.query_best("AAAA".as_bytes().as_bstr())?;
    assert_eq!(best.unwrap(), "AAAC".as_bytes().as_bstr());
    Ok(())
}
```

### Parallel batch queries

`HammingResonator` exposes `query_batch`, which runs queries in parallel via
[Rayon](https://docs.rs/rayon):

```rust
use bstr::{BString, ByteSlice};
use hamming_resonate::{HammingResonator, ResonateError};

fn main() -> Result<(), ResonateError> {
    let refs: Vec<BString> = ["AAAA", "AAAC", "TTTT"]
        .iter().map(|&s| BString::from(s)).collect();
    let index = HammingResonator::with_max_dist(refs, 1)?;

    let queries: Vec<BString> = ["AAAA", "TTTT"]
        .iter().map(|&s| BString::from(s)).collect();
    let results = index.query_batch(&queries);

    assert!(results[0].as_ref().unwrap().contains(&"AAAC".as_bytes().as_bstr()));
    assert!(results[1].as_ref().unwrap().contains(&"TTTT".as_bytes().as_bstr()));
    Ok(())
}
```

## API

### `HammingResonator`

| Method | Description |
|---|---|
| `with_max_dist(seqs, max_dist)` | Build the index. |
| `query(query)` | All references within `max_dist` mismatches, in arbitrary order. |
| `query_batch(queries)` | Parallel version of `query` over a slice of queries. |

### `HammingResonatorWeighted`

| Method | Description |
|---|---|
| `with_max_dist(seqs, max_dist)` | Build the index from `(sequence, score)` pairs. |
| `query_best(query)` | Highest-scoring reference within `max_dist`, or `None`. Ties broken by lowest original index. |

### Errors

All fallible operations return `ResonateError`:

| Variant | Cause |
|---|---|
| `EmptyInput` | Reference set is empty. |
| `InconsistentLength(expected, got)` | Sequences have differing lengths. |
| `InvalidBase(char, pos)` | A base other than A/C/G/T was encountered. |
| `SequenceTooShort { seq_len, max_dist, min }` | Sequence too short for the requested `max_dist` (need at least `max_dist + 1` bases). |
| `QueryLengthMismatch { got, expected }` | Query length differs from the indexed sequence length. |

## Algorithm

The index uses a **pigeonhole partition** strategy.  For `max_dist = d`, every
sequence is split into `d + 1` non-overlapping chunks.  By the pigeonhole
principle, any pair within Hamming distance `d` must share at least one chunk
exactly.  Lookup therefore:

1. Collects candidate indices by matching each query chunk against a per-chunk
   `HashMap`.
2. Deduplicates and verifies each candidate with an exact Hamming distance check.

Distance checks use a **nibble-packed** representation: each base occupies one
nibble with a one-hot encoding (`A=0x1`, `C=0x2`, `G=0x4`, `T=0x8`).  Any two
distinct bases XOR to exactly 2 set bits, so
`popcount(a XOR b) / 2` gives the base-level Hamming distance, computed by
[`hamming-bitwise-fast`](https://docs.rs/hamming-bitwise-fast) using SIMD
population-count intrinsics where available.

## License

Licensed under the MIT License.
