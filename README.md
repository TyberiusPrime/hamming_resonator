# hamming-resonate

Fast nearest-neighbors search over fixed-length byte sequences using Hamming distance.

For DNA usage.

Build an index / database from a set of reference sequences, then query it to find all
references within a given number of mismatches.  

Two variants are provided:

- **`HammingResonator`** — returns every reference within `max_dist`.
- **`HammingResonatorWeighted`** — returns the single best sequence within `max_dist`. 
  On distance ties, highest score wins. On score ties, lowest-index wins.

Sequences are byte strings. Comparisons are case sensitive.
All sequences in one index must be the same length. Query must also have the same length.



## Limitations

Tests have only been performed up to `max_dist` = 3.

Database size is at most 2^32 sequences.


## Quick start

```rust
use bstr::{BString, ByteSlice};
use hamming_resonate::{HammingResonator, ResonateError};

fn main() -> Result<(), ResonateError> {
    let refs: Vec<BString> = ["AAAA", "AAAC", "TTTT", "CCCC"]
        .iter().map(|&s| BString::from(s)).collect();

    let index = HammingResonator::with_max_dist(refs, 1)?;

    let mut hits: Vec<_> = index.query("AAAG".as_bytes().as_bstr())?
        .into_iter().map(|(h, _distance)| h.to_str().unwrap()).collect();
    hits.sort();
    // "AAAA" (1 mismatch at pos 3) and "AAAC" (0 mismatch) are within distance 1.
    assert_eq!(hits, ["AAAA", "AAAC"]);
    Ok(())
}
```

### Weighted variant

```rust
use bstr::{BString, ByteSlice};
use hamming_resonate::{HammingResonatorWeighted, ResonateError};

fn main() -> Result<(), ResonateError> {
    let refs: Vec<(BString, f32)> = vec![
        (BString::from("AAAT"), 9.0),
        (BString::from("AAAG"), 5.0),
    ];

    let index = HammingResonatorWeighted::with_max_dist(refs, 1)?;

    // Both "AAAA" (d=0) and "AAAG" (d=1) are within distance 1 of "AAAA",
    // but "AAAT" has a higher score so it wins.
    let best = index.query_best("AAAA".as_bytes().as_bstr())?;
    assert_eq!(best.unwrap(), "AAAT".as_bytes().as_bstr());

    // AAAG is distance 0 from AAAG, so it is returned
    let best = index.query_best("AAAG".as_bytes().as_bstr())?;
    assert_eq!(best.unwrap(), "AAAG".as_bytes().as_bstr());

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

    assert!(results[0].as_ref().unwrap().iter().map(|(seq, dist)| seq).any(|seq| seq == &"AAAC".as_bytes().as_bstr()));
    assert!(results[1].as_ref().unwrap().iter().map(|(seq, dist)| seq).any(|seq| seq == &"TTTT".as_bytes().as_bstr()));
    Ok(())
}
```


## Performance notes

Sorting your references sequences gives about a 7-10% speed boost in query performance.

But since it increases the build time substantially, it's not on by default.

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

## License

Licensed under the MIT License.
