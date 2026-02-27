# Implementation Plan: `hamming-resonate` crate

## API

```rust
// Unweighted: return all references within max_dist
let resonator = HammingResonator::new(vec_of_bstrings)?;
resonator.query(query_bstring) -> Result<Vec<&BStr>, ResonateError>

// Weighted: return the single highest-scoring reference within max_dist
let resonator = HammingResonatorWeighted::new(vec_of_bstrings_and_scores)?;
resonator.query_best(query_bstring) -> Result<Option<&BStr>, ResonateError>
```

---

## Crate layout

```
hamming-resonate/
├── Cargo.toml
├── README.md
└── src/
    ├── lib.rs          # public API re-exports
    ├── encode.rs       # DNA 2-bit encoding + hamming distance
    ├── index.rs        # pigeonhole partition index (shared core)
    ├── resonator.rs    # HammingResonator
    └── weighted.rs     # HammingResonatorWeighted
```

---

## `Cargo.toml`

```toml
[package]
name = "hamming-resonate"
version = "0.1.0"
edition = "2021"

[dependencies]
bstr = "1"
rayon = "1"
thiserror = "1"
bitnuc="1"
hamming-bitwise-fast = ""

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }

[[bench]]
name = "throughput"
harness = false
```

---

## Algorithm: Pigeonhole Partition Index

For a maximum Hamming distance `d`, split each sequence into `d + 1` chunks. By
the pigeonhole principle, if two sequences differ in at most `d` positions,
**at least one chunk must match exactly**. This allows us to:

1. Build an inverted index: `chunk_bytes → [ref_indices]` for each chunk position
2. At query time, look up each of the `d + 1` chunks, union the candidates, and verify with a full Hamming check

This gives **O(N·L)** build time and **O(d·L + candidates·L)** per query, where candidates ≈ `N / 4^(L/(d+1))` — typically a few hundred out of 1M for typical barcode lengths.

---

## Module: `encode.rs`

Encode `&BStr` into a `Vec<u8>` (one byte per base, 2-bit values) for fast comparison.

- `fn encode(seq: &BStr) -> Result<Vec<u8>, ResonateError>` — maps A→0, C→1, G→2, T→4, errors on anything else
- `#[inline] fn hamming(a: &[u8], b: &[u8]) -> u32` — use hamming-bitwisefast.

> One byte per base (not bit-packed) is intentional: bit-packing complicates chunk slicing with no benefit until L > 64, and the byte layout allows SIMD vectorization of the Hamming loop.

---

## Module: `index.rs` — `PartitionIndex`

The shared engine used by both public structs.

```rust
pub(crate) struct PartitionIndex {
    /// Encoded reference sequences (one byte per base)
    encoded: Vec<Vec<u8>>,
    /// chunk_maps[i]: chunk_bytes -> list of ref indices
    chunk_maps: Vec<HashMap<Vec<u8>, Vec<u32>>>,
    /// Byte range (start, end) for each chunk
    ranges: Vec<(usize, usize)>,
    pub(crate) seq_len: usize,
    pub(crate) max_dist: u32,
}
```

### Construction

```rust
pub(crate) fn build(encoded: Vec<Vec<u8>>, max_dist: u32) -> Result<Self, ResonateError>
```

- Number of chunks = `max_dist + 1`
- Chunk boundaries: `chunk_i = [i*L/(d+1) .. (i+1)*L/(d+1)]`
- Single pass over references: for each sequence and each chunk, insert into that chunk's `HashMap`

### Query

```rust
pub(crate) fn query_indices(&self, enc: &[u8]) -> Vec<u32>
```

- For each chunk, look up `enc[start..end]` in that chunk's map and collect candidate ref indices
- `sort_unstable` + `dedup` the merged candidates
- Filter: keep only indices where `hamming(enc, &self.encoded[idx]) <= self.max_dist`
- Return surviving indices

---

## Module: `resonator.rs` — `HammingResonator`

```rust
pub struct HammingResonator {
    originals: Vec<BString>,   // kept for returning &BStr references
    index: PartitionIndex,
}
```

### Constructors

```rust
impl HammingResonator {
    /// Build with default max_dist = 1
    pub fn new(seqs: Vec<BString>) -> Result<Self, ResonateError>

    /// Build with explicit max_dist
    pub fn with_max_dist(seqs: Vec<BString>, max_dist: u32) -> Result<Self, ResonateError>
}
```

Validation (shared with `HammingResonatorWeighted`):
- Non-empty input
- All sequences same length
- All sequences contain only valid bases (A/C/G/T, case-insensitive)
- `seq_len >= max_dist + 1` (need at least one base per chunk)

### Query

```rust
pub fn query(&self, query: &BStr) -> Result<Vec<&BStr>, ResonateError>
```

- Encode query, call `index.query_indices(enc)`
- Map indices back to `&self.originals[i].as_bstr()`
- Results are in **arbitrary order** (document this)

### Batch query

```rust
pub fn query_batch(&self, queries: &[BString]) -> Vec<Result<Vec<&BStr>, ResonateError>>
```

- Uses `rayon::par_iter()` — safe because `PartitionIndex` is read-only after build
- `HammingResonator` should be `Send + Sync` (document this explicitly)

---

## Module: `weighted.rs` — `HammingResonatorWeighted`

```rust
pub struct HammingResonatorWeighted {
    originals: Vec<BString>,
    scores: Vec<f64>,          // parallel to originals
    index: PartitionIndex,
}
```

### Constructors

```rust
impl HammingResonatorWeighted {
    /// Build with default max_dist = 1
    pub fn new(seqs: Vec<(BString, f64)>) -> Result<Self, ResonateError>

    /// Build with explicit max_dist
    pub fn with_max_dist(seqs: Vec<(BString, f64)>, max_dist: u32) -> Result<Self, ResonateError>
}
```

### Query

```rust
pub fn query_best(&self, query: &BStr) -> Result<Option<&BStr>, ResonateError>
```

- Call `index.query_indices(enc)` → candidate indices
- Return `originals[i]` where `scores[i]` is highest among candidates
- Tie-breaking: lowest original index (stable; document this)
- Returns `Ok(None)` if no candidates found within `max_dist`

---

## Error type

```rust
// src/error.rs (or inline in lib.rs)

#[derive(thiserror::Error, Debug)]
pub enum ResonateError {
    #[error("Empty reference set")]
    EmptyInput,

    #[error("All sequences must have the same length, got {0} and {1}")]
    InconsistentLength(usize, usize),

    #[error("Invalid base '{0}' at position {1}")]
    InvalidBase(char, usize),

    #[error("Sequence length {seq_len} too short for max_dist={max_dist} (need at least {min} bases)")]
    SequenceTooShort { seq_len: usize, max_dist: u32, min: usize },

    #[error("Query length {got} does not match index length {expected}")]
    QueryLengthMismatch { got: usize, expected: usize },

    #[error("Duplicate entries in reference set")]
    DuplicateEntries
}
```

---

## `lib.rs`

```rust
mod encode;
mod error;
mod index;
mod resonator;
mod weighted;

pub use bstr::{BStr, BString};
pub use error::ResonateError;
pub use resonator::HammingResonator;
pub use weighted::HammingResonatorWeighted;
```

---

## Tests

Write inline unit tests in each module and integration tests in `tests/integration.rs`.

| # | What to test |
|---|---|
| 1 | Exact match (distance 0) is always returned |
| 2 | d=1 match found; sequence at d=2 not returned when max_dist=1 |
| 3 | d=2 and d=3 matches found at correct positions |
| 4 | No matches → `Ok(vec![])` / `Ok(None)` |
| 5 | Invalid base in reference → `Err(InvalidBase(...))` |
| 6 | Invalid base in query → `Err(InvalidBase(...))` |
| 7 | Mismatched query length → `Err(QueryLengthMismatch {...})` |
| 8 | `query_best` returns highest-score match, not just first or closest |
| 9 | `query_best` tie-break: lowest index wins |
| 10 | All references identical — dedup works, single result returned |
| 11 | Case-insensitive bases (a/c/g/t == A/C/G/T) |
| 12 | `query_batch` results match serial `query` results |

---

## Benchmarks (`benches/throughput.rs`)

Use `criterion`. Two groups:

**Build**
- 1M sequences of length 16, `max_dist = 1 / 2 / 3`
- Measure `HammingResonator::with_max_dist(...)` wall time

**Query**
- Pre-built index of 1M seqs
- 10k queries, `max_dist = 1 / 2 / 3`
- Measure `query()` and `query_batch()` (parallel) separately

---

## Implementation order for Claude Code

1. `encode.rs` — encoding + Hamming function + unit tests
2. `index.rs` — `PartitionIndex` build + `query_indices` + unit tests
3. `error.rs` — `ResonateError`
4. `resonator.rs` — `HammingResonator` + integration tests
5. `weighted.rs` — `HammingResonatorWeighted` + integration tests
6. `lib.rs` — re-exports
7. `benches/throughput.rs` — criterion benchmarks
8. `README.md` — usage examples

---

## Open questions

| Question | Recommendation |
|---|---|
| Default `max_dist` | 1 (most common use case; override with `with_max_dist`) |
| Score type | `f64` for now; consider a generic `T: PartialOrd + Copy` in a follow-up |
| Alphabet | DNA only (ACGT, case-insensitive); arbitrary bytes would need a different chunking strategy |
| Case handling | Normalize to uppercase on encode; do not error on lowercase |
