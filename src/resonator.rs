use bstr::{BStr, BString, ByteSlice};
use rayon::prelude::*;

use crate::encode::{encode, EncodedSeq};
use crate::error::ResonateError;
use crate::index::PartitionIndex;

/// An index for finding all reference sequences within a given Hamming distance
/// of a query sequence.
///
/// # Ordering
///
/// Query results are in **arbitrary order** (determined by the pigeonhole index).
///
/// # Thread safety
///
/// `HammingResonator` is `Send + Sync`; the index is read-only after construction.
#[derive(Debug)]
pub struct HammingResonator {
    originals: Vec<BString>,
    index: PartitionIndex,
}

impl HammingResonator {
    /// Build with explicit `max_dist`.
    pub fn with_max_dist(seqs: Vec<BString>, max_dist: u32) -> Result<Self, ResonateError> {
        let encoded = validate_and_encode(&seqs, max_dist)?;
        let index = PartitionIndex::build(encoded, max_dist)?;
        Ok(Self {
            originals: seqs,
            index,
        })
    }

    /// Query the index and return all references within `max_dist`.
    pub fn query(&self, query: &BStr) -> Result<Vec<&BStr>, ResonateError> {
        let enc = encode(query)?;
        if enc.len() != self.index.seq_len {
            return Err(ResonateError::QueryLengthMismatch {
                got: enc.len(),
                expected: self.index.seq_len,
            });
        }
        let indices = self.index.query_indices(&enc);
        Ok(indices
            .into_iter()
            .map(|(i, d)| self.originals[i as usize].as_bstr())
            .collect())
    }

    /// Query the index for a batch of queries in parallel using Rayon.
    #[must_use]
    pub fn query_batch(&self, queries: &[BString]) -> Vec<Result<Vec<&BStr>, ResonateError>> {
        queries.par_iter().map(|q| self.query(q.as_ref())).collect()
    }
}

#[derive(Debug)]
pub struct EncodedSeqs {
    pub arena: Vec<u8>,
    pub entry_len: usize,  //how many bytes of sequence are there
    pub entry_size: usize, //how many bytes per entry (for weighted resonator)
    pub count: usize,
}

impl EncodedSeqs {
    pub fn get_entry(&self, idx: u32) -> &[u8] {
        let start = (idx as usize) * self.entry_size;
        let end = start + self.entry_len;
        &self.arena[start..end]
    }

    pub fn iter(&self) -> impl Iterator<Item = &[u8]> {
        EncodedSeqsIter {
            arena: &self.arena,
            entry_len: &self.entry_len,
            entry_size: &self.entry_size,
            count: &self.count,
            index: 0,
        }
    }
}

struct EncodedSeqsIter<'a> {
    arena: &'a [u8],
    entry_len: &'a usize,
    entry_size: &'a usize,
    count: &'a usize,
    index: usize,
}
impl<'a> Iterator for EncodedSeqsIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= *self.count {
            return None;
        }
        let start = self.index * *self.entry_size;
        let end = start + *self.entry_len;
        self.index += 1;
        Some(&self.arena[start..end])
    }
}

/// Shared validation: checks non-empty, consistent length, valid bases, `seq_len` >= `max_dist` + 1.
/// Returns encoded sequences ready for `PartitionIndex::build`.
pub(crate) fn validate_and_encode(
    seqs: &[BString],
    max_dist: u32,
) -> Result<EncodedSeqs, ResonateError> {
    if seqs.is_empty() {
        return Err(ResonateError::EmptyInput);
    }

    let expected_len = seqs[0].len();
    let min_len = (max_dist as usize) + 1;
    if expected_len < min_len {
        return Err(ResonateError::SequenceTooShort {
            seq_len: expected_len,
            max_dist,
            min: min_len,
        });
    }

    let total_bytes = seqs.len() * expected_len;
    let mut arena = Vec::with_capacity(total_bytes);
    arena.resize(total_bytes, 0);

    for (i, seq) in seqs.iter().enumerate() {
        if seq.len() != expected_len {
            return Err(ResonateError::InconsistentLength(expected_len, seq.len()));
        }
        let start = i * expected_len;
        let end = start + expected_len;
        let slice = &mut arena[start..end];
        for (j, &b) in seq.iter().enumerate() {
            slice[j] = match b | 0x20 {
                b'a' => 0u8,
                b'c' => 1u8,
                b'g' => 2u8,
                b't' => 4u8,
                _ => return Err(ResonateError::InvalidBase(b as char, j)),
            };
        }
    }

    Ok(EncodedSeqs {
        arena,
        count: seqs.len(),
        entry_len: expected_len,
        entry_size: expected_len,
    })
}
/// Shared validation: checks non-empty, consistent length, valid bases, `seq_len` >= `max_dist` + 1.
/// Returns encoded sequences ready for `PartitionIndex::build`.
pub(crate) fn validate_and_encode_with_scores(
    seqs: &[BString],
    scores: &[f64],
    max_dist: u32,
) -> Result<EncodedSeqs, ResonateError> {
    if seqs.is_empty() {
        return Err(ResonateError::EmptyInput);
    }

    let expected_len = seqs[0].len();
    let min_len = (max_dist as usize) + 1;
    if expected_len < min_len {
        return Err(ResonateError::SequenceTooShort {
            seq_len: expected_len,
            max_dist,
            min: min_len,
        });
    }
    if seqs.len() != scores.len() {
        return Err(ResonateError::LengthMismatch);
    }

    let bytes_per_entry = expected_len + 8; //for f64

    let total_bytes = seqs.len() * bytes_per_entry;
    let mut arena = Vec::with_capacity(total_bytes);
    arena.resize(total_bytes, 0);

    for (i, (seq, score)) in seqs.iter().zip(scores.iter()).enumerate() {
        if seq.len() != expected_len {
            return Err(ResonateError::InconsistentLength(expected_len, seq.len()));
        }
        let start = i * bytes_per_entry;
        let end = start + expected_len;
        let slice = &mut arena[start..end];
        for (j, &b) in seq.iter().enumerate() {
            slice[j] = match b | 0x20 {
                b'a' => 0u8,
                b'c' => 1u8,
                b'g' => 2u8,
                b't' => 4u8,
                _ => return Err(ResonateError::InvalidBase(b as char, j)),
            };
        }
        arena[end..end + 8].copy_from_slice(&score.to_le_bytes());
    }

    Ok(EncodedSeqs {
        arena,
        count: seqs.len(),
        entry_len: expected_len,
        entry_size: bytes_per_entry,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::ByteSlice;

    fn bstr(s: &str) -> BString {
        BString::from(s)
    }

    fn resonator(seqs: &[&str], max_dist: u32) -> HammingResonator {
        let v: Vec<BString> = seqs.iter().map(|&s| bstr(s)).collect();
        HammingResonator::with_max_dist(v, max_dist).unwrap()
    }

    #[test]
    fn exact_match_returned() {
        let r = resonator(&["ACGT", "TTTT"], 1);
        let hits = r.query("ACGT".as_bytes().as_bstr()).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0], "ACGT".as_bytes().as_bstr());
    }

    #[test]
    fn d1_match_found() {
        let r = resonator(&["AAAA", "AAAC"], 1);
        let hits = r.query("AAAA".as_bytes().as_bstr()).unwrap();
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn d2_not_returned_at_max_dist_1() {
        let r = resonator(&["AAAA", "AACC"], 1);
        let hits = r.query("AAAA".as_bytes().as_bstr()).unwrap();
        assert_eq!(hits.len(), 1); // only exact match
    }

    #[test]
    fn no_matches() {
        let r = resonator(&["GGGG"], 1);
        let hits = r.query("AAAA".as_bytes().as_bstr()).unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn invalid_base_in_query() {
        let r = resonator(&["ACGT"], 1);
        let err = r.query("ACGN".as_bytes().as_bstr()).unwrap_err();
        assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
    }

    #[test]
    fn mismatched_query_length() {
        let r = resonator(&["ACGT"], 1);
        let err = r.query("ACG".as_bytes().as_bstr()).unwrap_err();
        assert!(matches!(
            err,
            ResonateError::QueryLengthMismatch {
                got: 3,
                expected: 4
            }
        ));
    }

    #[test]
    fn query_batch_matches_serial() {
        let r = resonator(&["AAAA", "AAAC", "AAGC"], 1);
        let queries: Vec<BString> = vec![bstr("AAAA"), bstr("AAAC"), bstr("TTTT")];
        let batch = r.query_batch(&queries);
        let serial: Vec<_> = queries.iter().map(|q| r.query(q.as_ref())).collect();
        assert_eq!(batch.len(), serial.len());
        for (b, s) in batch.iter().zip(serial.iter()) {
            let mut bv: Vec<_> = b.as_ref().unwrap().iter().copied().collect();
            let mut sv: Vec<_> = s.as_ref().unwrap().iter().copied().collect();
            bv.sort();
            sv.sort();
            assert_eq!(bv, sv);
        }
    }
}
