use bstr::{BStr, BString, ByteSlice};
use rayon::prelude::*;

use crate::encode::EncodedSeqs;
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
    index: PartitionIndex<EncodedSeqs>,
}

impl HammingResonator {
    /// Build with explicit `max_dist`.
    pub fn with_max_dist(seqs: Vec<BString>, max_dist: u32) -> Result<Self, ResonateError> {
        let encoded = EncodedSeqs::new(&seqs, max_dist)?;
        let index = PartitionIndex::build(encoded, max_dist)?;
        Ok(Self {
            originals: seqs,
            index,
        })
    }

    /// Query the index and return all references within `max_dist`, and their distanc.
    pub fn query(&self, query: &BStr) -> Result<Vec<(&BStr, u32)>, ResonateError> {
        if query.len() != self.index.seq_len {
            return Err(ResonateError::QueryLengthMismatch {
                got: query.len(),
                expected: self.index.seq_len,
            });
        }
        let query = BString::from(query.to_ascii_uppercase());
        let indices = self.index.query_indices(query.as_bstr());
        Ok(indices
            .into_iter()
            .map(|(i, d, _ignored_score)| (self.originals[i as usize].as_bstr(), d))
            .collect())
    }

    /// Query the index for a batch of queries in parallel using Rayon.
    #[must_use]
    pub fn query_batch(
        &self,
        queries: &[BString],
    ) -> Vec<Result<Vec<(&BStr, u32)>, ResonateError>> {
        queries.par_iter().map(|q| self.query(q.as_ref())).collect()
    }
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
        assert_eq!(hits[0].0, "ACGT".as_bytes().as_bstr());
        assert_eq!(hits[0].1, 0);
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
