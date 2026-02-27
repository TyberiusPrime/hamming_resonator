use bstr::{BStr, BString};

use crate::encode::encode;
use crate::error::ResonateError;
use crate::index::PartitionIndex;
use crate::resonator::validate_and_encode;

/// An index for finding the highest-scoring reference sequence within a given
/// Hamming distance of a query sequence.
///
/// # Tie-breaking
///
/// When multiple references share the highest score, the one with the **lowest
/// original index** (i.e., its position in the input slice) is returned.
///
/// # Thread safety
///
/// `HammingResonatorWeighted` is `Send + Sync`; the index is read-only after construction.
#[derive(Debug)]
pub struct HammingResonatorWeighted {
    originals: Vec<BString>,
    scores: Vec<f64>,
    index: PartitionIndex,
}

impl HammingResonatorWeighted {
    /// Build with explicit `max_dist`.
    pub fn with_max_dist(
        seqs: Vec<(BString, f64)>,
        max_dist: u32,
    ) -> Result<Self, ResonateError> {
        let (bstrings, scores): (Vec<BString>, Vec<f64>) = seqs.into_iter().unzip();
        let encoded = validate_and_encode(&bstrings, max_dist)?;
        let index = PartitionIndex::build(encoded, max_dist)?;
        Ok(Self { originals: bstrings, scores, index })
    }

    /// Return the highest-scoring reference within `max_dist`, or `None` if no match is found.
    pub fn query_best(&self, query: &BStr) -> Result<Option<&BStr>, ResonateError> {
        let enc = encode(query)?;
        if enc.len() != self.index.seq_len {
            return Err(ResonateError::QueryLengthMismatch {
                got: enc.len(),
                expected: self.index.seq_len,
            });
        }

        let indices = self.index.query_indices(&enc);
        // query_indices returns sorted indices, so the first occurrence of the
        // maximum score wins the tie-break (lowest original index).
        let best = indices.into_iter().reduce(|best, idx| {
            if self.scores[idx as usize] > self.scores[best as usize] {
                idx
            } else {
                best
            }
        });

        Ok(best.map(|i| {
            let s: &BStr = self.originals[i as usize].as_ref();
            s
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::ByteSlice;

    fn bstr(s: &str) -> BString {
        BString::from(s)
    }

    fn weighted(seqs: &[(&str, f64)], max_dist: u32) -> HammingResonatorWeighted {
        let v: Vec<(BString, f64)> = seqs.iter().map(|&(s, w)| (bstr(s), w)).collect();
        HammingResonatorWeighted::with_max_dist(v, max_dist).unwrap()
    }

    #[test]
    fn returns_highest_score() {
        let r = weighted(&[("AAAA", 1.0), ("AAAC", 5.0), ("AACC", 10.0)], 1);
        // AACC is d=2 from AAAA, excluded; AAAC is d=1 and highest among candidates
        let hit = r.query_best("AAAA".as_bytes().as_bstr()).unwrap().unwrap();
        assert_eq!(hit, "AAAC".as_bytes().as_bstr());
    }

    #[test]
    fn tie_break_lowest_index() {
        let r = weighted(&[("AAAA", 5.0), ("AAAC", 5.0)], 1);
        let hit = r.query_best("AAAA".as_bytes().as_bstr()).unwrap().unwrap();
        // Both have score 5.0; index 0 wins
        assert_eq!(hit, "AAAA".as_bytes().as_bstr());
    }

    #[test]
    fn no_match_returns_none() {
        let r = weighted(&[("GGGG", 1.0)], 1);
        let hit = r.query_best("AAAA".as_bytes().as_bstr()).unwrap();
        assert!(hit.is_none());
    }

    #[test]
    fn invalid_base_in_query() {
        let r = weighted(&[("ACGT", 1.0)], 1);
        let err = r.query_best("ACGN".as_bytes().as_bstr()).unwrap_err();
        assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
    }
}
