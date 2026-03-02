use bstr::{BStr, BString};

use crate::encode::{EncodedSeqsAndScores};
use crate::error::ResonateError;
use crate::index::PartitionIndex;

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
    index: PartitionIndex<EncodedSeqsAndScores>,
}

impl HammingResonatorWeighted {
    /// Build with explicit `max_dist`.
    pub fn with_max_dist(seqs: Vec<(BString, f32)>, max_dist: u32) -> Result<Self, ResonateError> {
        let (bstrings, scores): (Vec<BString>, Vec<f32>) = seqs.into_iter().unzip();
        let encoded = EncodedSeqsAndScores::new(&bstrings, &scores, max_dist)?;
        let index = PartitionIndex::build(encoded, max_dist)?;
        Ok(Self {
            originals: bstrings,
            index,
        })
    }

    /// Return the highest-scoring reference within `max_dist`, or `None` if no match is found.
    pub fn query_best(&self, query: &BStr) -> Result<Option<&BStr>, ResonateError> {
        if query.len() != self.index.seq_len {
            return Err(ResonateError::QueryLengthMismatch {
                got: query.len(),
                expected: self.index.seq_len,
            });
        }

        let indices = self.index.query_indices(query);
        // query_indices returns sorted indices, so the first occurrence of the
        // maximum score wins the tie-break (lowest original index).
        // followed by lowest index
        let best = indices.into_iter().reduce(|best, (idx, d, score)| {
            if d < best.1
                || (d == best.1 && score > best.2)
            {
                (idx, d, score)
            } else {
                best
            }
        });

        Ok(best.map(|(i, _d, _score)| {
            let s: &BStr = self.originals[i as usize].as_ref();
            s
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::ByteSlice;

    

    fn weighted(seqs: &[(&str, f32)], max_dist: u32) -> HammingResonatorWeighted {
        let v: Vec<(BString, f32)> = seqs.iter().map(|&(s, w)| (BString::from(s), w)).collect();
        HammingResonatorWeighted::with_max_dist(v, max_dist).unwrap()
    }

    #[test]
    fn returns_highest_score() {
        let r = weighted(&[("AAAA", 1.0), ("AAAC", 5.0), ("AACC", 10.0)], 1);
        // AACC is d=2 from AAAA, excluded; AAAC is d=1 and highest among candidates
        let hit = r.query_best("AAAA".as_bytes().as_bstr()).unwrap().unwrap();
        assert_eq!(hit, "AAAA".as_bytes().as_bstr());
    }

    #[test]
    fn tie_break_lowest_index() {
        let r = weighted(&[("AAAA", 5.0), ("AAAC", 5.0)], 1);
        let hit = r
            .query_best("AAAA".into())
            .expect("query failed")
            .expect("No result");
        // Both have score 5.0; index 0 wins
        assert_eq!(hit, "AAAA".as_bytes().as_bstr());
    }

    #[test]
    fn no_match_returns_none() {
        let r = weighted(&[("GGGG", 1.0)], 1);
        let hit = r.query_best("AAAA".as_bytes().as_bstr()).unwrap();
        assert!(hit.is_none());
    }
}
