use std::collections::HashMap;
use crate::encode::hamming;
use crate::error::ResonateError;

#[derive(Debug)]
pub(crate) struct PartitionIndex {
    /// Encoded reference sequences (one byte per base).
    encoded: Vec<Vec<u8>>,
    /// chunk_maps[i]: chunk_bytes → list of ref indices
    chunk_maps: Vec<HashMap<Vec<u8>, Vec<u32>>>,
    /// Byte range (start, end) for each chunk position.
    ranges: Vec<(usize, usize)>,
    pub(crate) seq_len: usize,
    pub(crate) max_dist: u32,
}

impl PartitionIndex {
    pub(crate) fn build(
        encoded: Vec<Vec<u8>>,
        max_dist: u32,
    ) -> Result<Self, ResonateError> {
        let seq_len = encoded[0].len();
        let n_chunks = (max_dist + 1) as usize;

        let ranges: Vec<(usize, usize)> = (0..n_chunks)
            .map(|i| {
                let start = i * seq_len / n_chunks;
                let end = (i + 1) * seq_len / n_chunks;
                (start, end)
            })
            .collect();

        let mut chunk_maps: Vec<HashMap<Vec<u8>, Vec<u32>>> =
            (0..n_chunks).map(|_| HashMap::new()).collect();

        for (idx, enc) in encoded.iter().enumerate() {
            for (chunk_i, &(start, end)) in ranges.iter().enumerate() {
                chunk_maps[chunk_i]
                    .entry(enc[start..end].to_vec())
                    .or_default()
                    .push(idx as u32);
            }
        }

        Ok(Self { encoded, chunk_maps, ranges, seq_len, max_dist })
    }

    /// Return indices of all references within `max_dist` of `enc`.
    pub(crate) fn query_indices(&self, enc: &[u8]) -> Vec<u32> {
        let mut candidates: Vec<u32> = Vec::new();

        for (chunk_i, &(start, end)) in self.ranges.iter().enumerate() {
            if let Some(hits) = self.chunk_maps[chunk_i].get(&enc[start..end]) {
                candidates.extend_from_slice(hits);
            }
        }

        candidates.sort_unstable();
        candidates.dedup();

        candidates.retain(|&idx| {
            hamming(enc, &self.encoded[idx as usize]) <= self.max_dist
        });

        candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::encode;
    use bstr::ByteSlice;

    fn enc(s: &str) -> Vec<u8> {
        encode(s.as_bytes().as_bstr()).unwrap()
    }

    fn build(seqs: &[&str], max_dist: u32) -> PartitionIndex {
        let encoded: Vec<Vec<u8>> = seqs.iter().map(|s| enc(s)).collect();
        PartitionIndex::build(encoded, max_dist).unwrap()
    }

    #[test]
    fn exact_match_always_returned() {
        let idx = build(&["ACGT", "TTTT"], 1);
        let q = enc("ACGT");
        let hits = idx.query_indices(&q);
        assert!(hits.contains(&0));
    }

    #[test]
    fn d1_found_d2_not() {
        let idx = build(&["AAAA", "AAAC", "AACC"], 1);
        let q = enc("AAAA");
        let hits = idx.query_indices(&q);
        assert!(hits.contains(&0)); // d=0
        assert!(hits.contains(&1)); // d=1
        assert!(!hits.contains(&2)); // d=2, excluded
    }

    #[test]
    fn no_matches_returns_empty() {
        let idx = build(&["GGGG", "CCCC"], 1);
        let q = enc("AAAA");
        assert!(idx.query_indices(&q).is_empty());
    }

    #[test]
    fn dedup_candidates() {
        // When multiple chunks point to the same ref, it should appear once.
        let idx = build(&["AAAA"], 1);
        let q = enc("AAAA");
        let hits = idx.query_indices(&q);
        assert_eq!(hits, vec![0]);
    }
}
