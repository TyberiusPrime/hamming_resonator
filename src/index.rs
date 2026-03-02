use bstr::BStr;

use crate::encode::{EncSeqs, hamming_distance};
use crate::error::ResonateError;
use std::collections::HashMap;

#[derive(Debug)]
pub(crate) struct PartitionIndex<T: EncSeqs> {
    /// Full sequences (and possibly scores) in form for Hamming distance computation.
    pub(crate) encoded: T,
    /// Keys stay byte-per-base so chunk boundaries map cleanly to slice ranges.
    pub(crate) chunk_maps: Vec<HashMap<Vec<u8>, Vec<u32>>>,
    /// Byte-per-base coordinate range (start, end) for each chunk.
    pub(crate) ranges: Vec<(usize, usize)>,
    pub(crate) seq_len: usize,
    pub(crate) max_dist: u32,
}

impl<T: EncSeqs> PartitionIndex<T> {
    /// Build from byte-per-base encoded sequences.
    /// That is ascii...
    pub(crate) fn build(encoded: T, max_dist: u32) -> Result<Self, ResonateError> {
        let seq_len = encoded.get_entry_len();
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
            //            encoded_nibble.push(encode_nibble(enc));
        }
        Ok(Self {
            encoded,
            chunk_maps,
            ranges,
            seq_len,
            max_dist,
        })
    }

    /// Return indices and distance of all references within `max_dist` of `enc` (byte-per-base).
    pub(crate) fn query_indices(&self, query: &BStr) -> Vec<(u32, u32, T::Score)> {
        let mut candidates: Vec<_> = Vec::new();
        let max_hamming_dist = self.max_dist;
        for (chunk_i, &(start, end)) in self.ranges.iter().enumerate() {
            if let Some(hits) = self.chunk_maps[chunk_i].get(&query[start..end]) {
                for idx in hits.into_iter() {
                    if !candidates.iter().any(|(cidx, _, _)| *cidx == *idx) {
                        let (slice, score) = self.encoded.get_entry(*idx);
                        let dist = hamming_distance(query, slice);
                        if dist <= max_hamming_dist {
                            candidates.push((*idx, dist, score));
                        }
                    }
                }
                //candidates.extend_from_slice(hits);
            }
        }
    
        candidates
    }
  

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::EncodedSeqs;
    use bstr::BString;

    fn build(seqs: &[&str], max_dist: u32) -> PartitionIndex<EncodedSeqs> {
        let bstr_seqs: Vec<BString> = seqs.iter().map(|&s| BString::from(s)).collect();
        let encoded = EncodedSeqs::new(&bstr_seqs, max_dist).unwrap();
        PartitionIndex::build(encoded, max_dist).unwrap()
    }

    #[test]
    fn exact_match_always_returned() {
        let idx = build(&["ACGT", "TTTT"], 1);
        let q = "ACGT".into();
        let hits = idx.query_indices(q);
        assert!(hits.contains(&(0, 0, ())));
    }

    #[test]
    fn d1_found_d2_not() {
        let idx = build(&["AAAA", "AAAC", "AACC"], 1);
        let q = "AAAA".into();
        let hits = idx.query_indices(q);
        assert!(hits.contains(&(0, 0, ())));
        assert!(hits.contains(&(1, 1, ())));
        assert!(!hits.contains(&(2, 2, ())));
    }

    #[test]
    fn no_matches_returns_empty() {
        let idx = build(&["GGGG", "CCCC"], 1);
        let q = "AAAA".into();
        assert!(idx.query_indices(q).is_empty());
    }

    #[test]
    fn dedup_candidates() {
        let idx = build(&["AAAA"], 1);
        let q = "AAAA".into();
        let hits = idx.query_indices(q);
        assert_eq!(hits, vec![(0, 0, ())]);
    }
}
