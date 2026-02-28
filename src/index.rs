use crate::encode::{encode_nibble, hamming, hamming_nibble, EncodedSeq, NibbleSeq};
use crate::error::ResonateError;
use crate::resonator::EncodedSeqs;
use std::collections::HashMap;

#[derive(Debug)]
pub(crate) struct PartitionIndex {
    /// Full sequences in nibble-packed form for Hamming distance computation.
    encoded: EncodedSeqs,
    /// `chunk_maps`[i]: byte-per-base chunk slice → list of ref indices.
    /// Keys stay byte-per-base so chunk boundaries map cleanly to slice ranges.
    chunk_maps: Vec<HashMap<Vec<u8>, Vec<u32>>>,
    /// Byte-per-base coordinate range (start, end) for each chunk.
    ranges: Vec<(usize, usize)>,
    pub(crate) seq_len: usize,
    pub(crate) max_dist: u32,
}

impl PartitionIndex {
    /// Build from byte-per-base encoded sequences.
    pub(crate) fn build(encoded: EncodedSeqs, max_dist: u32) -> Result<Self, ResonateError> {
        let seq_len = encoded.entry_len;
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
    pub(crate) fn query_indices(&self, enc: &EncodedSeq) -> Vec<(u32, u32)> {
        // Nibble-pack the query once; reused for every candidate distance check.

        let mut candidates: Vec<u32> = Vec::new();
        for (chunk_i, &(start, end)) in self.ranges.iter().enumerate() {
            if let Some(hits) = self.chunk_maps[chunk_i].get(&enc[start..end]) {
                candidates.extend_from_slice(hits);
            }
        }

        candidates.sort_unstable();
        candidates.dedup();

        let mut candidates: Vec<_> = candidates
            .into_iter()
            .map(|idx| {
                let slice = self.encoded.get_entry(idx);
                (idx, hamming(enc, slice))
            })
            .collect();
        candidates.retain(|(idx, dist)| *dist <= self.max_dist);
        //hamming_nibble(&enc_nibble, &self.encoded_nibble[idx as usize]) <= self.max_dist

        candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{encode::encode, resonator::validate_and_encode};
    use bstr::{BString, ByteSlice};

    fn enc(s: &str) -> EncodedSeq {
        encode(s.as_bytes().as_bstr()).unwrap()
    }

    fn build(seqs: &[&str], max_dist: u32) -> PartitionIndex {
        let bstrSeqs: Vec<BString> = seqs.iter().map(|&s| BString::from(s)).collect();
        let encoded: EncodedSeqs = validate_and_encode(&bstrSeqs, max_dist).unwrap();
        PartitionIndex::build(encoded, max_dist).unwrap()
    }

    #[test]
    fn exact_match_always_returned() {
        let idx = build(&["ACGT", "TTTT"], 1);
        let q = enc("ACGT");
        let hits = idx.query_indices(&q);
        assert!(hits.contains(&(0,0)));
    }

    #[test]
    fn d1_found_d2_not() {
        let idx = build(&["AAAA", "AAAC", "AACC"], 1);
        let q = enc("AAAA");
        let hits = idx.query_indices(&q);
        assert!(hits.contains(&(0,0)));
        assert!(hits.contains(&(1,1)));
        assert!(!hits.contains(&(2,2)));
    }

    #[test]
    fn no_matches_returns_empty() {
        let idx = build(&["GGGG", "CCCC"], 1);
        let q = enc("AAAA");
        assert!(idx.query_indices(&q).is_empty());
    }

    #[test]
    fn dedup_candidates() {
        let idx = build(&["AAAA"], 1);
        let q = enc("AAAA");
        let hits = idx.query_indices(&q);
        assert_eq!(hits, vec![(0,0)]);
    }
}
