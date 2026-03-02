use crate::error::ResonateError;
use bstr::BString;

/// Count the number of differing base positions between two equal-length byte slices.
#[inline]
pub(crate) fn hamming_distance(a: &[u8], b: &[u8]) -> u32 {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| (x != y) as u32).sum()
}

#[derive(Debug)]
pub struct EncodedSeqs {
    pub arena: Vec<u8>,
    pub entry_len: usize,  //how many bytes of sequence are there
    pub entry_size: usize, //how many bytes per entry (for weighted resonator)
    pub count: usize,
}

impl EncodedSeqs {
    /// Shared validation:
    ///     checks non-empty,
    ///     consistent length,
    ///     valid bases,
    ///     `seq_len` >= `max_dist` + 1.
    /// Returns encoded sequences ready for `PartitionIndex::build`.
    pub(crate) fn new(seqs: &[BString], max_dist: u32) -> Result<Self, ResonateError> {
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
            arena[start..end].copy_from_slice(seq);
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
    pub(crate) fn new_with_scores(
        seqs: &[BString],
        scores: &[f64],
        max_dist: u32,
    ) -> Result<Self, ResonateError> {
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
            arena[start..end].copy_from_slice(seq);
            arena[end..end + 8].copy_from_slice(&score.to_le_bytes());
        }

        Ok(EncodedSeqs {
            arena,
            count: seqs.len(),
            entry_len: expected_len,
            entry_size: bytes_per_entry,
        })
    }

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
