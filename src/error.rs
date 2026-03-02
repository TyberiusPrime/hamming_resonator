#[derive(thiserror::Error, Debug, PartialEq)]
pub enum ResonateError {
    #[error("Empty reference set")]
    EmptyInput,

    #[error("All sequences must have the same length, got {0} and {1}")]
    InconsistentLength(usize, usize),

    #[error("Sequence length {seq_len} too short for max_dist={max_dist} (need at least {min} bases)")]
    SequenceTooShort { seq_len: usize, max_dist: u32, min: usize },

    #[error("Query length {got} does not match index length {expected}")]
    QueryLengthMismatch { got: usize, expected: usize },

    #[error("Duplicate entries in reference set")]
    DuplicateEntries,

    #[error("Seqs and scores had differing lengths")] // todo refactor out to be impossible
    LengthMismatch ,
}
