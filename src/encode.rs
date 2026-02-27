use bstr::BStr;
use crate::error::ResonateError;

/// Encode a DNA sequence (case-insensitive ACGT) into one byte per base.
/// A→0, C→1, G→2, T→4. Normalises to uppercase before mapping.
pub(crate) fn encode(seq: &BStr) -> Result<Vec<u8>, ResonateError> {
    seq.iter()
        .enumerate()
        .map(|(i, &b)| match b | 0x20 {
            b'a' => Ok(0u8),
            b'c' => Ok(1u8),
            b'g' => Ok(2u8),
            b't' => Ok(4u8),
            _ => Err(ResonateError::InvalidBase(b as char, i)),
        })
        .collect()
}

/// Count the number of differing base positions between two equal-length encoded slices.
#[inline]
pub(crate) fn hamming(a: &[u8], b: &[u8]) -> u32 {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| (x != y) as u32).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::ByteSlice;

    #[test]
    fn encode_uppercase() {
        let enc = encode("ACGT".as_bytes().as_bstr()).unwrap();
        assert_eq!(enc, [0, 1, 2, 4]);
    }

    #[test]
    fn encode_lowercase() {
        let enc = encode("acgt".as_bytes().as_bstr()).unwrap();
        assert_eq!(enc, [0, 1, 2, 4]);
    }

    #[test]
    fn encode_invalid_base() {
        let err = encode("ACGN".as_bytes().as_bstr()).unwrap_err();
        assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
    }

    #[test]
    fn hamming_identical() {
        assert_eq!(hamming(&[0, 1, 2, 4], &[0, 1, 2, 4]), 0);
    }

    #[test]
    fn hamming_one_diff() {
        assert_eq!(hamming(&[0, 1, 2, 4], &[1, 1, 2, 4]), 1);
    }

    #[test]
    fn hamming_all_diff() {
        assert_eq!(hamming(&[0, 1, 2, 4], &[4, 2, 1, 0]), 4);
    }
}
