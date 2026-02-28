use crate::error::ResonateError;
use bstr::BStr;
use hamming_bitwise_fast::hamming_bitwise_fast;

/// Byte-per-base encoded DNA sequence. Values: A=0, C=1, G=2, T=4.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct EncodedSeq(pub(crate) Vec<u8>);

impl std::ops::Deref for EncodedSeq {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl<const N: usize> PartialEq<[u8; N]> for EncodedSeq {
    fn eq(&self, other: &[u8; N]) -> bool {
        self.0 == other.as_slice()
    }
}
/// Count the number of differing base positions between two equal-length encoded slices.
#[inline]
pub(crate) fn hamming(a: &[u8], b: &[u8]) -> u32 {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| (x != y) as u32).sum()
}
/// Nibble-packed DNA sequence used for fast Hamming distance comparison.
/// Two bases per byte: even index → high nibble, odd index → low nibble.
/// One-hot encoding: A=0x1, C=0x2, G=0x4, T=0x8.
#[derive(Debug, Clone)]
pub(crate) struct NibbleSeq(Vec<u8>);

impl std::ops::Deref for NibbleSeq {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl<const N: usize> PartialEq<[u8; N]> for NibbleSeq {
    fn eq(&self, other: &[u8; N]) -> bool {
        self.0.as_slice() == other.as_slice()
    }
}

/// Encode a DNA sequence (case-insensitive ACGT) into one byte per base.
/// A→0, C→1, G→2, T→4. Normalises to uppercase before mapping.
pub(crate) fn encode(seq: &BStr) -> Result<EncodedSeq, ResonateError> {
    let mut data = vec![0u8; seq.len()];
    for (i, &b) in seq.iter().enumerate() {
        data[i] = match b | 0x20 {
            b'a' => 0u8,
            b'c' => 1u8,
            b'g' => 2u8,
            b't' => 4u8,
            _ => return Err(ResonateError::InvalidBase(b as char, i)),
        };
    }
    Ok(EncodedSeq(data))
}

/// Pack a byte-per-base `EncodedSeq` into nibble form: 2 bases per byte.
/// Even index → high nibble (bits 7:4), odd index → low nibble (bits 3:0).
/// Any two distinct bases XOR to exactly 2 set bits, so
/// `hamming_bitwise_fast(a, b) / 2` gives the correct base-level distance.
/// Odd-length sequences pad the final low nibble with 0x0; both sides carry the
/// same padding, so it contributes 0 to the XOR and does not affect the result.
pub(crate) fn encode_nibble(byte_enc: &EncodedSeq) -> NibbleSeq {
    let packed_len = byte_enc.len().div_ceil(2);
    let mut out = vec![0u8; packed_len];
    for (i, &b) in byte_enc.iter().enumerate() {
        let nibble: u8 = match b {
            0 => 0x1, // A
            1 => 0x2, // C
            2 => 0x4, // G
            4 => 0x8, // T
            _ => unreachable!(),
        };
        if i % 2 == 0 {
            out[i / 2] |= nibble << 4;
        } else {
            out[i / 2] |= nibble;
        }
    }
    NibbleSeq(out)
}

/// Count differing base positions between two equal-length nibble-packed sequences.
#[inline]
pub(crate) fn hamming_nibble(a: &NibbleSeq, b: &NibbleSeq) -> u32 {
    hamming_bitwise_fast(a, b) / 2
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::ByteSlice;

    #[test]
    fn encode_uppercase() {
        let enc = encode("ACGT".as_bytes().as_bstr()).unwrap();
        assert_eq!(enc, [0u8, 1, 2, 4]);
    }

    #[test]
    fn encode_lowercase() {
        let enc = encode("acgt".as_bytes().as_bstr()).unwrap();
        assert_eq!(enc, [0u8, 1, 2, 4]);
    }

    #[test]
    fn encode_invalid_base() {
        let err = encode("ACGN".as_bytes().as_bstr()).unwrap_err();
        assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
    }

    // Nibble encoding: ACGT → 0x12, 0x48 (even→high nibble, odd→low nibble)
    #[test]
    fn nibble_encode_acgt() {
        let byte_enc = encode("ACGT".as_bytes().as_bstr()).unwrap(); // [0,1,2,4]
        let nibble = encode_nibble(&byte_enc);
        // A(0→0x1) in high nibble, C(1→0x2) in low nibble → 0x12
        // G(2→0x4) in high nibble, T(4→0x8) in low nibble → 0x48
        assert_eq!(nibble, [0x12, 0x48]);
    }

    #[test]
    fn nibble_encode_odd_length() {
        let byte_enc = encode("ACG".as_bytes().as_bstr()).unwrap(); // [0,1,2]
        let nibble = encode_nibble(&byte_enc);
        // A→high, C→low → 0x12; G→high, pad 0→low → 0x40
        assert_eq!(nibble, [0x12, 0x40]);
    }

    #[test]
    fn hamming_nibble_identical() {
        let a = encode("ACGT".as_bytes().as_bstr()).unwrap();
        let na = encode_nibble(&a);
        assert_eq!(hamming_nibble(&na, &na), 0);
    }

    #[test]
    fn hamming_nibble_one_diff() {
        let a = encode("AAAA".as_bytes().as_bstr()).unwrap();
        let b = encode("AAAC".as_bytes().as_bstr()).unwrap();
        assert_eq!(hamming_nibble(&encode_nibble(&a), &encode_nibble(&b)), 1);
    }

    #[test]
    fn hamming_nibble_all_pairs() {
        // Verify every distinct-base pair contributes exactly 1 to the distance.
        let bases = [b'A', b'C', b'G', b'T'];
        for &x in &bases {
            for &y in &bases {
                let a = encode_nibble(&encode(&[x].as_bstr()).unwrap());
                let b = encode_nibble(&encode(&[y].as_bstr()).unwrap());
                let expected = if x == y { 0 } else { 1 };
                assert_eq!(hamming_nibble(&a, &b), expected, "{x} vs {y}");
            }
        }
    }

    #[test]
    fn hamming_nibble_all_diff() {
        let a = encode("ACGT".as_bytes().as_bstr()).unwrap();
        let b = encode("TGCA".as_bytes().as_bstr()).unwrap();
        assert_eq!(hamming_nibble(&encode_nibble(&a), &encode_nibble(&b)), 4);
    }

    #[test]
    fn hamming_nibble_odd_length() {
        // Odd-length: padding nibble must not contribute to the distance.
        let a = encode("ACG".as_bytes().as_bstr()).unwrap();
        let b = encode("ACT".as_bytes().as_bstr()).unwrap(); // 1 diff at pos 2
        assert_eq!(hamming_nibble(&encode_nibble(&a), &encode_nibble(&b)), 1);
    }
}
