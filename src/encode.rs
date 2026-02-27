use bstr::BStr;
use hamming_bitwise_fast::hamming_bitwise_fast;
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

/// Pack a byte-per-base slice into nibble form: 2 bases per byte.
/// Even index → high nibble (bits 7:4), odd index → low nibble (bits 3:0).
/// One-hot nibble values: A=0x1, C=0x2, G=0x4, T=0x8.
/// Any two distinct bases XOR to exactly 2 set bits, so
/// `hamming_bitwise_fast(a, b) / 2` gives the correct base-level distance.
/// Odd-length sequences pad the final low nibble with 0x0; both sides carry the
/// same padding, so it contributes 0 to the XOR and does not affect the result.
pub(crate) fn encode_nibble(byte_enc: &[u8]) -> Vec<u8> {
    let packed_len = (byte_enc.len() + 1) / 2;
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
    out
}

/// Count differing base positions between two equal-length nibble-packed slices.
#[inline]
pub(crate) fn hamming_nibble(a: &[u8], b: &[u8]) -> u32 {
    hamming_bitwise_fast(a, b) / 2
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
