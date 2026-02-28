use bstr::{BStr, BString, ByteSlice};
use hamming_resonate::{HammingResonator, //HammingResonatorWeighted,
ResonateError};

fn b(s: &str) -> BString {
    BString::from(s)
}

fn resonator(seqs: &[&str], max_dist: u32) -> HammingResonator {
    let v: Vec<BString> = seqs.iter().map(|&s| b(s)).collect();
    HammingResonator::with_max_dist(v, max_dist).unwrap()
}

// fn weighted(seqs: &[(&str, f64)], max_dist: u32) -> HammingResonatorWeighted {
//     let v: Vec<(BString, f64)> = seqs.iter().map(|&(s, w)| (b(s), w)).collect();
//     HammingResonatorWeighted::with_max_dist(v, max_dist).unwrap()
// }

fn hit_strings(hits: Vec<&BStr>) -> Vec<String> {
    let mut v: Vec<String> = hits.iter().map(|s| s.to_str().unwrap().to_owned()).collect();
    v.sort();
    v
}

// Test 1: Exact match (distance 0) is always returned
#[test]
fn t01_exact_match_always_returned() {
    let r = resonator(&["ACGT", "TTTT", "GGGG"], 2);
    let hits = r.query("ACGT".as_bytes().as_bstr()).unwrap();
    assert!(hits.contains(&"ACGT".as_bytes().as_bstr()));
}

// Test 2: d=1 match found; d=2 sequence not returned at max_dist=1
#[test]
fn t02_d1_found_d2_excluded() {
    let r = resonator(&["AAAA", "AAAC", "AACC"], 1);
    let hits = hit_strings(r.query("AAAA".as_bytes().as_bstr()).unwrap());
    assert!(hits.contains(&"AAAA".to_owned()));
    assert!(hits.contains(&"AAAC".to_owned()));
    assert!(!hits.contains(&"AACC".to_owned()));
}

// Test 3: d=2 and d=3 matches found at correct max_dist
#[test]
fn t03_d2_and_d3_at_correct_max_dist() {
    let r = resonator(&["AAAA", "AACC", "ACCC"], 3);
    let hits = hit_strings(r.query("AAAA".as_bytes().as_bstr()).unwrap());
    assert!(hits.contains(&"AAAA".to_owned())); // d=0
    assert!(hits.contains(&"AACC".to_owned())); // d=2
    assert!(hits.contains(&"ACCC".to_owned())); // d=3
}

// Test 4: No matches → Ok(vec![]) / Ok(None)
// #[test]
// fn t04_no_matches() {
//     let r = resonator(&["GGGG", "CCCC"], 1);
//     assert!(r.query("AAAA".as_bytes().as_bstr()).unwrap().is_empty());
//
//     let w = weighted(&[("GGGG", 1.0), ("CCCC", 1.0)], 1);
//     assert!(w.query_best("AAAA".as_bytes().as_bstr()).unwrap().is_none());
// }

// Test 5: Invalid base in reference → Err(InvalidBase)
#[test]
fn t05_invalid_base_in_reference() {
    let seqs = vec![b("ACGN")];
    let err = HammingResonator::with_max_dist(seqs, 1).unwrap_err();
    assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
}

// Test 6: Invalid base in query → Err(InvalidBase)
#[test]
fn t06_invalid_base_in_query() {
    let r = resonator(&["ACGT"], 1);
    let err = r.query("ACGN".as_bytes().as_bstr()).unwrap_err();
    assert!(matches!(err, ResonateError::InvalidBase('N', 3)));
}

// Test 7: Mismatched query length → Err(QueryLengthMismatch)
#[test]
fn t07_query_length_mismatch() {
    let r = resonator(&["ACGT"], 1);
    let err = r.query("ACG".as_bytes().as_bstr()).unwrap_err();
    assert!(matches!(
        err,
        ResonateError::QueryLengthMismatch { got: 3, expected: 4 }
    ));
}

// Test 8: query_best returns highest-score match, not just first or closest
// #[test]
// fn t08_query_best_highest_score() {
//     // AAAC is d=1 (score 5.0), AAAA is d=0 (score 1.0); best should be AAAC
//     let w = weighted(&[("AAAA", 1.0), ("AAAC", 5.0)], 1);
//     let hit = w.query_best("AAAA".as_bytes().as_bstr()).unwrap().unwrap();
//     assert_eq!(hit, "AAAC".as_bytes().as_bstr());
// }

// Test 9: query_best tie-break: lowest index wins
// #[test]
// fn t09_query_best_tiebreak_lowest_index() {
//     let w = weighted(&[("AAAA", 5.0), ("AAAC", 5.0)], 1);
//     let hit = w.query_best("AAAA".as_bytes().as_bstr()).unwrap().unwrap();
//     assert_eq!(hit, "AAAA".as_bytes().as_bstr());
// }

// Test 10: All references identical — single result returned
#[test]
fn t10_identical_references() {
    let r = resonator(&["AAAA", "AAAA", "AAAA"], 1);
    let hits = r.query("AAAA".as_bytes().as_bstr()).unwrap();
    // All three are valid matches, indices [0,1,2] all qualify
    // The key requirement: query returns references for all matching indices
    assert!(!hits.is_empty());
    // All results should be "AAAA"
    assert!(hits.iter().all(|h| *h == "AAAA".as_bytes().as_bstr()));
}

// Test 11: Case-insensitive bases (a/c/g/t == A/C/G/T)
#[test]
fn t11_case_insensitive() {
    let r = resonator(&["acgt"], 1);
    let hits = r.query("ACGT".as_bytes().as_bstr()).unwrap();
    assert_eq!(hits.len(), 1);

    let r2 = resonator(&["ACGT"], 1);
    let hits2 = r2.query("acgt".as_bytes().as_bstr()).unwrap();
    assert_eq!(hits2.len(), 1);
}

// Test 12: query_batch results match serial query results
#[test]
fn t12_query_batch_matches_serial() {
    let r = resonator(&["AAAA", "AAAC", "AAGC", "TTTT"], 1);
    let queries: Vec<BString> = vec![b("AAAA"), b("AAAC"), b("TTTT"), b("CCCC")];

    let batch = r.query_batch(&queries);
    let serial: Vec<_> = queries.iter().map(|q| r.query(q.as_ref())).collect();

    for (b_res, s_res) in batch.iter().zip(serial.iter()) {
        let mut bv = hit_strings(b_res.as_ref().unwrap().clone());
        let mut sv = hit_strings(s_res.as_ref().unwrap().clone());
        bv.sort();
        sv.sort();
        assert_eq!(bv, sv);
    }
}

// Validation errors
#[test]
fn empty_input_error() {
    let err = HammingResonator::with_max_dist(vec![], 1).unwrap_err();
    assert!(matches!(err, ResonateError::EmptyInput));
}

#[test]
fn inconsistent_length_error() {
    let seqs = vec![b("ACGT"), b("ACG")];
    let err = HammingResonator::with_max_dist(seqs, 1).unwrap_err();
    assert!(matches!(err, ResonateError::InconsistentLength(4, 3)));
}

#[test]
fn sequence_too_short_error() {
    // seq_len=2, max_dist=2 requires at least 3 bases
    let seqs = vec![b("AC")];
    let err = HammingResonator::with_max_dist(seqs, 2).unwrap_err();
    assert!(matches!(
        err,
        ResonateError::SequenceTooShort { seq_len: 2, max_dist: 2, min: 3 }
    ));
}
