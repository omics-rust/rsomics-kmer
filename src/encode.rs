use crate::{KmerError, Result};

/// A DNA k-mer encoded as two bits per base in the low `2 * k` bits.
///
/// Values passed to codec functions must be normalized: every bit above the
/// low `2 * k` bits is zero. All constructors in this crate produce normalized
/// values.
pub type Kmer = u64;

pub(crate) const MAX_K: usize = 32;

pub(crate) fn validate_k(k: usize) -> Result<()> {
    if !(1..=MAX_K).contains(&k) {
        return Err(KmerError::KOutOfRange(k));
    }
    Ok(())
}

const fn mask_for_valid_k(k: usize) -> Kmer {
    if k == MAX_K {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

fn validate_kmer(kmer: Kmer, k: usize) -> Result<Kmer> {
    validate_k(k)?;
    let mask = mask_for_valid_k(k);
    if kmer & !mask != 0 {
        return Err(KmerError::NonNormalizedKmer { kmer, k });
    }
    Ok(mask)
}

/// Returns the two-bit encoding for an ASCII DNA base.
///
/// Upper- and lowercase A/C/G/T are accepted. Every other byte returns
/// `None`.
pub const fn base_bits(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

const fn bits_base(bits: u64) -> u8 {
    match bits & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        _ => b'T',
    }
}

/// Encodes one non-empty A/C/G/T sequence into a normalized [`Kmer`].
///
/// # Errors
///
/// Returns [`KmerError::KOutOfRange`] unless the sequence length is in
/// `1..=32`, or [`KmerError::NonAcgt`] for the first non-ACGT byte.
pub fn encode(seq: &[u8]) -> Result<Kmer> {
    let k = seq.len();
    validate_k(k)?;
    let mut bits: u64 = 0;
    for (i, &b) in seq.iter().enumerate() {
        let v = base_bits(b).ok_or(KmerError::NonAcgt { pos: i, byte: b })?;
        bits = (bits << 2) | v;
    }
    Ok(bits)
}

/// Decodes a normalized two-bit k-mer into uppercase ASCII A/C/G/T.
///
/// # Panics
///
/// Panics when `k` is outside `1..=32` or `kmer` has bits set above its low
/// `2 * k`-bit representation. Use [`try_decode`] for runtime-validated input.
#[must_use]
pub fn decode(kmer: Kmer, k: usize) -> Vec<u8> {
    try_decode(kmer, k).expect("decode requires a normalized k-mer with k in 1..=32")
}

/// Decodes a normalized two-bit k-mer into uppercase ASCII A/C/G/T.
///
/// # Errors
///
/// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`, or
/// [`KmerError::NonNormalizedKmer`] when bits above the representation are set.
pub fn try_decode(kmer: Kmer, k: usize) -> Result<Vec<u8>> {
    validate_kmer(kmer, k)?;
    let mut out = vec![0u8; k];
    let mut bits = kmer;
    for slot in out.iter_mut().rev() {
        *slot = bits_base(bits);
        bits >>= 2;
    }
    Ok(out)
}

/// Returns the reverse complement of a normalized two-bit k-mer.
///
/// # Panics
///
/// Panics when `k` is outside `1..=32` or `kmer` is not normalized. Use
/// [`try_reverse_complement`] for runtime-validated input.
#[must_use]
pub fn reverse_complement(kmer: Kmer, k: usize) -> Kmer {
    try_reverse_complement(kmer, k)
        .expect("reverse_complement requires a normalized k-mer with k in 1..=32")
}

/// Returns the reverse complement of a normalized two-bit k-mer.
///
/// # Errors
///
/// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`, or
/// [`KmerError::NonNormalizedKmer`] when bits above the representation are set.
pub fn try_reverse_complement(kmer: Kmer, k: usize) -> Result<Kmer> {
    let mask = validate_kmer(kmer, k)?;
    let mut bits = kmer ^ mask;
    let mut rc: u64 = 0;
    for _ in 0..k {
        rc = (rc << 2) | (bits & 0b11);
        bits >>= 2;
    }
    Ok(rc)
}

/// Selects the numerically smaller of a normalized k-mer and its reverse complement.
///
/// For this two-bit representation, numeric ordering at a fixed `k` is
/// lexicographic ordering under `A < C < G < T`.
///
/// # Panics
///
/// Panics when `k` is outside `1..=32` or `kmer` is not normalized. Use
/// [`try_canonical`] for runtime-validated input.
#[must_use]
pub fn canonical(kmer: Kmer, k: usize) -> Kmer {
    try_canonical(kmer, k).expect("canonical requires a normalized k-mer with k in 1..=32")
}

/// Selects the numerically smaller of a normalized k-mer and its reverse complement.
///
/// # Errors
///
/// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`, or
/// [`KmerError::NonNormalizedKmer`] when bits above the representation are set.
pub fn try_canonical(kmer: Kmer, k: usize) -> Result<Kmer> {
    let rc = try_reverse_complement(kmer, k)?;
    Ok(kmer.min(rc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_round_trip() {
        let k = 8;
        let seq = b"ACGTACGT";
        let bits = encode(seq).unwrap();
        let back = decode(bits, k);
        assert_eq!(&back, seq);
    }

    #[test]
    fn encode_rejects_n() {
        let r = encode(b"ACGTN");
        assert!(matches!(r, Err(KmerError::NonAcgt { pos: 4, .. })));
    }

    #[test]
    fn encode_k_too_large_rejected() {
        let seq = vec![b'A'; 33];
        assert!(matches!(encode(&seq), Err(KmerError::KOutOfRange(33))));
    }

    #[test]
    fn rc_of_rc_is_identity() {
        let bits = encode(b"ACGTACGT").unwrap();
        let rc = reverse_complement(bits, 8);
        let back = reverse_complement(rc, 8);
        assert_eq!(back, bits);
    }

    #[test]
    fn rc_known_value() {
        assert_eq!(
            decode(reverse_complement(encode(b"AAAA").unwrap(), 4), 4),
            b"TTTT".to_vec()
        );
        assert_eq!(
            decode(reverse_complement(encode(b"ACGT").unwrap(), 4), 4),
            b"ACGT".to_vec()
        );
    }

    #[test]
    fn rc_k32_round_trip() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGA";
        let bits = encode(seq).unwrap();
        let rc = reverse_complement(bits, 32);
        assert_eq!(decode(rc, 32), b"TCGTACGTACGTACGTACGTACGTACGTACGT".to_vec());
        assert_eq!(reverse_complement(rc, 32), bits);
    }

    #[test]
    fn canonical_picks_lex_min_of_pair() {
        let fwd = encode(b"GGGG").unwrap();
        let rc = reverse_complement(fwd, 4);
        assert_eq!(canonical(fwd, 4), rc);
        assert_eq!(canonical(fwd, 4), canonical(rc, 4));
    }

    #[test]
    fn checked_codecs_cover_k_boundaries() {
        for k in [0, 33] {
            assert!(matches!(
                try_decode(0, k),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
            assert!(matches!(
                try_reverse_complement(0, k),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
            assert!(matches!(
                try_canonical(0, k),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
        }

        for k in [1, 31, 32] {
            let seq: Vec<u8> = (0..k).map(|i| b"ACGT"[(i * 3 + k) % 4]).collect();
            let encoded = encode(&seq).unwrap();
            let rc = try_reverse_complement(encoded, k).unwrap();
            assert_eq!(try_decode(encoded, k).unwrap(), seq);
            assert_eq!(try_reverse_complement(rc, k).unwrap(), encoded);
            assert_eq!(try_canonical(encoded, k).unwrap(), encoded.min(rc));
            assert_eq!(
                try_canonical(encoded, k).unwrap(),
                try_canonical(rc, k).unwrap()
            );
        }
    }

    #[test]
    fn non_normalized_values_are_rejected_consistently() {
        for k in [1, 31] {
            let non_normalized = 1u64 << (2 * k);
            assert!(matches!(
                try_decode(non_normalized, k),
                Err(KmerError::NonNormalizedKmer {
                    kmer,
                    k: actual_k
                }) if kmer == non_normalized && actual_k == k
            ));
            assert!(matches!(
                try_reverse_complement(non_normalized, k),
                Err(KmerError::NonNormalizedKmer { .. })
            ));
            assert!(matches!(
                try_canonical(non_normalized, k),
                Err(KmerError::NonNormalizedKmer { .. })
            ));
            assert!(std::panic::catch_unwind(|| decode(non_normalized, k)).is_err());
            assert!(std::panic::catch_unwind(|| reverse_complement(non_normalized, k)).is_err());
            assert!(std::panic::catch_unwind(|| canonical(non_normalized, k)).is_err());
        }
    }

    #[test]
    fn every_u64_is_normalized_at_k32() {
        for encoded in [0, 1, u64::MAX / 3, u64::MAX] {
            let decoded = try_decode(encoded, 32).unwrap();
            assert_eq!(decoded.len(), 32);
            assert_eq!(
                try_reverse_complement(try_reverse_complement(encoded, 32).unwrap(), 32).unwrap(),
                encoded
            );
        }
    }

    #[test]
    fn fail_loud_wrappers_reject_out_of_range_k() {
        for k in [0, 33] {
            assert!(std::panic::catch_unwind(|| decode(0, k)).is_err());
            assert!(std::panic::catch_unwind(|| reverse_complement(0, k)).is_err());
            assert!(std::panic::catch_unwind(|| canonical(0, k)).is_err());
        }
    }
}
