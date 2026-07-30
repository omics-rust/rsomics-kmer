use crate::encode::{Kmer, MAX_K, base_bits};
use crate::{KmerError, Result};

/// A linear-time scanner that reports the k-mer ending at each input position.
///
/// The iterator yields one item per input base. The first `k - 1` items are
/// `None`; a non-ACGT byte resets the rolling state, so every window containing
/// that byte is also reported as `None`.
pub struct RollingKmers<'a> {
    seq: &'a [u8],
    k: usize,
    mask: u64,
    pos: usize,
    current: u64,
    valid: usize,
}

impl<'a> RollingKmers<'a> {
    /// Creates a rolling scanner.
    ///
    /// # Panics
    ///
    /// Panics when `k` is outside `1..=32`. Use [`Self::try_new`] for a
    /// runtime-validated length.
    #[must_use]
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        Self::try_new(seq, k).expect("RollingKmers::new requires k in 1..=32")
    }

    /// Creates a rolling scanner after validating `k`.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`.
    pub fn try_new(seq: &'a [u8], k: usize) -> Result<Self> {
        if !(1..=MAX_K).contains(&k) {
            return Err(KmerError::KOutOfRange(k));
        }
        Ok(Self {
            seq,
            k,
            mask: if k == MAX_K {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            },
            pos: 0,
            current: 0,
            valid: 0,
        })
    }
}

impl Iterator for RollingKmers<'_> {
    type Item = Option<Kmer>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.seq.len() {
            return None;
        }
        let b = self.seq[self.pos];
        self.pos += 1;

        if let Some(bits) = base_bits(b) {
            self.current = ((self.current << 2) | bits) & self.mask;
            self.valid += 1;
        } else {
            self.current = 0;
            self.valid = 0;
        }

        if self.valid >= self.k {
            Some(Some(self.current))
        } else {
            Some(None)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.seq.len() - self.pos;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for RollingKmers<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::encode;

    #[test]
    fn rolling_matches_encode() {
        let seq = b"ACGTACGTACGT";
        let k = 4;
        let rolling: Vec<u64> = RollingKmers::try_new(seq, k).unwrap().flatten().collect();
        let naive: Vec<u64> = seq.windows(k).map(|w| encode(w).unwrap()).collect();
        assert_eq!(rolling, naive);
    }

    #[test]
    fn rolling_skips_n_bearing_windows() {
        let seq = b"ACGTNACGT";
        let k = 4;
        let results: Vec<Option<u64>> = RollingKmers::try_new(seq, k).unwrap().collect();
        assert_eq!(results.len(), 9);
        // First valid k-mer at index k-1=3 (ACGT)
        assert!(results[0].is_none());
        assert!(results[1].is_none());
        assert!(results[2].is_none());
        assert!(results[3].is_some()); // ACGT
        // N at index 4 resets; need 4 more valid bases to recover
        assert!(results[4].is_none()); // N resets
        assert!(results[5].is_none()); // A (valid=1)
        assert!(results[6].is_none()); // C (valid=2)
        assert!(results[7].is_none()); // G (valid=3)
        assert!(results[8].is_some()); // T (valid=4 → ACGT)
    }

    #[test]
    fn rolling_empty_seq() {
        let results: Vec<_> = RollingKmers::try_new(b"", 4).unwrap().collect();
        assert!(results.is_empty());
    }

    #[test]
    fn rolling_seq_shorter_than_k() {
        for k in [1, 31, 32] {
            let seq = if k == 1 {
                b"".as_slice()
            } else {
                b"ACG".as_slice()
            };
            let results: Vec<_> = RollingKmers::try_new(seq, k).unwrap().collect();
            assert_eq!(results.len(), seq.len());
            assert!(results.iter().all(Option::is_none));
        }
    }

    #[test]
    fn rolling_k32_matches_encode() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGA";
        let rolling: Vec<u64> = RollingKmers::try_new(seq, 32).unwrap().flatten().collect();
        assert_eq!(rolling, vec![encode(seq).unwrap()]);
    }

    #[test]
    fn rolling_rejects_zero_k() {
        assert!(matches!(
            RollingKmers::try_new(b"ACGT", 0),
            Err(KmerError::KOutOfRange(0))
        ));
    }

    #[test]
    fn rolling_rejects_k_above_encoding_capacity() {
        assert!(matches!(
            RollingKmers::try_new(b"ACGT", 33),
            Err(KmerError::KOutOfRange(33))
        ));
    }

    #[test]
    fn source_compatible_new_returns_the_iterator() {
        let mut rolling: RollingKmers<'_> = RollingKmers::new(b"A", 1);
        assert_eq!(rolling.next(), Some(Some(0)));
        assert_eq!(rolling.next(), None);
    }

    #[test]
    fn source_compatible_new_fails_loudly_for_invalid_k() {
        for k in [0, 33] {
            assert!(std::panic::catch_unwind(|| RollingKmers::new(b"ACGT", k)).is_err());
        }
    }

    #[test]
    fn rolling_matches_window_encoding_across_boundaries_and_ambiguity_runs() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTNNacgtacgtacgtacgtacgtacgtacgtacgt";
        for k in [1, 31, 32] {
            let rolling: Vec<_> = RollingKmers::try_new(seq, k).unwrap().collect();
            assert_eq!(rolling.len(), seq.len());
            for (end, actual) in rolling.into_iter().enumerate() {
                let expected = if end + 1 < k {
                    None
                } else {
                    encode(&seq[end + 1 - k..=end]).ok()
                };
                assert_eq!(actual, expected, "k={k}, end={end}");
            }
        }
    }

    #[test]
    fn exact_size_tracks_remaining_input_positions() {
        let mut rolling = RollingKmers::try_new(b"ACGTN", 2).unwrap();
        for remaining in (0..=5).rev() {
            assert_eq!(rolling.len(), remaining);
            if remaining > 0 {
                rolling.next();
            }
        }
    }
}
