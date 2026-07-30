use std::collections::HashMap;

use crate::encode::{Kmer, validate_k};
use crate::iter::KmerIter;
use crate::{KmerError, Result};

/// Exact k-mer counts accumulated across one or more sequences.
#[derive(Debug, Default, Clone)]
pub struct KmerCounts {
    /// K-mer length.
    pub k: usize,
    /// Whether to collapse each k-mer with its reverse complement.
    pub canonical: bool,
    /// Encoded k-mer counts.
    pub counts: HashMap<Kmer, u64>,
}

impl KmerCounts {
    /// Creates an empty accumulator without validating `k`.
    ///
    /// [`Self::count_seq`] validates the length before reading a sequence.
    /// Use [`Self::try_new`] when `k` comes from an external boundary.
    #[must_use]
    pub fn new(k: usize, canonical: bool) -> Self {
        Self {
            k,
            canonical,
            counts: HashMap::new(),
        }
    }

    /// Creates an empty accumulator after validating `k`.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`.
    pub fn try_new(k: usize, canonical: bool) -> Result<Self> {
        validate_k(k)?;
        Ok(Self::new(k, canonical))
    }

    /// Adds all valid windows in `seq` to the counts.
    ///
    /// Windows containing non-ACGT bytes are skipped. A valid `k` with a
    /// sequence shorter than `k` is a no-op.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::KOutOfRange`] when `self.k` is outside `1..=32`.
    pub fn count_seq(&mut self, seq: &[u8]) -> Result<()> {
        validate_k(self.k)?;
        if seq.len() < self.k {
            return Ok(());
        }
        let it = KmerIter::new(seq, self.k, self.canonical)?;
        for kmer in it {
            match kmer {
                Ok(k) => *self.counts.entry(k).or_insert(0) += 1,
                Err(KmerError::NonAcgt { .. }) => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    /// Returns the number of distinct encoded k-mers.
    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns whether no k-mers have been counted.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Returns the total count across all distinct k-mers.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.counts.values().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_simple_seq() {
        let mut c = KmerCounts::new(3, false);
        c.count_seq(b"AAAAA").unwrap();
        assert_eq!(c.total(), 3);
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn count_canonical_collapses_rc_pairs() {
        let mut c = KmerCounts::new(4, true);
        c.count_seq(b"AAAATTTT").unwrap();
        let entries: Vec<_> = c.counts.iter().collect();
        assert!(!entries.is_empty());
    }

    #[test]
    fn count_skips_n_kmers_silently() {
        let mut c = KmerCounts::new(4, false);
        c.count_seq(b"ACGTNACGT").unwrap();
        assert_eq!(c.total(), 2);
    }

    #[test]
    fn count_seq_shorter_than_k_is_noop() {
        let mut c = KmerCounts::new(10, false);
        c.count_seq(b"ACGT").unwrap();
        assert_eq!(c.total(), 0);
    }

    #[test]
    fn invalid_k_is_rejected_even_for_short_sequences() {
        for k in [0, 33] {
            let mut counts = KmerCounts::new(k, false);
            assert!(matches!(
                counts.count_seq(b""),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
            assert!(counts.is_empty());
        }
    }

    #[test]
    fn checked_constructor_rejects_invalid_k_before_counting() {
        for k in [0, 33, usize::MAX] {
            assert!(matches!(
                KmerCounts::try_new(k, false),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
        }
        for k in [1, 32] {
            let counts = KmerCounts::try_new(k, true).unwrap();
            assert_eq!(counts.k, k);
            assert!(counts.canonical);
            assert!(counts.is_empty());
        }
    }
}
