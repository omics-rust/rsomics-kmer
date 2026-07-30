use crate::encode::{Kmer, encode};
use crate::{KmerError, Result};

/// An iterator over every fixed-width sequence window.
///
/// Unlike [`crate::RollingKmers`], this iterator yields only the
/// `n - k + 1` complete windows, and reports a non-ACGT window as an error.
pub struct KmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
    current: Option<Kmer>,
    canonical: bool,
}

impl<'a> KmerIter<'a> {
    /// Creates a window iterator, optionally yielding canonical encodings.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::KOutOfRange`] when `k` is outside `1..=32`, or
    /// [`KmerError::SeqTooShort`] when `seq` is shorter than `k`.
    pub fn new(seq: &'a [u8], k: usize, canonical: bool) -> Result<Self> {
        if !(1..=32).contains(&k) {
            return Err(KmerError::KOutOfRange(k));
        }
        if seq.len() < k {
            return Err(KmerError::SeqTooShort { len: seq.len(), k });
        }
        Ok(Self {
            seq,
            k,
            pos: 0,
            current: None,
            canonical,
        })
    }
}

impl Iterator for KmerIter<'_> {
    type Item = Result<Kmer>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos + self.k > self.seq.len() {
            return None;
        }
        let window = &self.seq[self.pos..self.pos + self.k];
        self.pos += 1;
        let kmer = match encode(window) {
            Ok(k) => k,
            Err(e) => {
                self.current = None;
                return Some(Err(e));
            }
        };
        let out = if self.canonical {
            crate::encode::canonical(kmer, self.k)
        } else {
            kmer
        };
        self.current = Some(out);
        Some(Ok(out))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::{canonical, decode};

    #[test]
    fn iter_yields_n_minus_k_plus_1_kmers() {
        let seq = b"ACGTACGT";
        let it = KmerIter::new(seq, 4, false).unwrap();
        let kmers: Vec<_> = it.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(kmers.len(), seq.len() - 4 + 1);
        let first = decode(kmers[0], 4);
        let last = decode(kmers[4], 4);
        assert_eq!(&first, b"ACGT");
        assert_eq!(&last, b"ACGT");
    }

    #[test]
    fn iter_canonical_is_idempotent_under_rc() {
        let fwd: Vec<u64> = KmerIter::new(b"ACGTACGT", 4, true)
            .unwrap()
            .collect::<Result<_>>()
            .unwrap();
        let rev_comp = b"ACGTACGT";
        let rc: Vec<u64> = KmerIter::new(rev_comp, 4, true)
            .unwrap()
            .collect::<Result<_>>()
            .unwrap();
        let mut a = fwd;
        let mut b = rc;
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn iter_too_short_seq_errors() {
        assert!(matches!(
            KmerIter::new(b"AC", 4, false),
            Err(KmerError::SeqTooShort { len: 2, k: 4 })
        ));
    }

    #[test]
    fn iter_propagates_non_acgt_error() {
        let seq = b"ACGTNAC";
        let mut it = KmerIter::new(seq, 4, false).unwrap();
        assert!(it.next().unwrap().is_ok());
        let third = it.next().unwrap();
        assert!(matches!(third, Err(KmerError::NonAcgt { .. })));
    }

    #[test]
    fn iterator_boundaries_validate_and_round_trip() {
        for k in [0, 33] {
            assert!(matches!(
                KmerIter::new(b"ACGT", k, true),
                Err(KmerError::KOutOfRange(actual)) if actual == k
            ));
        }
        for k in [1, 31, 32] {
            let seq: Vec<_> = (0..k).map(|i| b"ACGT"[i % 4]).collect();
            let values: Vec<_> = KmerIter::new(&seq, k, true)
                .unwrap()
                .collect::<Result<_>>()
                .unwrap();
            assert_eq!(values.len(), 1);
            assert_eq!(values[0], canonical(encode(&seq).unwrap(), k));
        }
    }
}
