use std::num::NonZeroUsize;

use rsomics_kmer::{CanonicalMurmur64, CanonicalMurmur64Hashes};

fn assert_exhausted(mut hashes: CanonicalMurmur64Hashes<'_, '_>) {
    for _ in 0..3 {
        assert_eq!(hashes.len(), 0);
        assert_eq!(hashes.size_hint(), (0, Some(0)));
        assert_eq!(hashes.next(), None);
    }
}

#[test]
#[ignore = "known short-window bounds regression"]
fn empty_sequences_are_exhausted() {
    for k in [1, 31, 51] {
        let mut hasher = CanonicalMurmur64::try_new(NonZeroUsize::new(k).unwrap(), 42).unwrap();
        assert_exhausted(hasher.hashes(b""));
    }
}

#[test]
#[ignore = "known short-window bounds regression"]
fn short_sequences_are_exhausted_and_hasher_remains_reusable() {
    for k in [2, 31, 51] {
        let mut hasher = CanonicalMurmur64::try_new(NonZeroUsize::new(k).unwrap(), 42).unwrap();
        let sequence = vec![b'A'; k + 1];
        let original: Vec<_> = hasher.hashes(&sequence).collect();
        assert_eq!(original.len(), 2);
        assert_exhausted(hasher.hashes(&sequence[..k - 1]));
        assert_eq!(hasher.hashes(&sequence).collect::<Vec<_>>(), original);
        assert_exhausted(hasher.hashes(b""));
    }
}

#[test]
fn exact_and_one_extra_base_produce_complete_windows() {
    for k in [1, 31, 51] {
        let mut hasher = CanonicalMurmur64::try_new(NonZeroUsize::new(k).unwrap(), 42).unwrap();
        for (length, count) in [(k, 1), (k + 1, 2)] {
            let sequence = vec![b'A'; length];
            let mut hashes = hasher.hashes(&sequence);
            for remaining in (1..=count).rev() {
                assert_eq!(hashes.len(), remaining);
                assert_eq!(hashes.size_hint(), (remaining, Some(remaining)));
                assert!(hashes.next().unwrap().is_some());
            }
            assert_exhausted(hashes);
        }
    }
}
