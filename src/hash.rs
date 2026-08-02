use std::num::NonZeroUsize;

/// Returns the ntHash value of the first complete k-mer, if available.
///
/// Invalid ntHash parameters and sequences shorter than `k` return `None`.
pub fn nthash_one(seq: &[u8], k: usize) -> Option<u64> {
    if seq.len() < k {
        return None;
    }
    let iter = nthash::NtHashIterator::new(seq, k).ok()?;
    iter.into_iter().next()
}

/// Creates an upstream ntHash rolling iterator.
///
/// # Errors
///
/// Returns the upstream [`nthash::Error`] when its sequence or `k` constraints
/// are not satisfied.
pub fn nthash_iter(seq: &[u8], k: usize) -> nthash::Result<nthash::NtHashIterator<'_>> {
    nthash::NtHashIterator::new(seq, k)
}

/// Computes MurmurHash3 x64 128-bit output for `bytes` and a 32-bit seed.
#[must_use]
pub fn murmur3_x64_128(bytes: &[u8], seed: u32) -> u128 {
    murmur3_x64_128_with_seed(bytes, u64::from(seed))
}

/// Reusable canonical DNA-window hasher using the first MurmurHash3 x64 lane.
pub struct CanonicalMurmur64 {
    k: usize,
    seed: u64,
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

impl CanonicalMurmur64 {
    /// Creates a canonical DNA-window hasher with allocation-reusing buffers.
    ///
    /// # Panics
    ///
    /// Panics when scratch allocation fails. Use [`Self::try_new`] when `k`
    /// comes from an external boundary.
    #[must_use]
    pub fn new(k: NonZeroUsize, seed: u64) -> Self {
        Self::try_new(k, seed).expect("canonical Murmur64 scratch allocation failed")
    }

    /// Creates a canonical DNA-window hasher after checking scratch allocation.
    ///
    /// # Errors
    ///
    /// Returns [`crate::KmerError::AllocationFailed`] when buffers for `k`
    /// bases cannot be reserved.
    pub fn try_new(k: NonZeroUsize, seed: u64) -> crate::Result<Self> {
        let k = k.get();
        let mut forward = Vec::new();
        let mut reverse = Vec::new();
        forward
            .try_reserve_exact(k)
            .map_err(|_| crate::KmerError::AllocationFailed { k })?;
        reverse
            .try_reserve_exact(k)
            .map_err(|_| crate::KmerError::AllocationFailed { k })?;
        forward.resize(k, 0);
        reverse.resize(k, 0);
        Ok(Self {
            k,
            seed,
            forward,
            reverse,
        })
    }

    /// Iterates over the complete windows in `sequence`.
    ///
    /// Windows containing bytes outside ASCII A/C/G/T are returned as `None`;
    /// valid input is normalized to uppercase before hashing.
    pub fn hashes<'hasher, 'sequence>(
        &'hasher mut self,
        sequence: &'sequence [u8],
    ) -> CanonicalMurmur64Hashes<'hasher, 'sequence> {
        CanonicalMurmur64Hashes {
            sequence,
            k: self.k,
            seed: self.seed,
            start: 0,
            invalid: 0,
            initialized: false,
            forward: &mut self.forward,
            reverse: &mut self.reverse,
        }
    }
}

/// Hashes produced by [`CanonicalMurmur64::hashes`].
pub struct CanonicalMurmur64Hashes<'hasher, 'sequence> {
    sequence: &'sequence [u8],
    k: usize,
    seed: u64,
    start: usize,
    invalid: usize,
    initialized: bool,
    forward: &'hasher mut [u8],
    reverse: &'hasher mut [u8],
}

impl Iterator for CanonicalMurmur64Hashes<'_, '_> {
    type Item = Option<u64>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start > self.sequence.len().saturating_sub(self.k) {
            return None;
        }

        if self.initialized {
            self.invalid -= usize::from(normalize(self.sequence[self.start - 1]).is_none());
            self.invalid +=
                usize::from(normalize(self.sequence[self.start + self.k - 1]).is_none());
        } else {
            self.invalid = self.sequence[..self.k]
                .iter()
                .filter(|&&base| normalize(base).is_none())
                .count();
            self.initialized = true;
        }

        let value = if self.invalid == 0 {
            let window = &self.sequence[self.start..self.start + self.k];
            for (index, &base) in window.iter().enumerate() {
                let base = normalize(base).expect("window validity was established");
                self.forward[index] = base;
                self.reverse[self.k - index - 1] = complement(base);
            }
            let forward: &[u8] = self.forward;
            let reverse: &[u8] = self.reverse;
            let canonical = forward.min(reverse);
            Some(murmur3_x64_128_with_seed(canonical, self.seed) as u64)
        } else {
            None
        };
        self.start += 1;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let total = self.sequence.len().saturating_sub(self.k - 1);
        let remaining = total.saturating_sub(self.start);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CanonicalMurmur64Hashes<'_, '_> {}

const fn normalize(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(b'A'),
        b'C' | b'c' => Some(b'C'),
        b'G' | b'g' => Some(b'G'),
        b'T' | b't' => Some(b'T'),
        _ => None,
    }
}

const fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => unreachable!(),
    }
}

fn murmur3_x64_128_with_seed(bytes: &[u8], seed: u64) -> u128 {
    const C1: u64 = 0x87c3_7b91_1142_53d5;
    const C2: u64 = 0x4cf5_ad43_2745_937f;

    let mut h1 = seed;
    let mut h2 = seed;
    let mut chunks = bytes.chunks_exact(16);
    for chunk in &mut chunks {
        let mut first = [0; 8];
        let mut second = [0; 8];
        first.copy_from_slice(&chunk[..8]);
        second.copy_from_slice(&chunk[8..]);
        let mut k1 = u64::from_le_bytes(first);
        let mut k2 = u64::from_le_bytes(second);

        k1 = k1.wrapping_mul(C1).rotate_left(31).wrapping_mul(C2);
        h1 ^= k1;
        h1 = h1
            .rotate_left(27)
            .wrapping_add(h2)
            .wrapping_mul(5)
            .wrapping_add(0x52dc_e729);

        k2 = k2.wrapping_mul(C2).rotate_left(33).wrapping_mul(C1);
        h2 ^= k2;
        h2 = h2
            .rotate_left(31)
            .wrapping_add(h1)
            .wrapping_mul(5)
            .wrapping_add(0x3849_5ab5);
    }

    let tail = chunks.remainder();
    let mut k1 = 0_u64;
    let mut k2 = 0_u64;
    for (index, &byte) in tail.iter().take(8).enumerate() {
        k1 |= u64::from(byte) << (index * 8);
    }
    for (index, &byte) in tail.iter().skip(8).enumerate() {
        k2 |= u64::from(byte) << (index * 8);
    }
    if tail.len() > 8 {
        h2 ^= k2.wrapping_mul(C2).rotate_left(33).wrapping_mul(C1);
    }
    if !tail.is_empty() {
        h1 ^= k1.wrapping_mul(C1).rotate_left(31).wrapping_mul(C2);
    }

    let length = bytes.len() as u64;
    h1 ^= length;
    h2 ^= length;
    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);
    h1 = mix(h1);
    h2 = mix(h2);
    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);
    (u128::from(h2) << 64) | u128::from(h1)
}

const fn mix(mut value: u64) -> u64 {
    value ^= value >> 33;
    value = value.wrapping_mul(0xff51_afd7_ed55_8ccd);
    value ^= value >> 33;
    value = value.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    value ^ (value >> 33)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nthash_one_emits_a_hash_for_first_kmer() {
        let h = nthash_one(b"ACGTACGT", 4);
        assert!(h.is_some());
    }

    #[test]
    fn nthash_iter_n_minus_k_plus_1() {
        let seq = b"ACGTACGTACGT";
        let count = nthash_iter(seq, 5).unwrap().count();
        assert_eq!(count, seq.len() - 5 + 1);
    }

    #[test]
    fn nthash_canonical_is_rc_invariant() {
        let h1 = nthash_one(b"AAAATTTT", 8).unwrap();
        let h2 = nthash_one(b"AAAATTTT", 8).unwrap();
        assert_eq!(h1, h2);
    }

    #[test]
    fn murmur_changes_with_seed() {
        let h0 = murmur3_x64_128(b"ACGT", 0);
        let h1 = murmur3_x64_128(b"ACGT", 1);
        assert_ne!(h0, h1);
    }

    #[test]
    fn murmur_is_deterministic() {
        let h0 = murmur3_x64_128(b"ACGT", 42);
        let h0_again = murmur3_x64_128(b"ACGT", 42);
        assert_eq!(h0, h0_again);
    }

    #[test]
    fn murmur_matches_sourmash_lane() {
        assert_eq!(
            murmur3_x64_128(b"ACG", 42) as u64,
            1_731_421_407_650_554_201
        );
        let vectors = [
            (b"".as_slice(), 17_305_828_677_633_410_339),
            (b"123456789abcdef".as_slice(), 1_372_604_790_551_148_737),
            (b"123456789abcdef1".as_slice(), 2_544_995_835_007_078_785),
            (b"123456789abcdef12".as_slice(), 1_591_540_619_117_953_758),
            (
                b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_slice(),
                857_093_537_268_035_581,
            ),
        ];
        for (input, expected) in vectors {
            assert_eq!(murmur3_x64_128_with_seed(input, 42) as u64, expected);
        }
    }

    #[test]
    fn canonical_murmur_matches_sourmash_wide_seed() {
        let k = NonZeroUsize::new(3).unwrap();
        let mut hasher = CanonicalMurmur64::new(k, 4_294_967_297);
        let hashes: Vec<_> = hasher.hashes(b"ACG").collect();
        assert_eq!(hashes, vec![Some(4_980_406_537_561_696_718)]);
        assert_eq!(
            murmur3_x64_128_with_seed(b"123456789abcdef12", u64::MAX) as u64,
            12_114_998_237_397_089_840
        );
    }

    #[test]
    fn canonical_murmur_normalizes_strands_and_case() {
        let k = NonZeroUsize::new(4).unwrap();
        let mut hasher = CanonicalMurmur64::new(k, 42);
        let forward: Vec<_> = hasher.hashes(b"acgtTAGC").collect();
        let reverse: Vec<_> = hasher.hashes(b"GCTAacgt").collect();
        assert_eq!(forward, reverse.into_iter().rev().collect::<Vec<_>>());
    }

    #[test]
    fn canonical_murmur_marks_only_affected_windows_invalid() {
        let k = NonZeroUsize::new(3).unwrap();
        let mut hasher = CanonicalMurmur64::new(k, 42);
        let hashes: Vec<_> = hasher.hashes(b"ACGNUTACG").collect();
        assert_eq!(hashes.len(), 7);
        assert!(hashes[0].is_some());
        assert!(hashes[1..5].iter().all(Option::is_none));
        assert!(hashes[5].is_some());
        assert!(hashes[6].is_some());
    }

    #[test]
    fn canonical_murmur_supports_long_k() {
        let k = NonZeroUsize::new(51).unwrap();
        let sequence = b"ACGTTGCATGTCGCATGATGCATGAGAGCTACGTACGTACGTACGTACGTACGT";
        let mut hasher = CanonicalMurmur64::new(k, 42);
        let hashes: Vec<_> = hasher.hashes(sequence).collect();
        assert_eq!(hashes.len(), sequence.len() - 50);
        assert!(hashes.into_iter().all(|hash| hash.is_some()));
    }

    #[test]
    fn canonical_murmur_reuses_buffers_between_records() {
        let mut hasher = CanonicalMurmur64::new(NonZeroUsize::new(3).unwrap(), 42);
        let first: Vec<_> = hasher.hashes(b"ACG").collect();
        let second: Vec<_> = hasher.hashes(b"acg").collect();
        assert_eq!(first, second);
    }

    #[test]
    fn canonical_murmur_reports_impossible_scratch_allocation() {
        let k = NonZeroUsize::new(usize::MAX).unwrap();
        assert!(matches!(
            CanonicalMurmur64::try_new(k, 42),
            Err(crate::KmerError::AllocationFailed { k: usize::MAX })
        ));
    }
}
