#![allow(clippy::cast_possible_truncation, clippy::must_use_candidate)]
#![warn(missing_docs)]

//! Checked two-bit DNA k-mer encoding, iteration, hashing, and counting.
//!
//! ```
//! use rsomics_kmer::{RollingKmers, decode, encode, try_decode};
//!
//! let encoded = encode(b"ACGT").unwrap();
//! assert_eq!(decode(encoded, 4), b"ACGT");
//! assert_eq!(try_decode(encoded, 4).unwrap(), b"ACGT");
//!
//! let kmers: Vec<_> = RollingKmers::new(b"ACGTAC", 4).flatten().collect();
//! assert_eq!(kmers.len(), 3);
//! ```

/// K-mer counting.
pub mod count;
/// Two-bit DNA encoding and canonicalization.
pub mod encode;
/// K-mer hashing adapters.
pub mod hash;
/// Window-by-window k-mer iteration.
pub mod iter;
/// Linear-time rolling k-mer iteration.
pub mod roll;

pub use count::KmerCounts;
pub use encode::{
    Kmer, base_bits, canonical, decode, encode, reverse_complement, try_canonical, try_decode,
    try_reverse_complement,
};
pub use hash::{
    CanonicalMurmur64, CanonicalMurmur64Hashes, murmur3_x64_128, nthash_iter, nthash_one,
};
pub use iter::KmerIter;
pub use roll::RollingKmers;

/// Errors from checked k-mer operations.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum KmerError {
    /// `k` cannot be represented by a non-empty 64-bit two-bit k-mer.
    #[error("k must be in 1..=32 (got {0})")]
    KOutOfRange(usize),
    /// A sequence is shorter than the requested k-mer length.
    #[error("sequence shorter than k: len={len}, k={k}")]
    SeqTooShort {
        /// Sequence length in bases.
        len: usize,
        /// Requested k-mer length.
        k: usize,
    },
    /// A sequence contains a byte outside ASCII A/C/G/T, case-insensitively.
    #[error("non-ACGT base at position {pos}: {byte:?}")]
    NonAcgt {
        /// Zero-based byte position.
        pos: usize,
        /// Rejected byte.
        byte: u8,
    },
    /// Bits above the low `2 * k` representation are set.
    #[error("encoded k-mer {kmer:#018x} has bits outside its {k}-base representation")]
    NonNormalizedKmer {
        /// Rejected encoded value.
        kmer: Kmer,
        /// Declared k-mer length.
        k: usize,
    },
    /// Scratch storage for a requested k-mer length cannot be reserved.
    #[error("cannot allocate canonical hash buffers for k={k}")]
    AllocationFailed {
        /// Requested k-mer length.
        k: usize,
    },
}

/// Result type for checked k-mer operations.
pub type Result<T> = std::result::Result<T, KmerError>;
