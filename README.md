# rsomics-kmer

`rsomics-kmer` provides checked two-bit DNA k-mer primitives for the
`rsomics-*` product family:

- encoding, decoding, reverse complement, and lexicographic canonicalization;
- fixed-window and linear-time rolling iteration;
- exact k-mer counting;
- ntHash and MurmurHash3 adapters;
- allocation-reusing canonical Murmur64 windows for arbitrary non-zero `k`.

The codec supports `k` in `1..=32`. A `Kmer` stores A/C/G/T as `00/01/10/11`
in the low `2 * k` bits. Values passed to codec functions must be normalized:
bits above that representation must be zero. Values returned by this crate are
always normalized.

## API example

```rust
use rsomics_kmer::{
    RollingKmers, canonical, decode, encode, reverse_complement, try_decode,
};

let encoded = encode(b"ACGT").unwrap();
assert_eq!(decode(encoded, 4), b"ACGT");

let reverse = reverse_complement(encoded, 4);
assert_eq!(canonical(encoded, 4), encoded.min(reverse));

let kmers: Vec<_> = RollingKmers::new(b"ACGTAC", 4).flatten().collect();
assert_eq!(kmers.len(), 3);

assert!(RollingKmers::try_new(b"ACGT", 0).is_err());
assert!(try_decode(0, 33).is_err());

let mut hasher = rsomics_kmer::CanonicalMurmur64::try_new(
    std::num::NonZeroUsize::new(51).unwrap(),
    42,
)?;
let hashes: Vec<_> = hasher
    .hashes(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
    .collect();
assert_eq!(hashes.len(), 2);
assert_eq!(hasher.hashes(b"ACG").next(), None);
# Ok::<(), rsomics_kmer::KmerError>(())
```

`RollingKmers` yields one `Option<Kmer>` per input base. Positions before the
first complete window are `None`; non-ACGT bytes reset the rolling state, so
windows spanning ambiguity are also `None`. Flattening the iterator yields only
valid k-mers.

`CanonicalMurmur64` yields one item per complete window. Inputs shorter than
`k`, including empty input, yield no items. Complete windows containing
non-ACGT bytes yield `None`; valid windows are normalized to uppercase and
canonicalized across both strands before hashing.

The infallible `decode`, `reverse_complement`, `canonical`, and
`RollingKmers::new` APIs fail loudly when their documented representation or
length preconditions are violated. Their `try_*` counterparts return
`KmerError` for runtime-validated input.

## Validation and benchmark

```text
cargo test --all-targets
cargo bench --bench codec
```

The benchmark covers the rolling scanner on a one-megabase sequence with
ambiguity runs, batched reverse-complement/canonicalization, and the canonical
Murmur64 consumer path with reusable scratch buffers.

## License

Licensed under either Apache-2.0 or MIT, at your option.
