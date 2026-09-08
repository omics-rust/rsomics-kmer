use std::error::Error;
use std::fs::{self, File};
use std::hint::black_box;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::time::Instant;

const ORDERS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [2, 1, 0],
    [1, 2, 0],
    [0, 2, 1],
    [2, 0, 1],
    [1, 0, 2],
];
const NAMES: [&str; 3] = ["baseline", "reference", "candidate"];

enum Hashes<'hasher, 'sequence> {
    Baseline(baseline::CanonicalMurmur64Hashes<'hasher, 'sequence>),
    Reference(reference::CanonicalMurmur64Hashes<'hasher, 'sequence>),
    Candidate(candidate::CanonicalMurmur64Hashes<'hasher, 'sequence>),
}

fn fixture(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| if i % 4_093 == 0 { b'N' } else { b"ACGT"[i % 4] })
        .collect()
}

#[inline(never)]
fn checksum(hashes: &mut dyn Iterator<Item = Option<u64>>) -> u64 {
    hashes.flatten().fold(0, |acc, hash| acc ^ hash)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let output = PathBuf::from(args.next().ok_or("expected an output directory")?);
    if args.next().is_some() {
        return Err("unexpected additional argument".into());
    }
    let sequence = fixture(1_048_576);
    fs::write(output.join("input.bin"), &sequence)?;
    let k = NonZeroUsize::new(31).unwrap();
    let mut baseline = baseline::CanonicalMurmur64::try_new(k, 42)?;
    let mut reference = reference::CanonicalMurmur64::try_new(k, 42)?;
    let mut candidate = candidate::CanonicalMurmur64::try_new(k, 42)?;
    let expected = checksum(&mut baseline.hashes(&sequence));
    let mut observations = Vec::with_capacity(66 * 3);

    for round in 0..66 {
        for (position, &variant) in ORDERS[round % ORDERS.len()].iter().enumerate() {
            let mut selected = match variant {
                0 => Hashes::Baseline(baseline.hashes(black_box(&sequence))),
                1 => Hashes::Reference(reference.hashes(black_box(&sequence))),
                2 => Hashes::Candidate(candidate.hashes(black_box(&sequence))),
                _ => unreachable!(),
            };
            let hashes: &mut dyn Iterator<Item = Option<u64>> = match &mut selected {
                Hashes::Baseline(hashes) => hashes,
                Hashes::Reference(hashes) => hashes,
                Hashes::Candidate(hashes) => hashes,
            };
            let started = Instant::now();
            let sum = checksum(black_box(hashes));
            let elapsed = started.elapsed().as_nanos();
            black_box(sum);
            if sum != expected {
                return Err(format!("{} produced a different checksum", NAMES[variant]).into());
            }
            observations.push((round, position, variant, elapsed, sum));
        }
    }

    let mut csv = BufWriter::new(File::create_new(output.join("observations.csv"))?);
    writeln!(csv, "phase,round,position,variant,wall_ns,checksum")?;
    for (round, position, variant, elapsed, sum) in observations {
        let (phase, index) = if round < 6 {
            ("warmup", round)
        } else {
            ("measurement", round - 6)
        };
        writeln!(
            csv,
            "{phase},{index},{position},{},{elapsed},{sum}",
            NAMES[variant]
        )?;
    }
    csv.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schedule_balances_every_position() {
        let mut counts = [[0; 3]; 3];
        for order in ORDERS {
            let mut variants = order;
            variants.sort();
            assert_eq!(variants, [0, 1, 2]);
            for (position, variant) in order.into_iter().enumerate() {
                counts[position][variant] += 1;
            }
        }
        assert_eq!(counts, [[2; 3]; 3]);
    }

    #[test]
    fn fixture_retains_the_existing_codec_workload() {
        assert_eq!(fixture(8), b"NCGTACGT");
        let sequence = fixture(4_095);
        assert_eq!(&sequence[4_092..], b"ANG");
    }

    #[test]
    fn checksum_excludes_invalid_windows() {
        assert_eq!(checksum(&mut [Some(4), None, Some(8)].into_iter()), 12);
    }
}
