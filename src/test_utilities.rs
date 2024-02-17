//! Test cases and test utility functions.
//!

use std::path::PathBuf;

use crate::{
    create_granges_with_seqlens,
    error::GRangesError,
    prelude::{GRanges, VecRangesIndexed},
    ranges::{
        coitrees::COITrees,
        vec::{VecRanges, VecRangesEmpty},
        RangeEmpty,
    },
    Position,
};
use indexmap::IndexMap;
use rand::{seq::SliceRandom, thread_rng, Rng};

/// Get the path to the `grange` command line tool after a build.
/// This is used for integration tests and benchmarks.
pub fn granges_binary_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push(if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    });
    path.push("granges");
    path
}

// Stochastic test ranges defaults
//
// This is the random number of range to use in tests.
// The tradeoff is catching stochastic errors vs test time.
pub const NRANDOM_RANGES: usize = 10000;

// range length
pub const MIN_LEN: Position = 1;
pub const MAX_LEN: Position = 10000;

// number of chromosome sequences
pub const NCHROM: usize = 22;

// chromosome sizes
pub const MIN_CHROM_LEN: Position = 50_000_000;
pub const MAX_CHROM_LEN: Position = 250_000_000;

/// Build a random range start/end on a sequence of `max_len`.
/// 0-indexed, right exclusive
pub fn random_range(chrom_len: Position) -> (Position, Position) {
    let mut rng = thread_rng();
    let len = rng.gen_range(MIN_LEN..MAX_LEN);
    let start = rng.gen_range(0..chrom_len - len + 1);
    (start, start + len)
}

/// Build random sequence lengths
pub fn random_seqlen() -> Position {
    let mut rng = thread_rng();
    rng.gen_range(MIN_CHROM_LEN..=MAX_CHROM_LEN)
}

/// Sample a random chromosome
pub fn random_chrom() -> String {
    let mut rng = thread_rng();
    format!("chr{}", rng.gen_range(1..NCHROM + 1))
}

/// Build random [`VecRanges`]
pub fn random_vecranges(n: usize) -> VecRanges<RangeEmpty> {
    let seqlen = random_seqlen();
    let mut vr = VecRanges::new(seqlen);
    for _i in 0..n {
        let (start, end) = random_range(seqlen);
        let range = RangeEmpty::new(start, end);
        vr.push_range(range);
    }
    vr
}

/// Build a random [`GRanges<VecRangesEmpty, ()`] using a set of
/// sequence lengths.
pub fn random_granges(
    seqlens: &IndexMap<String, Position>,
    num: usize,
) -> Result<GRanges<VecRangesEmpty, ()>, GRangesError> {
    let mut rng = thread_rng();

    let mut gr: GRanges<VecRangesEmpty, ()> = GRanges::new_vec(seqlens);

    let seqnames: Vec<String> = seqlens.keys().cloned().collect();
    for _ in 0..num {
        let seqname = seqnames.choose(&mut rng).unwrap();
        let chrom_len = *seqlens
            .get(seqname)
            .ok_or_else(|| GRangesError::MissingSequence(seqname.clone()))?;
        let (start, end) = random_range(chrom_len);
        gr.push_range_empty(seqname, start, end)?;
    }
    Ok(gr)
}

/// Build random [`COITrees`] from a random [`VecRanges`].
pub fn random_coitrees() -> COITrees<()> {
    let vr = random_vecranges(100);
    let cr: COITrees<()> = vr.into();
    cr
}

/// Range test case #1
///
/// Ranges:
///   - chr1:
///      (0, 5, Some(1.1))
///      (4, 7, Some(8.1))
///      (10, 17, Some(10.1))
///   - chr2:
///      (10, 20, Some(3.7))
///      (18, 32, Some(1.1))
///
/// Seqlens: { "chr1" => 30, "chr2" => 100 }
///
/// Sum of all elements: 24.1
pub fn granges_test_case_01() -> GRanges<VecRangesIndexed, Vec<f64>> {
    create_granges_with_seqlens!(VecRangesIndexed, Vec<f64>, {
        "chr1" => [(0, 5, 1.1), (4, 7, 8.1), (10, 17, 10.1)],
        "chr2" => [(10, 20, 3.7), (18, 32, 1.1)]
    }, seqlens: { "chr1" => 30, "chr2" => 100 })
}
