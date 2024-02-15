//! Test cases and test utility functions.
//!

use indexmap::IndexMap;
use rand::{thread_rng, Rng};
use crate::{Position, ranges::{RangeIndexed, vec::VecRanges, RangeEmpty}};

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
    let rand_len = rng.gen_range(MIN_CHROM_LEN..=MAX_CHROM_LEN);
    rand_len
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



 
