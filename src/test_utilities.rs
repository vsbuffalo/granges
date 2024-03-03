//! Test cases and test utility functions.
//!

use std::{io::BufRead, path::PathBuf};

use crate::{
    commands::granges_random_bed,
    create_granges_with_seqlens,
    error::GRangesError,
    granges::GRangesEmpty,
    io::{parsers::bed::Bed5Addition, InputStream},
    prelude::{GRanges, VecRangesIndexed},
    ranges::{
        coitrees::COITrees,
        vec::{VecRanges, VecRangesEmpty},
        RangeEmpty,
    },
    Position,
};
#[cfg(feature = "ndarray")]
use genomap::GenomeMap;
use indexmap::IndexMap;
#[cfg(feature = "ndarray")]
use ndarray::{Array1, Array2};
use rand::{distributions::Uniform, seq::SliceRandom, thread_rng, Rng};
use tempfile::{Builder, NamedTempFile};

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
pub const MAX_LEN: Position = 1000;

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

/// Build a random [`GRangesEmpty`] using a set of
/// sequence lengths.
pub fn random_granges(
    seqlens: &IndexMap<String, Position>,
    num: usize,
) -> Result<GRangesEmpty<VecRangesEmpty>, GRangesError> {
    let mut rng = thread_rng();

    let mut gr = GRangesEmpty::new_vec(seqlens);

    let seqnames: Vec<String> = seqlens.keys().cloned().collect();
    for _ in 0..num {
        let seqname = seqnames.choose(&mut rng).unwrap();
        let chrom_len = *seqlens
            .get(seqname)
            .ok_or_else(|| GRangesError::MissingSequence(seqname.clone()))?;
        let (start, end) = random_range(chrom_len);
        gr.push_range(seqname, start, end)?;
    }
    Ok(gr)
}

/// Generate random strings, e.g. for mock feature names.
fn generate_random_string(n: usize) -> String {
    let mut rng = rand::thread_rng();
    let letters: Vec<char> = ('a'..='z').collect();
    let letters_dist = Uniform::from(0..letters.len());

    (0..n).map(|_| letters[rng.sample(letters_dist)]).collect()
}

/// Generate a random float value, e.g. for a mock BED "score".
fn generate_random_uniform(start: f64, end: f64) -> f64 {
    let mut rng = rand::thread_rng();
    let uniform = Uniform::new(start, end); // Specify the range
    rng.sample(uniform)
}

/// Build a random [`GRanges`] using a set of sequence lengths,
/// with BED5 like data.
pub fn random_granges_mock_bed5(
    seqlens: &IndexMap<String, Position>,
    num: usize,
) -> Result<GRanges<VecRangesIndexed, Vec<Bed5Addition>>, GRangesError> {
    let mut rng = thread_rng();

    let mut gr = GRanges::new_vec(seqlens);

    let seqnames: Vec<String> = seqlens.keys().cloned().collect();
    for _ in 0..num {
        let seqname = seqnames.choose(&mut rng).unwrap();
        let chrom_len = *seqlens
            .get(seqname)
            .ok_or_else(|| GRangesError::MissingSequence(seqname.clone()))?;
        let (start, end) = random_range(chrom_len);
        let bed5_cols = Bed5Addition {
            name: generate_random_string(8),
            score: Some(generate_random_uniform(0.0, 1.0)),
        };
        gr.push_range(seqname, start, end, bed5_cols)?;
    }
    Ok(gr)
}

/// Build random [`COITrees`] from a random [`VecRanges`].
pub fn random_coitrees() -> COITrees<()> {
    let vr = random_vecranges(100);
    let cr: COITrees<()> = vr.into();
    cr
}

/// Get a temporary file with a .bed suffix.
pub fn temp_bedfile() -> NamedTempFile {
    Builder::new()
        .suffix(".bed")
        .tempfile()
        .expect("Failed to create temp file")
}

/// Create a random BED3 file based on the hg38 sequence lengths, and write to disk.
pub fn random_bed3file(length: usize) -> NamedTempFile {
    let temp_bedfile = temp_bedfile();
    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        length,
        Some(temp_bedfile.path()),
        true,
        false,
    )
    .expect("could not generate random BED file");
    temp_bedfile
}

/// Create a random BED5 file based on the hg38 sequence lengths, and write to disk.
///
/// The feature names are random lowercase characters, and the scores are
/// random floats.
pub fn random_bed5file(length: usize) -> NamedTempFile {
    let temp_bedfile = temp_bedfile();
    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        length,
        Some(temp_bedfile.path()),
        true,
        true,
    )
    .expect("could not generate random BED file");
    temp_bedfile
}

// "head" a file -- useful for peaking at temp files in failed tests.
pub fn head_file(filepath: impl Into<PathBuf>) {
    let filepath = filepath.into();
    let file = InputStream::new(filepath);
    let mut reader = file.reader().unwrap();
    for _ in 0..10 {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        dbg!(&line);
    }
}

/// Range test case #1
///
/// This also has a corresponding reference sequence
/// test file in: `tests_data/sequences/test_case_01.fa.gz`
///
/// Ranges:
///   - chr1:
///      (0, 5, 1.1)
///      (4, 7, 8.1)
///      (10, 17, 10.1)
///   - chr2:
///      (10, 20, 3.7)
///      (18, 32, 1.1)
///
/// Seqlens: { "chr1" => 30, "chr2" => 100 }
///
/// Sum of all elements: 24.1
///
pub fn granges_test_case_01() -> GRanges<VecRangesIndexed, Vec<f64>> {
    create_granges_with_seqlens!(VecRangesIndexed, Vec<f64>, {
        "chr1" => [(0, 5, 1.1), (4, 7, 8.1), (10, 17, 10.1)],
        "chr2" => [(10, 20, 3.7), (18, 32, 1.1)]
    }, seqlens: { "chr1" => 30, "chr2" => 100 })
}

/// Range test case #2
///
/// Ranges:
///   - chr1:
///      (30, 50, Some(1.1))
///   - chr2:
///      (100, 200, Some(3.7))
///      (250, 300, Some(1.1))
///
/// Seqlens: { "chr1" => 50, "chr2" => 300 }
pub fn granges_test_case_02() -> GRanges<VecRangesIndexed, Vec<f64>> {
    create_granges_with_seqlens!(VecRangesIndexed, Vec<f64>, {
        "chr1" => [(30, 50, 1.1)],
        "chr2" => [(100, 200, 3.7), (250, 300, 1.1)]
    }, seqlens: { "chr1" => 50, "chr2" => 300 })
}

/// Range test case #3 (modified from #1)
///
/// This is a test case with Bed5 data. This is to test
/// merging iterators, with grouping by feature name.
///
/// This matches `tests_data/test_case_03.bed`.
///
// Seqlens: { "chr1" => 30, "chr2" => 100 }
///
/// Sum of all elements: 24.1
///
pub fn granges_test_case_03() -> GRanges<VecRangesIndexed, Vec<Bed5Addition>> {
    create_granges_with_seqlens!(VecRangesIndexed, Vec<Bed5Addition>, {
        "chr1" => [(0, 5, Bed5Addition { name: "a".to_string(), score: Some(1.1) }), 
                   (4, 7, Bed5Addition { name: "b".to_string(), score: Some(8.1) }), 
                   (10, 17, Bed5Addition { name: "c".to_string(), score: Some(10.1) })],
        "chr2" => [(10, 20, Bed5Addition { name: "d".to_string(), score: Some(3.7)}), 
                   (18, 32, Bed5Addition { name: "d".to_string(), score: Some(1.1) })]
    }, seqlens: { "chr1" => 30, "chr2" => 100 })
}

#[cfg(feature = "ndarray")]
/// Mock numeric `NumericSequences1` data
pub fn random_array1_sequences(n: usize) -> GenomeMap<Array1<f64>> {
    let mut seqs = GenomeMap::new();
    // length of sequence
    let max_len = 100_000;
    for seqnum in 1..=n {
        let chrom = format!("chr{}", seqnum);
        let array: Array1<f64> = Array1::from_iter(1..=max_len).map(|&x| x as f64);
        seqs.insert(&chrom, array).unwrap();
    }
    seqs
}

#[cfg(feature = "ndarray")]
/// Mock numeric `NumericSequences2` data
pub fn random_array2_sequences(n: usize) -> GenomeMap<Array2<f64>> {
    let mut seqs = GenomeMap::new();
    // length of sequence
    let max_len = 100_000;
    for seqnum in 1..=n {
        let chrom = format!("chr{}", seqnum);

        // Create a 2-column Array2<f64> with max_len rows
        let rows = max_len as usize;
        let mut array = Array2::<f64>::zeros((rows, 2));

        // Fill the array with values from 0 to 2 * max_len
        let mut value = 0.0;
        for row in array.rows_mut() {
            for elem in row {
                *elem = value;
                value += 1.0;
            }
        }

        seqs.insert(&chrom, array).unwrap();
    }
    seqs
}

#[cfg(test)]
mod tests {
    use crate::{
        granges::GRanges, io::Bed5Iterator, seqlens, test_utilities::granges_test_case_01,
    };

    #[test]
    fn test_test_case_01() {
        let genome = seqlens!("chr1" => 30, "chr2" => 100);
        // yea we're testing test cases. yo dawg.
        let iter = Bed5Iterator::new("tests_data/test_case_01.bed").unwrap();

        // extract out the score, passing the errors through
        // NOTE: this is because we don't have a try_fold iterator
        // method setup yet -- we need that, and when we get it,
        // we should adjust this
        let iter_simplified = iter.map(|result| {
            match result {
                // this unwrap is for '.'
                Ok(record) => Ok(record.into_map_data(|bed5_cols| bed5_cols.score.unwrap())),
                Err(e) => Err(e),
            }
        });

        // build a new GRanges
        let gr = GRanges::from_iter(iter_simplified, &genome).unwrap();
        assert_eq!(gr, granges_test_case_01());
    }
}
