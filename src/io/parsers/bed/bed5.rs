//! BED5 Parsers, which are built off of the [`GenomicRangeRecordEmpty`]
//! and [`Bed5Addition`].

use super::bed_missing;
use crate::{io::TsvRecordIterator, ranges::GenomicRangeRecord, GRangesError};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// The additional two BED5 columns.
///
/// # Fields
/// * `name`: the feature name.
/// * `score`: a score.
// TODO/RENAME: maybe since this goes against spec it should
// be called Bed5AdditionPermissive?
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Bed5Addition {
    pub name: String,
    #[serde(deserialize_with = "bed_missing")]
    pub score: Option<f64>,
}

/// An iterator over BED5 entries, which contain the three
/// range entries (sequence name, start and end positions),
/// a feature name, and a score.
///
/// Note that the [`Bed5Addition`] is *permissive* in the sense
/// it allows for more input types than the official [BED format
/// specification](https://samtools.github.io/hts-specs/BEDv1.pdf).
// TODO strict type?
#[derive(Debug)]
pub struct Bed5Iterator {
    iter: TsvRecordIterator<GenomicRangeRecord<Bed5Addition>>,
}

impl Bed5Iterator {
    /// Creates a parsing iterator over a BED5 file.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let iter = TsvRecordIterator::new(filepath)?;

        Ok(Self { iter })
    }
}

impl Iterator for Bed5Iterator {
    type Item = Result<GenomicRangeRecord<Bed5Addition>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}
