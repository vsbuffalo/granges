//! BED4 Parsers, which are built off of the [`GenomicRangeRecordEmpty`]
//! and [`Bed4Addition`].

use crate::{io::TsvRecordIterator, ranges::GenomicRangeRecord, GRangesError};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// The additional two BED5 columns.
///
/// # Fields
/// * `name`: the feature name.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Bed4Addition {
    pub name: String,
}

/// An iterator over BED4 entries, which contain the three
/// range entries (sequence name, start and end positions),
/// and a feature name
///
/// Note that the [`Bed5Addition`] is *permissive* in the sense
/// it allows for more input types than the official [BED format
/// specification](https://samtools.github.io/hts-specs/BEDv1.pdf).
// TODO strict type?
#[derive(Debug)]
pub struct Bed4Iterator {
    iter: TsvRecordIterator<GenomicRangeRecord<Bed4Addition>>,
}

impl Bed4Iterator {
    /// Creates a parsing iterator over a BED4 file.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let iter = TsvRecordIterator::new(filepath)?;

        Ok(Self { iter })
    }
}

impl Iterator for Bed4Iterator {
    type Item = Result<GenomicRangeRecord<Bed4Addition>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}
