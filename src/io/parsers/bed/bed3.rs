//! BED3 Parsers, which are built off of the [`GenomicRangeRecordEmpty`].
//!

use crate::{io::TsvRecordIterator, ranges::GenomicRangeRecordEmpty, GRangesError};
use std::path::PathBuf;

/// An iterator over BED3 entries (which just contain ranges no data).
#[derive(Debug)]
pub struct Bed3Iterator {
    iter: TsvRecordIterator<GenomicRangeRecordEmpty>,
}

impl Bed3Iterator {
    /// Creates a parsing iterator over a BED5 file.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let iter = TsvRecordIterator::new(filepath)?;
        Ok(Self { iter })
    }
}

impl Iterator for Bed3Iterator {
    type Item = Result<GenomicRangeRecordEmpty, GRangesError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}
