use crate::{traits::RangeContainer, Position, error::GRangesError};

use super::{RangeIndexed, validate_range, RangeEmpty};
pub type VecRangesIndexed = VecRanges<RangeIndexed>;
pub type VecRangesEmpty = VecRanges<RangeEmpty>;

#[derive(Clone, Default)]
pub struct VecRanges<R: Clone> {
    pub (crate) ranges: Vec<R>,
    pub length: Position,
}

impl<R: Clone> VecRanges<R> {
    pub fn validate_range(&self, start: Position, end: Position) -> Result<(), GRangesError> {
        let range = start..end;
        validate_range(&range, self.length)
    }

    /// Create a new empty [`VecRanges`] container.
    pub fn new(length: Position) -> Self {
        Self {
            ranges: Vec::new(),
            length,
        }
    }

    /// Add a new range to the [`VecRanges`] container.
    pub fn push_range(&mut self, range: R) {
        self.ranges.push(range)
    }

    /// Return the number of ranges in this [`VecRanges`] container.
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    /// Return whether the [`VecRanges`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<R: Clone> RangeContainer for VecRanges<R> {
    fn len(&self) -> usize {
        self.ranges.len()
    }
}
