use super::{validate_range, RangeEmpty, RangeIndexed};
use crate::traits::{RangesIntoIterable, RangesIterable};
use crate::{error::GRangesError, traits::RangeContainer, Position};

pub type VecRangesIndexed = VecRanges<RangeIndexed>;
pub type VecRangesEmpty = VecRanges<RangeEmpty>;

#[derive(Clone, Debug)]
pub struct VecRanges<R: Clone> {
    pub(crate) ranges: Vec<R>,
    pub length: Position,
}

impl<R: Clone> VecRanges<R> {
    /// Create a new empty [`VecRanges`] container.
    pub fn new(length: Position) -> Self {
        Self {
            ranges: Vec::new(),
            length,
        }
    }

    /// Validate a range, raising an error if it is invalid for some reason.
    pub fn validate_range(&self, start: Position, end: Position) -> Result<(), GRangesError> {
        validate_range(start, end, self.length)
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

impl RangesIntoIterable<RangeIndexed> for VecRanges<RangeIndexed> {
    fn into_iter_ranges(self) -> Box<dyn Iterator<Item = RangeIndexed>> {
        let iter = self.ranges.into_iter();
        Box::new(iter)
    }
}

impl RangesIterable<RangeIndexed> for VecRanges<RangeIndexed> {
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = RangeIndexed> + '_> {
        let iter = self.ranges.iter();
        // NOTE: RangeIndexed is copyable but there still is overhead here
        let converted_iter = iter.map(|interval| interval.to_owned());
        Box::new(converted_iter)
    }
}
