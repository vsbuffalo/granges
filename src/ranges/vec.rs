use super::operations::adjust_range;
use super::{validate_range, RangeEmpty, RangeIndexed};
use crate::traits::{GenericRange, IntoIterableRangesContainer, IterableRangeContainer};
use crate::PositionOffset;
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

impl<R: GenericRange> VecRanges<R> {
    /// Sort all the ranges.
    pub fn sort(&mut self) {
        self.ranges.sort_by(|a, b| {
            a.start()
                .cmp(&b.start())
                .then_with(|| a.end().cmp(&b.end()))
                .then_with(|| a.index().cmp(&b.index()))
        });
    }

    /// Adjust all the ranges in this [`VecRanges`] range container.
    pub fn adjust_ranges(&mut self, start_delta: PositionOffset, end_delta: PositionOffset) {
        let mut ranges = std::mem::replace(&mut self.ranges, Vec::new());

        ranges = ranges
            .into_iter()
            .filter_map(|range| adjust_range(range, start_delta, end_delta, self.length))
            .collect();

        self.ranges = ranges;
    }
}

impl<R: Clone> RangeContainer for VecRanges<R> {
    fn len(&self) -> usize {
        self.ranges.len()
    }
    fn sequence_length(&self) -> Position {
        self.length
    }
}

impl IntoIterableRangesContainer<RangeIndexed> for VecRanges<RangeIndexed> {
    fn into_iter_ranges(self) -> Box<dyn Iterator<Item = RangeIndexed>> {
        let iter = self.ranges.into_iter();
        Box::new(iter)
    }
}

impl IterableRangeContainer for VecRanges<RangeIndexed> {
    type RangeType = RangeIndexed;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = RangeIndexed> + '_> {
        let iter = self.ranges.iter();
        // NOTE: RangeIndexed is copyable but there still is overhead here
        let converted_iter = iter.map(|interval| interval.to_owned());
        Box::new(converted_iter)
    }
}

impl IntoIterableRangesContainer<RangeEmpty> for VecRanges<RangeEmpty> {
    fn into_iter_ranges(self) -> Box<dyn Iterator<Item = RangeEmpty>> {
        let iter = self.ranges.into_iter();
        Box::new(iter)
    }
}

impl IterableRangeContainer for VecRanges<RangeEmpty> {
    type RangeType = RangeEmpty;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = RangeEmpty> + '_> {
        let iter = self.ranges.iter();
        // NOTE: RangeIndexed is copyable but there still is overhead here
        let converted_iter = iter.map(|interval| interval.to_owned());
        Box::new(converted_iter)
    }
}
