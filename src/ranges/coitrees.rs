//! The [`COITrees<M>`] type.
//!
//! This wraps the functionality og the [`coitrees`] library by Daniel C. Jones.
//!
use coitrees::{BasicCOITree, GenericInterval, Interval, IntervalNode, IntervalTree};

use crate::{
    error::GRangesError,
    traits::IterableRangeContainer,
    traits::{GenericRange, RangeContainer},
    Position,
};

use super::{validate_range, vec::VecRanges, RangeEmpty, RangeIndexed};

pub type COITreesIndexed = COITrees<usize>;
pub type COITreesEmpty = COITrees<()>;

/// A [`coitrees::BasicCOITree`] interval tree for a single sequence's ranges.
///
/// This is generic over the interval type, to handle the case where one
/// may want to do overlap operations on ranges without associated data in
/// a data container (e.g. ranges that just define megabase windwows).
///
pub struct COITrees<M: Clone> {
    pub(crate) ranges: BasicCOITree<M, usize>,
    /// The sequence length, used to validate new ranges.
    pub length: Position,
}

impl<M: Clone> std::fmt::Debug for COITrees<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("COITrees")
            .field("number of ranges:", &self.ranges.len())
            .field("length", &self.length)
            .finish_non_exhaustive()
    }
}

impl<M: Clone> COITrees<M> {
    /// Validate a range, raising an error if it is invalid for some reason.
    pub fn validate_range(&self, start: Position, end: Position) -> Result<(), GRangesError> {
        validate_range(start, end, self.length)
    }

    /// Query this range container for a particular range, and call a visit function on all
    /// overlapping ranges.
    pub fn query<F>(&self, start: Position, end: Position, visit: F)
    where
        F: FnMut(&IntervalNode<M, usize>),
    {
        // Note the terminology change to match coitrees (and uses i32s)
        let first = start.try_into().expect("could not covert");
        let end: i32 = end.try_into().expect("could not covert");
        // internally coitrees uses 0-indexed, right-inclusive "last"
        self.ranges.query(first, end - 1, visit)
    }

    /// Returns the number of ranges that overlap the specified range.
    pub fn count_overlaps(&self, start: Position, end: Position) -> usize {
        let first = start.try_into().expect("could not covert");
        let end: i32 = end.try_into().expect("could not covert");
        self.ranges.query_count(first, end - 1)
    }

    /// Return the number of ranges in this [`COITrees`] container.
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    /// Return whether the [`COITrees`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Convert a [`VecRanges`] range container to a [`COITrees`] range container.
impl<R: Clone + GenericInterval<M>, M: Clone> From<VecRanges<R>> for COITrees<M> {
    fn from(value: VecRanges<R>) -> Self {
        let ranges = BasicCOITree::new(&value.ranges);
        let length = value.length;
        Self { ranges, length }
    }
}

/// Convert a [`COITrees`] range container to a [`VecRanges`] range container.
impl From<COITrees<usize>> for VecRanges<RangeIndexed> {
    fn from(value: COITrees<usize>) -> Self {
        let length = value.length;
        let mut ranges: VecRanges<RangeIndexed> = VecRanges::new(length);
        for interval in value.ranges.iter() {
            ranges.push_range(interval.into());
        }
        ranges
    }
}

/// [`RangeContainer`] trait implementations.
impl<R: Clone> RangeContainer for COITrees<R> {
    fn len(&self) -> usize {
        self.ranges.len()
    }
    fn sequence_length(&self) -> Position {
        self.length
    }
}

/// Convert between [`coitrees::Interval`] with index metadata to a [`RangeEmpty`].
impl From<Interval<&()>> for RangeEmpty {
    fn from(value: Interval<&()>) -> Self {
        let first: Position = value.first.try_into().unwrap();
        let last: Position = value.last.try_into().unwrap();
        RangeEmpty {
            start: first,
            end: last + 1,
        }
    }
}

/// Convert between [`coitrees::Interval`] with index metadata to a [`RangeIndexed`].
impl From<Interval<&usize>> for RangeIndexed {
    fn from(value: Interval<&usize>) -> Self {
        let first: Position = value.first.try_into().unwrap();
        let last: Position = value.last.try_into().unwrap();
        RangeIndexed {
            start: first,
            end: last + 1,
            index: *value.metadata(),
        }
    }
}

/// # Developer Notes
///
/// Internally, the [`coitrees`] iterator is over their interval type.
/// Their iterator does not consume, but we need an owned (or copyable,
/// as is case here) type to map out. Thus, we need this odd variant of
/// an iterator that doesn't return references and does not consume.
impl IterableRangeContainer for COITrees<usize> {
    type RangeType = RangeIndexed;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = RangeIndexed> + '_> {
        let iter = self.ranges.iter();
        let converted_iter = iter.map(RangeIndexed::from);
        Box::new(converted_iter)
    }
}

impl IterableRangeContainer for COITrees<()> {
    type RangeType = RangeEmpty;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = RangeEmpty> + '_> {
        let iter = self.ranges.iter();
        let converted_iter = iter.map(RangeEmpty::from);
        Box::new(converted_iter)
    }
}

impl GenericInterval<()> for RangeEmpty {
    fn first(&self) -> i32 {
        self.start.try_into().unwrap()
    }
    fn last(&self) -> i32 {
        let end: i32 = self.end.try_into().unwrap();
        // convert right-exclusive to inclusive
        end - 1
    }
    fn metadata(&self) -> &() {
        &()
    }
}

impl GenericInterval<usize> for RangeIndexed {
    fn first(&self) -> i32 {
        self.start.try_into().unwrap()
    }
    fn last(&self) -> i32 {
        let end: i32 = self.end.try_into().unwrap();
        // convert right-exclusive to inclusive
        end - 1
    }
    fn metadata(&self) -> &usize {
        &self.index
    }
}

// Note: this is predominantly used in the join logic, where overlapping
// ranges as [`coitrees::IntervalNode`] are processed.
impl GenericRange for IntervalNode<(), usize> {
    fn start(&self) -> Position {
        self.first().try_into().unwrap()
    }
    fn end(&self) -> Position {
        (self.last() + 1).try_into().unwrap()
    }
    fn index(&self) -> Option<usize> {
        None
    }
}

// Note: this is predominantly used in the join logic, where overlapping
// ranges as [`coitrees::IntervalNode`] are processed.
impl GenericRange for IntervalNode<usize, usize> {
    fn start(&self) -> Position {
        self.first().try_into().unwrap()
    }
    fn end(&self) -> Position {
        (self.last() + 1).try_into().unwrap()
    }
    fn index(&self) -> Option<usize> {
        Some(*self.metadata())
    }
}

impl<M: Clone + PartialEq> PartialEq for COITrees<M> {
    fn eq(&self, other: &Self) -> bool {
        if self.ranges.len() != other.ranges.len() {
            return false;
        }
        let self_ranges = self.ranges.iter();
        let other_ranges = other.ranges.iter();
        self_ranges.zip(other_ranges).all(|(r1, r2)| {
            // NOTE: we really should just implement PartialEq here in coitrees
            (r1.first == r2.first) && (r1.last == r2.last) && (r1.metadata == r2.metadata)
        })
    }
}

#[cfg(test)]
mod tests {
    use coitrees::{GenericInterval, Interval};

    use crate::prelude::*;
    use crate::ranges::{RangeEmpty, RangeIndexed};
    use crate::test_utilities::granges_test_case_01;

    #[test]
    fn test_ranges_iterable_coitrees() {
        let gr = granges_test_case_01().into_coitrees().unwrap();
        let mut chr1_iter = gr.get_ranges("chr1").unwrap().iter_ranges();
        assert_eq!(chr1_iter.next().unwrap(), RangeIndexed::new(0, 5, 0));
        assert_eq!(chr1_iter.next().unwrap(), RangeIndexed::new(4, 7, 1));
        assert_eq!(chr1_iter.next().unwrap(), RangeIndexed::new(10, 17, 2));
        assert!(chr1_iter.next().is_none());
    }

    #[test]
    fn test_from_interval_to_range_empty() {
        let interval: Interval<&()> = Interval::new(0, 10, &());
        let range_empty: RangeEmpty = interval.into();
        assert_eq!(range_empty.start, 0);
        assert_eq!(range_empty.end, 11);
        assert_eq!(range_empty.start(), 0);
        assert_eq!(range_empty.end(), 11);
    }

    #[test]
    fn test_from_interval_to_range_indexed() {
        let interval: Interval<&usize> = Interval::new(0, 10, &0);
        let range_empty: RangeIndexed = interval.into();
        assert_eq!(range_empty.start, 0);
        assert_eq!(range_empty.end, 11);
        assert_eq!(range_empty.start(), 0);
        assert_eq!(range_empty.end(), 11);
    }

    #[test]
    fn test_generic_interval_for_range_empty() {
        let range = RangeEmpty::new(0, 10);
        assert_eq!(range.first(), 0);
        assert_eq!(range.last(), 9);
    }

    #[test]
    fn test_generic_interval_for_range_indexed() {
        let range = RangeIndexed::new(0, 10, 1);
        assert_eq!(range.first(), 0);
        assert_eq!(range.last(), 9);
    }
}
