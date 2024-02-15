use coitrees::{Interval, BasicCOITree, IntervalTree, IntervalNode, GenericInterval};

use crate::{Position, traits::RangeContainer, error::GRangesError};

use super::{vec::VecRanges, RangeIndexed, validate_range};

type COITreeIntervalIndexed = Interval<usize>;

impl GenericInterval<usize> for RangeIndexed {
    fn first(&self) -> i32 {
        self.start().try_into().unwrap()
    }
    fn last(&self) -> i32 {
        self.end().try_into().unwrap()
    }
    fn metadata(&self) -> &usize {
       self.index()
    }
}

/// A [`coitrees::BasicCOITree`] interval tree for a single sequence's ranges.
///
/// This is generic over the interval type, to handle the case where one 
/// may want to do overlap operations on ranges without associated data in 
/// a data container (e.g. ranges that just define megabase windwows).
pub struct COITreeRangeContainer<R: Clone> {
    ranges: BasicCOITree<R, usize>,
    /// The sequence length, used to validate new ranges.
    length: Position,
}

impl<R: Clone> COITreeRangeContainer<R> {
    pub fn validate_range(&self, start: Position, end: Position) -> Result<(), GRangesError> {
        let range = start..end;
        validate_range(&range, self.length)
    }

    pub fn query<F>(&self, start: Position, end: Position, visit: F) 
    where F: FnMut(&IntervalNode<R, usize>) {
        // Note the terminology change to match coitrees (and uses i32s)
        let first = start.try_into().expect("could not covert");
        let end: i32  = end.try_into().expect("could not covert");
        // internally coitrees uses 0-indexed, right-inclusive "last"
        self.ranges.query(first, end - 1, visit)
    }

    /// Return the number of ranges in this [`COITreeRangeContainer`] container.
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    /// Return whether the [`COITreeRangeContainer`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

}

impl<R: Clone + GenericInterval<R>> From<VecRanges<R>> for COITreeRangeContainer<R> {
    fn from(value: VecRanges<R>) -> Self {
        let ranges = BasicCOITree::new(&value.ranges);
        let length = value.length;
        Self {
            ranges,
            length
        }
    }
}

impl<R: Clone> RangeContainer for COITreeRangeContainer<R> {
    fn len(&self) -> usize {
        self.ranges.len()
    }
}
