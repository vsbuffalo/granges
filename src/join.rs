//! [`LeftGroupedJoin`], [`JoinData`], and [`JoinDataIterator`] types for overlaps.
//!
#![allow(clippy::all)]

use crate::{
    traits::{GenericRange, IndexedDataContainer},
    Position,
};

/// [`LeftGroupedJoin`] contains information about the right ranges
/// and their degree of overlap with a focal left range. This information
/// is designed to facilitate downstream statistical sumamries of the
/// corresponding data in overlapping ranges.
#[derive(Clone, Debug, PartialEq)]
pub struct LeftGroupedJoin {
    /// The data index for the left range.
    left: Option<usize>,

    /// A `Vec` of the indices for the overlapping right ranges.
    /// This is `None` if the right ranges do not have a data container.
    rights: Option<Vec<usize>>,

    /// The length of the left range.
    left_length: Position,

    /// The lengths of the right ranges.
    right_lengths: Vec<Position>,

    /// The lengths of the overlaps between the left and right ranges.
    overlaps: Vec<Position>,
    // TODO: we may want some simple summary of whether an overlapping range is
    // up or downstream. I think the cleanest summary is a signed integer
    // representing the side and degree of non-overlap. E.g. a range
    // that overlaps another but overhangs the 3' side of the focal left
    // range by 10bp is +10; if it were 5', it would be -10.

    // left_data: Option<&'a DL>,
    // right_data: Option<&'a DR>,
}

impl LeftGroupedJoin {
    /// Create a new [`LeftGroupedJoin`].
    pub fn new<R: GenericRange>(left_range: &R) -> Self {
        Self {
            left: left_range.index(),
            rights: None,
            left_length: left_range.width(),
            right_lengths: Vec::new(),
            overlaps: Vec::new(),
            // left_data,
            // right_data,
        }
    }
    /// Add a right (overlapping) range to this [`LeftGroupedJoin`].
    ///
    // Note: in principle, this can be called on *non-overlapping* right ranges too,
    // for a full-outer join.
    pub fn add_right<R: GenericRange, Q: GenericRange>(&mut self, left: &R, right: &Q) {
        if let Some(right_index) = right.index() {
            // the right range has data -- add to vec, initializing if not there
            self.rights.get_or_insert_with(Vec::new).push(right_index)
        }
        self.right_lengths.push(right.width());
        self.overlaps.push(left.overlap_width(right));
    }
    /// Return whether this left range has any [`LeftGroupedJoin`].
    pub fn has_overlaps(&self) -> bool {
        !self.overlaps.is_empty()
    }

    /// Retrieve the number of right overlaps.
    pub fn num_overlaps(&self) -> usize {
        self.overlaps.len()
    }
}

/// [`JoinData`] contains a [`Vec<LeftGroupedJoin>`] of all overlap joins,
/// and owns the left data container from the join. It stores a reference
/// to the right data container.
#[derive(Clone, Debug)]
pub struct JoinData<'a, DL, DR> {
    pub joins: Vec<LeftGroupedJoin>,
    pub left_data: DL,
    pub right_data: &'a DR,
}

impl<'a, DL, DR> JoinData<'a, DL, DR> {
    /// Create a new [`JoinData`].
    pub fn new(left_data: DL, right_data: &'a DR) -> Self {
        let joins = Vec::new();
        JoinData {
            joins,
            left_data,
            right_data,
        }
    }

    /// Push the [`LeftGroupedJoin`] to joins.
    pub fn push(&mut self, join: LeftGroupedJoin) {
        self.joins.push(join)
    }

    /// Get the total number of joins.
    pub fn len(&self) -> usize {
        self.joins.len()
    }

    /// Return whether the [`JoinData`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    ///// Create an iterator over the joins.
    //pub fn iter(&'a self) -> JoinDataIterator<'a, DL, DR> {
    //    JoinDataIterator {
    //        inner: self.joins.iter(),
    //        left_data: self.left_data.as_ref(),
    //        right_data: self.right_data,
    //    }
    //}
}

pub struct CombinedJoinData<DL, DR> {
    pub join: LeftGroupedJoin, // Information on the join
    pub left_data: DL,         // The left data element
    pub right_data: Vec<DR>,   // The right data elements
}

impl<'a, DL, DR> JoinData<'a, DL, DR>
where
    DL: IndexedDataContainer + 'a,
    DR: IndexedDataContainer + 'a,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn apply<F, V>(&self, func: F) -> Vec<V>
    where
        F: Fn(
            CombinedJoinData<
                <DL as IndexedDataContainer>::OwnedItem,
                <DR as IndexedDataContainer>::OwnedItem,
            >,
        ) -> V,
    {
        // Cloning `left_data` and `right_data` to ensure they live long enough.
        // This might not be the most efficient but ensures lifetime correctness.

        self.joins
            .iter()
            .map(|join| {
                let left_data = self.left_data.get_owned(join.left.unwrap());
                let right_data = join
                    .rights.as_ref()
                    .unwrap()
                    .iter()
                    .map(|idx| self.right_data.get_owned(*idx))
                    .collect();
                // Now `func` is applied to each `CombinedJoinData`
                func(CombinedJoinData {
                    join: join.clone(),
                    left_data,
                    right_data,
                })
            })
            .collect()
    }
}

/// [`JoinDataLeftEmpty`] contains a [`Vec<LeftGroupedJoin>`] of all overlap joins,
/// and stores a reference to the right data container.
#[derive(Clone, Debug)]
pub struct JoinDataLeftEmpty<'a, DR> {
    pub joins: Vec<LeftGroupedJoin>,
    pub right_data: &'a DR,
}

impl<'a, DR> JoinDataLeftEmpty<'a, DR> {
    /// Create a new [`JoinData`].
    pub fn new(right_data: &'a DR) -> Self {
        let joins = Vec::new();
        JoinDataLeftEmpty { joins, right_data }
    }

    /// Push the [`LeftGroupedJoin`] to joins.
    pub fn push(&mut self, join: LeftGroupedJoin) {
        self.joins.push(join)
    }

    /// Get the total number of joins.
    pub fn len(&self) -> usize {
        self.joins.len()
    }

    /// Return whether the [`JoinData`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct CombinedJoinDataLeftEmpty<DR> {
    pub join: LeftGroupedJoin, // Information on the join
    pub right_data: Vec<DR>,   // The right data element
}

impl<'a, DR> JoinDataLeftEmpty<'a, DR>
where
    DR: IndexedDataContainer + 'a,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn apply<F, V>(&self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataLeftEmpty<<DR as IndexedDataContainer>::OwnedItem>) -> V,
    {
        // Cloning `left_data` and `right_data` to ensure they live long enough.
        // This might not be the most efficient but ensures lifetime correctness.

        self.joins
            .iter()
            .map(|join| {
                let right_data = join
                    .rights.as_ref()
                    .unwrap()
                    .iter()
                    .map(|idx| self.right_data.get_owned(*idx))
                    .collect();
                // Now `func` is applied to each `CombinedJoinDataLeftEmpty`
                func(CombinedJoinDataLeftEmpty {
                    join: join.clone(),
                    right_data,
                })
            })
            .collect()
    }
}

/// [`JoinDataRightEmpty`] contains a [`Vec<LeftGroupedJoin>`] of all overlap joins,
/// and owns the left data.
#[derive(Clone, Debug)]
pub struct JoinDataRightEmpty<DR> {
    pub joins: Vec<LeftGroupedJoin>,
    pub left_data: DR,
}

impl<'a, DL> JoinDataRightEmpty<DL> {
    /// Create a new [`JoinData`].
    pub fn new(left_data: DL) -> Self {
        let joins = Vec::new();
        JoinDataRightEmpty { joins, left_data }
    }

    /// Push the [`LeftGroupedJoin`] to joins.
    pub fn push(&mut self, join: LeftGroupedJoin) {
        self.joins.push(join)
    }

    /// Get the total number of joins.
    pub fn len(&self) -> usize {
        self.joins.len()
    }

    /// Return whether the [`JoinData`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct CombinedJoinDataRightEmpty<DL> {
    pub join: LeftGroupedJoin, // Information on the join
    pub left_data: DL,         // The right data element
}

impl<'a, DL> JoinDataRightEmpty<DL>
where
    DL: IndexedDataContainer,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn apply<F, V>(&self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataRightEmpty<<DL as IndexedDataContainer>::OwnedItem>) -> V,
    {
        // Cloning `left_data` and `right_data` to ensure they live long enough.
        // This might not be the most efficient but ensures lifetime correctness.

        self.joins
            .iter()
            .map(|join| {
                let left_data = self.left_data.get_owned(join.left.unwrap());
                // Now `func` is applied to each `CombinedJoinDataRightEmpty`
                func(CombinedJoinDataRightEmpty {
                    join: join.clone(),
                    left_data,
                })
            })
            .collect()
    }
}

/// [`JoinDataBothEmpty`] contains a [`Vec<LeftGroupedJoin>`] of all overlap joins
/// without any owned or references to data containers.
#[derive(Clone, Debug)]
pub struct JoinDataBothEmpty {
    pub joins: Vec<LeftGroupedJoin>,
}

impl JoinDataBothEmpty {
    /// Create a new [`JoinData`].
    pub fn new() -> Self {
        let joins = Vec::new();
        JoinDataBothEmpty { joins }
    }

    /// Push the [`LeftGroupedJoin`] to joins.
    pub fn push(&mut self, join: LeftGroupedJoin) {
        self.joins.push(join)
    }

    /// Get the total number of joins.
    pub fn len(&self) -> usize {
        self.joins.len()
    }

    /// Return whether the [`JoinData`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct CombinedJoinDataBothEmpty {
    pub join: LeftGroupedJoin,
}

impl JoinDataBothEmpty {
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn apply<F, V>(&self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataBothEmpty) -> V,
    {
        // Cloning `left_data` and `right_data` to ensure they live long enough.
        // This might not be the most efficient but ensures lifetime correctness.

        self.joins
            .iter()
            .map(|join| func(CombinedJoinDataBothEmpty { join: join.clone() }))
            .collect()
    }
}

///// An iterator over the [`LeftGroupedJoin`] types that represent
///// information about overlaps right ranges have with a particular left range.
/////
///// This also contains references to the left and right data containers, for
///// better ergonomics in downstream data processing.
//pub struct JoinDataIterator<'a, DL, DR> {
//    inner: std::slice::Iter<'a, LeftGroupedJoin>,
//    pub left_data: Option<&'a DL>,
//    pub right_data: Option<&'a DR>,
//}
//
//impl<'a, DL, DR> Iterator for JoinDataIterator<'a, DL, DR> {
//    type Item = &'a LeftGroupedJoin;
//
//    fn next(&mut self) -> Option<Self::Item> {
//        self.inner.next()
//    }
//}

#[cfg(test)]
mod tests {
    use crate::ranges::RangeIndexed;

    use super::{JoinData, LeftGroupedJoin};

    #[test]
    fn test_join_data_new() {
        let left_data = vec![1, 2];
        let right_data = vec![4, 8];
        let mut jd = JoinData::new(left_data, &right_data);
        assert_eq!(jd.len(), 0);

        let left = RangeIndexed::new(0, 10, 1);
        let mut join = LeftGroupedJoin::new(&left);
        let right = RangeIndexed::new(8, 10, 1);
        join.add_right(&left, &right);
        jd.push(join);
        assert_eq!(jd.len(), 1);
    }
}
