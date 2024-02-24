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

/// [`JoinData`] contains a [`Vec<LeftGroupedJoin>`] of all overlap
/// joins, as well as references to the left and right data containers.
#[derive(Clone, Debug)]
pub struct JoinData<'a, DL, DR> {
    pub joins: Vec<LeftGroupedJoin>,
    pub left_data: Option<DL>,
    pub right_data: Option<&'a DR>,
}

impl<'a, DL, DR> JoinData<'a, DL, DR> {
    /// Create a new [`JoinData`].
    pub fn new(left_data: Option<DL>, right_data: Option<&'a DR>) -> Self {
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

    /// Create an iterator over the joins.
    pub fn iter(&'a self) -> JoinDataIterator<'a, DL, DR> {
        JoinDataIterator {
            inner: self.joins.iter(),
            left_data: self.left_data.as_ref(),
            right_data: self.right_data,
        }
    }
}

impl<'a, DL, DR> JoinData<'a, DL, DR>
where
    DL: IndexedDataContainer + 'a,
    DR: IndexedDataContainer + 'a,
{
    pub fn apply_into_vec<F, V>(&self, func: F) -> Vec<V>
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
                let left_data = join.left.and_then(|idx| {
                    self.left_data
                        .as_ref()
                        .and_then(|data| Some(data.get_owned(idx)))
                });

                let right_data = join.rights.as_ref().and_then(|indices| {
                    Some(
                        indices
                            .iter()
                            .filter_map(|&idx| {
                                self.right_data
                                    .as_ref()
                                    .and_then(|data| Some(data.get_owned(idx)))
                            })
                            .collect::<Vec<_>>(),
                    )
                });

                // Now `func` is applied to each `CombinedJoinData`
                func(CombinedJoinData {
                    left_data,
                    right_data,
                })
            })
            .collect()
    }
}

/// Represents a combined view of a single join operation along with references to
/// associated left and right data.
///
/// This struct is particularly useful for iterating over join results while maintaining
/// access to the original data elements that were involved in each join. It encapsulates
/// a reference to the join information (`join`), which details how two data elements are
/// related (e.g., through overlap or proximity). Additionally, it holds optional references
/// to the data elements themselves (`left_data` and `right_data`), allowing for easy retrieval
/// and inspection of the data involved in the join.
pub struct CombinedJoinData<DL, DR> {
    // pub join: LeftGroupedJoin,      // The join information
    pub left_data: Option<DL>,       // The left data element
    pub right_data: Option<Vec<DR>>, // The right data elements
}

/// An iterator over the [`LeftGroupedJoin`] types that represent
/// information about overlaps right ranges have with a particular left range.
///
/// This also contains references to the left and right data containers, for
/// better ergonomics in downstream data processing.
pub struct JoinDataIterator<'a, DL, DR> {
    inner: std::slice::Iter<'a, LeftGroupedJoin>,
    pub left_data: Option<&'a DL>,
    pub right_data: Option<&'a DR>,
}

impl<'a, DL, DR> Iterator for JoinDataIterator<'a, DL, DR> {
    type Item = &'a LeftGroupedJoin;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

#[cfg(test)]
mod tests {
    use crate::ranges::RangeIndexed;

    use super::{JoinData, LeftGroupedJoin};

    #[test]
    fn test_join_data_new() {
        let left_data = vec![1, 2];
        let right_data = vec![4, 8];
        let mut jd = JoinData::new(Some(&left_data), Some(&right_data));
        assert_eq!(jd.len(), 0);

        let left = RangeIndexed::new(0, 10, 1);
        let mut join = LeftGroupedJoin::new(&left);
        let right = RangeIndexed::new(8, 10, 1);
        join.add_right(&left, &right);
        jd.push(join);
        assert_eq!(jd.len(), 1);
    }

    #[test]
    fn test_join_iter() {
        let left_data = vec![1, 2];
        let right_data = vec![4, 8];

        let mut jd = JoinData::new(Some(&left_data), Some(&right_data));

        let left = RangeIndexed::new(0, 10, 1);
        let mut join = LeftGroupedJoin::new(&left);
        let right = RangeIndexed::new(8, 10, 1);
        join.add_right(&left, &right);
        jd.push(join);

        let right = RangeIndexed::new(9, 11, 1);
        let mut join = LeftGroupedJoin::new(&left);
        join.add_right(&left, &right);
        jd.push(join);

        let mut iter = jd.iter();
        assert_eq!(iter.next().unwrap().num_overlaps(), 1);
        assert_eq!(iter.next().unwrap().num_overlaps(), 1);
        assert_eq!(iter.next(), None);
    }
}
