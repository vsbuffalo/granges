//! [`LeftGroupedJoin`], [`JoinData`], and [`JoinDataIterator`] types for overlaps.
//!
#![allow(clippy::all)]

use std::collections::HashSet;

use crate::{
    traits::{GenericRange, IndexedDataContainer, JoinDataOperations},
    Position,
};

/// This is a generic range used just in join logic, to avoid
/// having to handle bringing range types into [`LeftGroupedJoin`],
/// which clog up the API a bit.
/// These can represent indexed and empty ranges.
#[derive(Clone, Debug, PartialEq)]
pub struct RangeTuple((Position, Position, Option<usize>));

impl GenericRange for RangeTuple {
    fn start(&self) -> Position {
        self.0 .0
    }
    fn end(&self) -> Position {
        self.0 .1
    }
    fn index(&self) -> Option<usize> {
        self.0 .2
    }
}

/// This is a special "reduced" range that stores indices
/// to multiple data elements.
#[derive(Clone, Debug, PartialEq)]
pub struct RangeReduced((Position, Position, HashSet<Option<usize>>));

impl RangeReduced {
    pub fn indices(&self) -> &HashSet<Option<usize>> {
        &self.0 .2
    }
}

/// Create a vector of "reduced" or flattened ranges, that stack the indices of each range.
///
///
/// This first finds every range position, sorts, and dedups these. Then,
/// it iterates through each range of adjacent positions. Each of these
/// new ranges is guaranteed to be covered by >= 1 range.
pub fn reduce_ranges<R: GenericRange>(ranges: &Vec<R>) -> Vec<RangeReduced> {
    let mut range_ends: Vec<Position> = ranges
        .iter()
        .flat_map(|range| vec![range.start(), range.end()].into_iter())
        .collect();
    range_ends.sort_unstable();
    range_ends.dedup();

    let mut ranges_reduced = Vec::new();
    for range in range_ends.windows(2) {
        let mut indices = HashSet::new();
        if let [start, end] = range {
            for range in ranges {
                if range.start() < *end && range.end() > *start {
                    indices.insert(range.index());
                }
            }
            if !indices.is_empty() {
                ranges_reduced.push(RangeReduced((*start, *end, indices)));
            }
        }
    }
    ranges_reduced
}

impl GenericRange for RangeReduced {
    fn start(&self) -> Position {
        self.0 .0
    }
    fn end(&self) -> Position {
        self.0 .1
    }
    // Note: [`RangeReduced`] do not have valid indices,
    // so this returns [`None`]. (They do have a [`Vec<Option<usize>>`]
    // of indices, but this has to be accessed through the type.)
    fn index(&self) -> Option<usize> {
        None
    }
}

/// [`LeftGroupedJoin`] contains information about the right ranges
/// and their degree of overlap with a focal left range. This information
/// is designed to facilitate downstream statistical sumamries of the
/// corresponding data in overlapping ranges.
#[derive(Clone, Debug, PartialEq)]
pub struct LeftGroupedJoin {
    /// The left range.
    pub left: RangeTuple,

    /// A `Vec` of the right overlapping ranges (unsorted).
    // NOTE: previously lengths, overlap width, and their data
    // indices were stored. For just one extra u32 we can store
    // all the data in the original structure.
    pub rights: Vec<RangeTuple>,
}

impl LeftGroupedJoin {
    /// Create a new [`LeftGroupedJoin`].
    pub fn new<R: GenericRange>(left_range: &R) -> Self {
        Self {
            left: RangeTuple(left_range.as_tuple()),
            rights: Vec::new(),
        }
    }
    /// Add a right (overlapping) range to this [`LeftGroupedJoin`].
    ///
    // Note: in principle, this can be called on *non-overlapping* right ranges too,
    // for a full-outer join.
    pub fn add_right<R: GenericRange>(&mut self, right: &R) {
        self.rights.push(RangeTuple(right.as_tuple()))
    }

    /// Sort the right ranges, for faster downstream processing.
    pub fn sort_ranges(&mut self) {
        self.rights.sort_by(|a, b| {
            a.start()
                .cmp(&b.start())
                .then_with(|| a.end().cmp(&b.end()))
                .then_with(|| a.index().cmp(&b.index()))
        });
    }

    /// "Reduce" the ranges into a minimum set, with all
    /// their indices gathered in a [`Vec<Option<usize>>`].
    /// This returns a [`Vec<RangeReduced>`].
    pub fn reduce_ranges(&mut self) -> Vec<RangeReduced> {
        // we need to trim these by the left range to get the
        // proper overlaps within this left range
        let rights: Vec<_> = self
            .rights
            .iter()
            .map(|range| {
                let (start, end) = range.overlap_range(&self.left).unwrap();
                RangeTuple((start, end, range.index()))
            })
            .collect();
        reduce_ranges(&rights)
    }

    /// Return whether this left range has any [`LeftGroupedJoin`].
    pub fn has_overlaps(&self) -> bool {
        !self.overlaps().is_empty()
    }

    /// Retrieve the number of right overlaps.
    pub fn num_overlaps(&self) -> usize {
        self.overlaps().len()
    }

    /// Retrieve the right overlaps.
    pub fn overlaps(&self) -> Vec<Position> {
        self.rights
            .iter()
            .map(|r| r.overlap_width(&self.left))
            .collect()
    }

    /// Get the left index.
    pub fn left_index(&self) -> Option<usize> {
        self.left.index()
    }

    /// Get the right indices.
    pub fn right_indices(&self) -> Vec<Option<usize>> {
        self.rights.iter().map(|r| r.index()).collect()
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

    /// Create an iterator over the joins.
    pub fn iter(&'a self) -> JoinDataIterator<'a, DL, DR> {
        JoinDataIterator {
            inner: self.joins.iter(),
            left_data: Some(&self.left_data),
            right_data: Some(self.right_data),
        }
    }
}

pub struct CombinedJoinData<DL, DR> {
    pub join: LeftGroupedJoin, // Information on the join
    pub left_data: DL,         // The left data element
    pub right_data: Vec<DR>,   // The right data elements
}

impl<DL, DR> JoinDataOperations<DL, DR> for CombinedJoinData<DL, DR> {
    type LeftDataElementType = DL;
    type RightDataElementType = DR;
    fn join(&self) -> &LeftGroupedJoin {
        &self.join
    }
    fn left_data(&self) -> Option<&Self::LeftDataElementType> {
        Some(&self.left_data)
    }
    fn right_data(&self) -> Option<&Vec<Self::RightDataElementType>> {
        Some(&self.right_data)
    }
}

impl<'a, DL, DR> JoinData<'a, DL, DR>
where
    DL: IndexedDataContainer + 'a,
    DR: IndexedDataContainer + 'a,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn map<F, V>(self, func: F) -> Vec<V>
    where
        F: Fn(
            CombinedJoinData<
                <DL as IndexedDataContainer>::OwnedItem,
                <DR as IndexedDataContainer>::OwnedItem,
            >,
        ) -> V,
    {
        self.joins
            .into_iter()
            .map(|join| {
                let left_data = self.left_data.get_owned(join.left_index().unwrap());
                let right_data = join
                    .right_indices()
                    .iter()
                    .map(|idx| self.right_data.get_owned(idx.unwrap()))
                    .collect();

                func(CombinedJoinData {
                    join,
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

    /// Create an iterator over the joins.
    pub fn iter(&'a self) -> JoinDataIterator<'a, (), DR> {
        JoinDataIterator {
            inner: self.joins.iter(),
            left_data: None,
            right_data: Some(self.right_data),
        }
    }
}

pub struct CombinedJoinDataLeftEmpty<DR> {
    pub join: LeftGroupedJoin, // Information on the join
    pub right_data: Vec<DR>,   // The right data element
}

impl<DR> JoinDataOperations<(), DR> for CombinedJoinDataLeftEmpty<DR> {
    type LeftDataElementType = ();
    type RightDataElementType = DR;
    fn join(&self) -> &LeftGroupedJoin {
        &self.join
    }
    fn left_data(&self) -> Option<&Self::LeftDataElementType> {
        None
    }
    fn right_data(&self) -> Option<&Vec<Self::RightDataElementType>> {
        Some(&self.right_data)
    }
}

impl<'a, DR> JoinDataLeftEmpty<'a, DR>
where
    DR: IndexedDataContainer + 'a,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn map<F, V>(self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataLeftEmpty<<DR as IndexedDataContainer>::OwnedItem>) -> V,
    {
        // TODO/OPTIMIZE: would consuming here (and analagous funcs) be better/faster?
        // Would require a bit of a redesign.
        self.joins
            .into_iter()
            .map(|join| {
                let right_indices = join.right_indices();
                let right_data = right_indices
                    .iter()
                    .map(|idx| self.right_data.get_owned(idx.unwrap()))
                    .collect();

                func(CombinedJoinDataLeftEmpty { join, right_data })
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

    /// Create an iterator over the joins.
    pub fn iter(&'a self) -> JoinDataIterator<'a, DL, ()> {
        JoinDataIterator {
            inner: self.joins.iter(),
            left_data: Some(&self.left_data),
            right_data: None,
        }
    }
}

pub struct CombinedJoinDataRightEmpty<DL> {
    pub join: LeftGroupedJoin, // Information on the join
    pub left_data: DL,         // The right data element
}

impl<DL> JoinDataOperations<DL, ()> for CombinedJoinDataRightEmpty<DL> {
    type LeftDataElementType = DL;
    type RightDataElementType = ();
    fn join(&self) -> &LeftGroupedJoin {
        &self.join
    }
    fn left_data(&self) -> Option<&Self::LeftDataElementType> {
        Some(&self.left_data)
    }
    fn right_data(&self) -> Option<&Vec<Self::RightDataElementType>> {
        None
    }
}

impl<'a, DL> JoinDataRightEmpty<DL>
where
    DL: IndexedDataContainer,
{
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn map<F, V>(self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataRightEmpty<<DL as IndexedDataContainer>::OwnedItem>) -> V,
    {
        self.joins
            .into_iter()
            .map(|join| {
                let left_data = self.left_data.get_owned(join.left_index().unwrap());

                func(CombinedJoinDataRightEmpty { join, left_data })
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

    /// Create an iterator over the joins.
    pub fn iter(&self) -> JoinDataIterator<'_, (), ()> {
        JoinDataIterator {
            inner: self.joins.iter(),
            left_data: None,
            right_data: None,
        }
    }
}

pub struct CombinedJoinDataBothEmpty {
    pub join: LeftGroupedJoin,
}

impl JoinDataOperations<(), ()> for CombinedJoinDataBothEmpty {
    type LeftDataElementType = ();
    type RightDataElementType = ();
    fn join(&self) -> &LeftGroupedJoin {
        &self.join
    }
    fn left_data(&self) -> Option<&Self::LeftDataElementType> {
        None
    }
    fn right_data(&self) -> Option<&Vec<Self::RightDataElementType>> {
        None
    }
}

impl JoinDataBothEmpty {
    /// Apply `func` to each element, putting the results into the returned
    /// `Vec<U>`.
    pub fn map<F, V>(self, func: F) -> Vec<V>
    where
        F: Fn(CombinedJoinDataBothEmpty) -> V,
    {
        self.joins
            .into_iter()
            .map(|join| func(CombinedJoinDataBothEmpty { join }))
            .collect()
    }
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
    use super::*;
    use crate::ranges::RangeIndexed;

    #[test]
    fn test_join_data_new() {
        let left_data = vec![1, 2];
        let right_data = vec![4, 8];
        let mut jd = JoinData::new(left_data, &right_data);
        assert_eq!(jd.len(), 0);

        let left = RangeIndexed::new(0, 10, 1);
        let mut join = LeftGroupedJoin::new(&left);
        let right = RangeIndexed::new(8, 10, 1);
        join.add_right(&right);
        jd.push(join);
        assert_eq!(jd.len(), 1);
    }

    #[test]
    fn test_single_range_indexed() {
        let ranges = vec![RangeIndexed {
            start: 1,
            end: 5,
            index: 0,
        }];
        let reduced = reduce_ranges(&ranges);

        assert_eq!(reduced.len(), 1);
        assert_eq!(reduced[0].0 .0, 1);
        assert_eq!(reduced[0].0 .1, 5);
        assert!(reduced[0].0 .2.contains(&Some(0)));
    }

    #[test]
    fn test_overlapping_ranges_indexed() {
        let ranges = vec![
            RangeIndexed {
                start: 1,
                end: 4,
                index: 0,
            },
            RangeIndexed {
                start: 3,
                end: 6,
                index: 1,
            },
        ];
        let reduced = reduce_ranges(&ranges);

        assert_eq!(reduced.len(), 3);

        let mut iter = reduced.iter();
        let first_range = iter.next().unwrap();
        assert_eq!(first_range.start(), 1);
        assert_eq!(first_range.end(), 3);
        assert_eq!(first_range.indices().len(), 1);
        assert!(first_range.indices().contains(&Some(0)));

        let second_range = iter.next().unwrap();
        assert_eq!(second_range.start(), 3);
        assert_eq!(second_range.end(), 4);
        assert_eq!(second_range.indices().len(), 2);
        assert!(second_range.indices().contains(&Some(0)));
        assert!(second_range.indices().contains(&Some(1)));

        let third_range = iter.next().unwrap();
        assert_eq!(third_range.start(), 4);
        assert_eq!(third_range.end(), 6);
        assert_eq!(third_range.indices().len(), 1);
        assert!(third_range.indices().contains(&Some(1)));
    }

    #[test]
    fn test_non_overlapping_ranges_indexed() {
        let ranges = vec![
            RangeIndexed {
                start: 1,
                end: 3,
                index: 0,
            },
            RangeIndexed {
                start: 4,
                end: 6,
                index: 1,
            },
        ];
        let reduced = reduce_ranges(&ranges);

        assert_eq!(reduced.len(), 2);

        let mut iter = reduced.iter();
        let first_range = iter.next().unwrap();
        assert_eq!(first_range.start(), 1);
        assert_eq!(first_range.end(), 3);
        assert_eq!(first_range.indices().len(), 1);
        assert!(first_range.indices().contains(&Some(0)));

        let second_range = iter.next().unwrap();
        assert_eq!(second_range.start(), 4);
        assert_eq!(second_range.end(), 6);
        assert_eq!(second_range.indices().len(), 1);
        assert!(second_range.indices().contains(&Some(1)));
    }
}
