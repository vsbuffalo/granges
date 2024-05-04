//! Range types and Range Containers.
//!
//!

use serde::{Deserialize, Serialize};
use std::ops::Range;

use crate::{
    error::GRangesError,
    traits::{AdjustableGenericRange, GenericRange, GenericRangeOperations, IndexedDataContainer},
    Position,
};

pub mod coitrees;
pub mod operations;
pub mod vec;

#[derive(Clone, Debug, Default, PartialEq, Deserialize, Serialize)]
pub struct RangeEmpty {
    pub start: Position,
    pub end: Position,
}

unsafe impl Sync for RangeEmpty {}
unsafe impl Send for RangeEmpty {}

impl RangeEmpty {
    /// Create a new 0-indexed right-exclusive range.
    pub fn new(start: Position, end: Position) -> Self {
        assert!(end > start);
        Self { start, end }
    }
}

impl GenericRange for RangeEmpty {
    fn start(&self) -> Position {
        self.start
    }
    fn end(&self) -> Position {
        self.end
    }
    fn index(&self) -> Option<usize> {
        None
    }
}

impl GenericRangeOperations for RangeEmpty {
    /// Create flanking regions for this [`RangeEmpty`] range.
    fn flanking_ranges<R: GenericRange>(
        &self,
        left_flank: Option<Position>,
        right_flank: Option<Position>,
        seqlen: Position,
    ) -> Vec<Self> {
        let mut flanking = Vec::new();
        if let Some(left) = left_flank {
            let flank_start = std::cmp::max(self.start.saturating_sub(left), 0);
            let flank_end = std::cmp::min(self.start, seqlen);
            if flank_end > flank_start {
                let left_flank_region = RangeEmpty::new(flank_start, flank_end);
                flanking.push(left_flank_region);
            }
        }
        if let Some(right) = right_flank {
            let flank_start = std::cmp::max(self.end, 0);
            let flank_end = std::cmp::min(self.end + right, seqlen);
            if flank_end > flank_start {
                let right_flank_region = RangeEmpty::new(flank_start, flank_end);
                flanking.push(right_flank_region);
            }
        }
        flanking
    }
}

impl AdjustableGenericRange for RangeEmpty {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

impl From<RangeIndexed> for RangeEmpty {
    fn from(value: RangeIndexed) -> Self {
        RangeEmpty {
            start: value.start,
            end: value.end,
        }
    }
}

/// [`RangeIndexed`] is a range with a valid
/// index to a data element in the data container.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RangeIndexed {
    pub start: Position,
    pub end: Position,
    pub index: usize,
}

unsafe impl Sync for RangeIndexed {}
unsafe impl Send for RangeIndexed {}

impl RangeIndexed {
    /// Create a new 0-indexed right-exclusive range.
    pub fn new(start: Position, end: Position, index: usize) -> Self {
        assert!(end > start, "{}-{}", start, end);
        Self { start, end, index }
    }
}

impl GenericRange for RangeIndexed {
    fn start(&self) -> Position {
        self.start
    }
    fn end(&self) -> Position {
        self.end
    }
    fn index(&self) -> Option<usize> {
        Some(self.index)
    }
}

impl GenericRangeOperations for RangeIndexed {
    /// Create flanking regions for this [`RangeIndexed`] range. Note that
    /// the index **will be copied** to the new [`RangeIndexed`] flanking ranges
    /// returned by this method.
    fn flanking_ranges<R: GenericRange>(
        &self,
        left_flank: Option<Position>,
        right_flank: Option<Position>,
        seqlen: Position,
    ) -> Vec<Self> {
        let mut flanking = Vec::new();
        if let Some(left) = left_flank {
            let flank_start = std::cmp::max(self.start.saturating_sub(left), 0);
            let flank_end = std::cmp::min(self.start, seqlen);
            if flank_end > flank_start {
                let left_flank_region = RangeIndexed::new(flank_start, flank_end, self.index);
                flanking.push(left_flank_region);
            }
        }
        if let Some(right) = right_flank {
            let flank_start = std::cmp::max(self.end, 0);
            let flank_end = std::cmp::min(self.end + right, seqlen);
            if flank_end > flank_start {
                let right_flank_region = RangeIndexed::new(flank_start, flank_end, self.index);
                flanking.push(right_flank_region);
            }
        }
        flanking
    }
}

impl AdjustableGenericRange for RangeIndexed {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

/// Represents a genomic range entry with some borrowed data.
/// This is used primarily as a temporary store for deserializing
/// a genomic range.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct GenomicRangeRecordEmptyBorrowed<'a> {
    pub seqname: &'a str,
    pub start: Position,
    pub end: Position,
}

/// Represents a genomic range entry with some borrowed data.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct GenomicRangeRecordBorrowed<'a, U> {
    pub seqname: &'a str,
    pub start: Position,
    pub end: Position,
    pub data: U,
}

/// Represents a genomic range entry with some data.
///
/// This is used as a type for holding a range with associated data directly
/// from a parser. Thus it owns the memory it uses for the seqname and the
/// data.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct GenomicRangeRecord<U> {
    pub seqname: String,
    pub start: Position,
    pub end: Position,
    pub data: U,
}

impl<U> GenomicRangeRecord<U> {
    pub fn new(seqname: String, start: Position, end: Position, data: U) -> Self {
        assert!(end > start);
        Self {
            seqname,
            start,
            end,
            data,
        }
    }

    /// Consume this [`GenomicRangeRecord`], apply a function to its data
    /// and return the new [`GenomicRangeRecord`] with different data.
    pub fn into_map_data<F, V>(self, func: F) -> GenomicRangeRecord<V>
    where
        F: Fn(U) -> V,
    {
        GenomicRangeRecord {
            seqname: self.seqname,
            start: self.start,
            end: self.end,
            data: func(self.data),
        }
    }

    /// Consume this [`GenomicRangeRecord`], drop its data,
    /// and turn it into a [`GenomicRangeRecordEmpty`].
    pub fn into_empty(self) -> GenomicRangeRecordEmpty {
        GenomicRangeRecordEmpty {
            seqname: self.seqname,
            start: self.start,
            end: self.end,
        }
    }
}

impl<U: Clone> GenericRange for GenomicRangeRecord<U> {
    fn start(&self) -> Position {
        self.start
    }
    fn end(&self) -> Position {
        self.end
    }
    fn index(&self) -> Option<usize> {
        None
    }
}

impl<U: Clone> AdjustableGenericRange for GenomicRangeRecord<U> {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

impl<U: Clone> GenericRangeOperations for GenomicRangeRecord<U> {
    /// Create flanking regions for this [`GenomicRangeRecord<U>`] range.
    fn flanking_ranges<R: GenericRange>(
        &self,
        left_flank: Option<Position>,
        right_flank: Option<Position>,
        seqlen: Position,
    ) -> Vec<Self> {
        let mut flanking = Vec::new();
        if let Some(left) = left_flank {
            let flank_start = std::cmp::max(self.start.saturating_sub(left), 0);
            let flank_end = std::cmp::min(self.start, seqlen);
            if flank_end > flank_start {
                let left_flank_region = GenomicRangeRecord::new(
                    self.seqname.clone(),
                    flank_start,
                    flank_end,
                    self.data.clone(),
                );
                flanking.push(left_flank_region);
            }
        }
        if let Some(right) = right_flank {
            let flank_start = std::cmp::max(self.end, 0);
            let flank_end = std::cmp::min(self.end + right, seqlen);
            if flank_end > flank_start {
                let right_flank_region = GenomicRangeRecord::new(
                    self.seqname.clone(),
                    flank_start,
                    flank_end,
                    self.data.clone(),
                );
                flanking.push(right_flank_region);
            }
        }
        flanking
    }
}

/// Represents a genomic range entry without data, e.g. from a BED3 parser.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GenomicRangeRecordEmpty {
    pub seqname: String,
    pub start: Position,
    pub end: Position,
}

impl GenomicRangeRecordEmpty {
    pub fn new(seqname: String, start: Position, end: Position) -> Self {
        assert!(end > start);
        Self {
            seqname,
            start,
            end,
        }
    }
}

impl GenericRange for GenomicRangeRecordEmpty {
    fn start(&self) -> Position {
        self.start
    }
    fn end(&self) -> Position {
        self.end
    }
    fn index(&self) -> Option<usize> {
        None
    }
}

impl AdjustableGenericRange for GenomicRangeRecordEmpty {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

impl GenericRangeOperations for GenomicRangeRecordEmpty {
    /// Create flanking regions for this [`GenomicRangeRecordEmpty`] range.
    fn flanking_ranges<R: GenericRange>(
        &self,
        left_flank: Option<Position>,
        right_flank: Option<Position>,
        seqlen: Position,
    ) -> Vec<Self> {
        let mut flanking = Vec::new();
        if let Some(left) = left_flank {
            let flank_start = std::cmp::max(self.start.saturating_sub(left), 0);
            let flank_end = std::cmp::min(self.start, seqlen);
            if flank_end > flank_start {
                let left_flank_region =
                    GenomicRangeRecordEmpty::new(self.seqname.clone(), flank_start, flank_end);
                flanking.push(left_flank_region);
            }
        }
        if let Some(right) = right_flank {
            let flank_start = std::cmp::max(self.end, 0);
            let flank_end = std::cmp::min(self.end + right, seqlen);
            if flank_end > flank_start {
                let right_flank_region =
                    GenomicRangeRecordEmpty::new(self.seqname.clone(), flank_start, flank_end);
                flanking.push(right_flank_region);
            }
        }
        flanking
    }
}

/// Represents a range entry, with indices to sequence name and possibly data.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomicRangeIndexedRecord {
    pub seqname_index: usize,
    pub start: Position,
    pub end: Position,
    pub index: Option<usize>,
}

impl GenomicRangeIndexedRecord {
    pub fn new(seqname_index: usize, start: Position, end: Position, index: Option<usize>) -> Self {
        assert!(end > start);
        Self {
            seqname_index,
            start,
            end,
            index,
        }
    }
    /// Using the *corresponding ordered* sequence names the `seqname_index` indices
    /// correspond to, get the sequence name.
    pub fn seqname<'a>(&'a self, seqnames: &'a [String]) -> &String {
        &seqnames[self.seqname_index]
    }
    pub fn to_record<'a, T>(
        self,
        seqnames: &'a[String],
        data: &'a T,
    ) -> GenomicRangeRecordBorrowed<'a, <T as IndexedDataContainer>::Item<'a>>
    where
        T: IndexedDataContainer,
    {
        let data = data.get_value(self.index().unwrap());

        GenomicRangeRecordBorrowed {
            seqname: &seqnames[self.seqname_index],
            start: self.start,
            end: self.end,
            data,
        }
    }
    pub fn to_record_empty<'a, T>(self, seqnames: &'a[String]) -> 
        GenomicRangeRecordEmptyBorrowed<'a> {
        GenomicRangeRecordEmptyBorrowed {
            seqname: &seqnames[self.seqname_index],
            start: self.start,
            end: self.end,
        }
    }
}

impl GenericRange for GenomicRangeIndexedRecord {
    fn start(&self) -> Position {
        self.start
    }
    fn end(&self) -> Position {
        self.end
    }
    fn index(&self) -> Option<usize> {
        self.index
    }
}

impl AdjustableGenericRange for GenomicRangeIndexedRecord {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

/// Validates whether a given range is valid for accessing a sequence of a given `length`,
/// raising a [`GRangesError`] if not.
///
/// # Arguments
///
/// * `range` - The range to validate.
/// * `length` - The length of the sequence.
///
pub fn validate_range(
    start: Position,
    end: Position,
    length: Position,
) -> Result<(), GRangesError> {
    if start > end {
        return Err(GRangesError::InvalidGenomicRange(start, end));
    }

    if end >= length {
        return Err(GRangesError::InvalidGenomicRangeForSequence(
            start, end, length,
        ));
    }
    Ok(())
}

/// Try converting genome positions to an right-exclusive [`std::ops::Range`], with
/// checking that the range is valid. This is predominantly used for building [`std::ops::Range`]
/// items that are used to slice (e.g. nucleotide) sequences.
///
/// # Arugments
/// * `first` - the 0-indexed first basepair's position.
/// * `last` - the 0-indexed, last basepair's position (right-inclusive end position).
/// * `len` - the length of the sequence.
///
/// # Returns
/// [`Some(Range<usize>)`] if the conversion could be be done successfully, otherwise
/// returns an error.
pub fn try_range(
    start: Position,
    end: Position,
    length: Position,
) -> Result<Range<usize>, GRangesError> {
    if start >= end {
        return Err(GRangesError::InvalidGenomicRange(start, end));
    }
    if end > length {
        return Err(GRangesError::InvalidGenomicRangeForSequence(
            start, end, length,
        ));
    }
    let start_usize: usize = start.try_into().unwrap();
    let end_usize: usize = end.try_into().unwrap();
    Ok(start_usize..end_usize)
}

#[cfg(test)]
mod tests {
    use super::{validate_range, RangeEmpty};
    use crate::prelude::*;

    #[test]
    fn test_invalid_range_start_end() {
        let result = validate_range(5, 1, 10);
        assert!(matches!(
            result,
            Err(GRangesError::InvalidGenomicRange(5, 1))
        ));
    }

    #[test]
    fn test_valid_range_length() {
        let result = validate_range(1, 10, 11);
        assert!(result.is_ok());
    }

    #[test]
    fn test_invalid_range_length() {
        let result = validate_range(1, 10, 10);
        assert!(matches!(
            result,
            Err(GRangesError::InvalidGenomicRangeForSequence(1, 10, 10))
        ));
    }

    #[test]
    fn test_overlap_range() {
        let range_a = RangeEmpty::new(5, 8);
        let range_b = RangeEmpty::new(4, 6);
        assert_eq!(range_a.overlap_range(&range_b), Some((5, 6)));
    }

    #[test]
    fn test_width() {
        let range_a = RangeEmpty::new(5, 8);
        assert_eq!(range_a.width(), 3);
    }

    #[test]
    fn test_overlap_width() {
        let range_a = RangeEmpty::new(0, 2);
        let range_b = RangeEmpty::new(4, 6);
        assert_eq!(range_a.overlap_width(&range_b), 0);

        let range_a = RangeEmpty::new(0, 2);
        let range_b = RangeEmpty::new(2, 6);
        assert_eq!(range_a.overlap_width(&range_b), 0);

        let range_a = RangeEmpty::new(1, 3);
        let range_b = RangeEmpty::new(2, 5);
        assert_eq!(range_a.overlap_width(&range_b), 1);

        let range_a = RangeEmpty::new(1, 10);
        let range_b = RangeEmpty::new(2, 5);
        assert_eq!(range_a.overlap_width(&range_b), 3);
    }

    #[test]
    fn test_distance_or_overlap() {
        // | 0 | 1 |   |   |
        // |   |   |   |   | 4 | 5 |
        let range_a = RangeEmpty::new(0, 2);
        let range_b = RangeEmpty::new(4, 6);
        assert_eq!(range_a.distance_or_overlap(&range_b), 2);

        // | 0 | 1 |   |   |
        // |   |   | 2 | 3 | 4 | 5 |
        let range_a = RangeEmpty::new(0, 2);
        let range_b = RangeEmpty::new(2, 6);
        assert_eq!(range_a.distance_or_overlap(&range_b), 0);

        // | 0 | 1 | 2 |   |
        // |   |   | 2 | 3 | 4 |   |
        let range_a = RangeEmpty::new(1, 3);
        let range_b = RangeEmpty::new(2, 5);
        assert_eq!(range_a.distance_or_overlap(&range_b), -1);

        // |   | 1 | 2 | 3 | 4 | ...
        // |   |   | 2 | 3 | 4 |   |
        let range_a = RangeEmpty::new(1, 10);
        let range_b = RangeEmpty::new(2, 5);
        assert_eq!(range_a.distance_or_overlap(&range_b), -3);
    }
}
