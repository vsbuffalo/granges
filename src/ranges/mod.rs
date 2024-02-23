//! Range types and Range Containers.
//!
//!

use crate::{
    error::GRangesError,
    traits::{
        AdjustableGenericRange, GenericRange, GenericRangeOperations, IndexedDataContainer,
        TsvSerialize,
    },
    Position,
};

pub mod coitrees;
pub mod operations;
pub mod vec;

#[derive(Clone, Debug, Default, PartialEq)]
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

/// Represents a genomic range entry with some data.
///
/// This is used as a type for holding a range with associated data directly
/// from a parser.
#[derive(Debug, Clone, PartialEq)]
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

impl TsvSerialize for GenomicRangeRecord<()> {
    fn to_tsv(&self) -> String {
        format!("{}\t{}\t{}", self.seqname, self.start, self.end)
    }
}

// impl<U: TsvSerialize> TsvSerialize for GenomicRangeRecord<Option<U>> {
//     fn to_tsv(&self) -> String {
//         match &self.data {
//             None => {
//                 format!("{}\t{}\t{}", self.seqname, self.start, self.end,)
//             }
//             Some(data) => {
//                 format!(
//                     "{}\t{}\t{}\t{}",
//                     self.seqname,
//                     self.start,
//                     self.end,
//                     data.to_tsv()
//                 )
//             }
//         }
//     }
// }

impl<U: TsvSerialize> TsvSerialize for GenomicRangeRecord<U> {
    /// Serialize this [`GenomicRangeRecord<U>`] to TSV, using the recursive [`TsvSerialize`]
    /// trait method.
    ///
    /// Note that this does *not* append the newline — it is up to
    /// the `U.to_tsv()` trait method implementation. This is design is
    /// intentional, since one common use case is when `U` is an [`Option<String>`],
    /// which is the type combination commonly used in lazy parsing. In this case,
    /// when the data is `None`, that indicates no more columns.
    fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}{}",
            self.seqname,
            self.start,
            self.end,
            self.data.to_tsv()
        )
    }
}

/// Represents a genomic range entry without data, e.g. from a BED3 parser.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomicRangeEmptyRecord {
    pub seqname: String,
    pub start: Position,
    pub end: Position,
}

impl GenomicRangeEmptyRecord {
    pub fn new(seqname: String, start: Position, end: Position) -> Self {
        assert!(end > start);
        Self {
            seqname,
            start,
            end,
        }
    }
}

impl GenericRange for GenomicRangeEmptyRecord {
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

impl AdjustableGenericRange for GenomicRangeEmptyRecord {
    fn set_start(&mut self, start: Position) {
        self.start = start
    }
    fn set_end(&mut self, end: Position) {
        self.end = end
    }
}

impl GenericRangeOperations for GenomicRangeEmptyRecord {
    /// Create flanking regions for this [`GenomicRangeEmptyRecord`] range.
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
                    GenomicRangeEmptyRecord::new(self.seqname.clone(), flank_start, flank_end);
                flanking.push(left_flank_region);
            }
        }
        if let Some(right) = right_flank {
            let flank_start = std::cmp::max(self.end, 0);
            let flank_end = std::cmp::min(self.end + right, seqlen);
            if flank_end > flank_start {
                let right_flank_region =
                    GenomicRangeEmptyRecord::new(self.seqname.clone(), flank_start, flank_end);
                flanking.push(right_flank_region);
            }
        }
        flanking
    }
}

impl TsvSerialize for GenomicRangeEmptyRecord {
    fn to_tsv(&self) -> String {
        format!("{}\t{}\t{}", self.seqname, self.start, self.end)
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
    pub fn seqname(&self, seqnames: &[String]) -> String {
        seqnames[self.seqname_index].clone()
    }
    pub fn to_record<'a, T>(
        self,
        seqnames: &[String],
        data: &'a T,
    ) -> GenomicRangeRecord<<T as IndexedDataContainer<'a>>::Item>
    where
        T: IndexedDataContainer<'a> + TsvSerialize,
    {
        let data = data.get_value(self.index().unwrap());

        GenomicRangeRecord {
            seqname: seqnames[self.seqname_index].clone(),
            start: self.start,
            end: self.end,
            data,
        }
    }
    pub fn to_record_empty<T>(self, seqnames: &[String]) -> GenomicRangeRecord<()> {
        GenomicRangeRecord {
            seqname: seqnames[self.seqname_index].clone(),
            start: self.start,
            end: self.end,
            data: (),
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

/// Validates whether a given range is valid for accessing a sequence of a given `length`.
///
/// # Arguments
///
/// * `range` - The range to validate.
/// * `length` - The length of the sequence.
///
/// # Returns
///
/// * `bool` - `true` if the range is valid for the sequence; otherwise, `false`.
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
}
