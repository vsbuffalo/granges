//! Traits used by the GRanges library.
//!

use std::path::PathBuf;

use crate::{
    error::GRangesError,
    granges::GRanges,
    io::parsers::{FilteredRanges, UnwrappedRanges},
    ranges::GenomicRangeRecord,
    Position,
};

/// Traits for [`GRanges`] types that can be modified.
//pub trait GenomicRangesModifiableRanges<C: RangeContainer> {
//    fn adjust_ranges(self, start_delta: PositionOffset, end_delta: PositionOffset) -> Self;
//}

/// The [`AsGRangesRef`] trait improves the ergonomics of working
/// with both [`GRanges`] and [`GRangesEmpty`] function arguments.
///
/// [`GRangesEmpty`]: crate::granges::GRanges
pub trait AsGRangesRef<'a, C, T> {
    fn as_granges_ref(&'a self) -> &'a GRanges<C, T>;
}

/// The [`GenomicRangesTsvSerialize`] trait defines how to convert a [`GRanges<R, T>`]
/// object, for some mix of generic types, to a TSV file.
pub trait GenomicRangesTsvSerialize<'a, C: RangeContainer> {
    /// Output the TSV version of this [`GRanges`] object.
    fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError>;
}

/// The [`GenericRange`] trait defines common functionality for all range types.
pub trait GenericRange: Clone {
    fn start(&self) -> Position;
    fn end(&self) -> Position;
    fn index(&self) -> Option<usize>;
    fn width(&self) -> Position {
        self.end() - self.start()
    }
    /// Calculate how many basepairs overlap this range and other.
    fn overlap_width<R: GenericRange>(&self, other: &R) -> Position {
        let overlap_start = std::cmp::max(self.start(), other.start());
        let overlap_end = std::cmp::min(self.end(), other.end());
        if overlap_start >= overlap_end {
            return 0;
        }
        overlap_end.saturating_sub(overlap_start)
    }

    /// Return a tuple of the range created by an overlap with another range; `None` if no overlap.
    fn overlap_range<R: GenericRange>(&self, other: &R) -> Option<(Position, Position)> {
        let overlap_start = std::cmp::max(self.start(), other.start());
        let overlap_end = std::cmp::min(self.end(), other.end());

        if overlap_start <= overlap_end {
            Some((overlap_start, overlap_end))
        } else {
            None
        }
    }

    /// Return a tuple version of this range.
    fn as_tuple(&self) -> (Position, Position, Option<usize>) {
        (self.start(), self.end(), self.index())
    }
}

/// The [`GenericRangeOperations`] trait extends additional functionality to [`GenericRange`],
/// such as the creation of flanking regions.
pub trait GenericRangeOperations: GenericRange {
    /// Creates the appropriate flanking ranges.
    ///
    /// If the range has an index, this index will be duplicated for any flanking ranges. This
    /// is because it may be the case that the user wishes to use the data for the original
    /// range for these flanking ranges.
    ///
    /// # Stability
    /// The interface of this is likely to change, in order to handle
    /// flanking regions based on strand.
    fn flanking_ranges<R: GenericRange>(
        &self,
        left_flank: Option<Position>,
        right_flank: Option<Position>,
        seqlen: Position,
    ) -> Vec<Self>;
}

/// The [`AdjustableGenericRange`] trait extends additional functionality to adjustable generic ranges.
pub trait AdjustableGenericRange: GenericRange {
    /// Set the start to the specified position.
    fn set_start(&mut self, start: Position);
    /// Set the end to the specified position.
    fn set_end(&mut self, end: Position);
}

/// Defines functionality common to all range containers, e.g. [`VecRanges<R>`] and
/// [`COITrees`].
///
/// [`VecRanges<R>`]: crate::ranges::vec::VecRanges
/// [`COITrees`]: crate::ranges::coitrees::COITrees
pub trait RangeContainer {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn sequence_length(&self) -> Position;
}

/// Marker trait for data container.
/// This is presently not used for much, but is useful to indicate certain
/// type combinations are data containers.
pub trait DataContainer {}

/// Defines functionality for filtering [`GenomicRangeRecord`] entries based on their
/// sequence names. [`GenomicRangeRecord`] are predominantly used in reading in data,
/// so these trait methods simplify excluding or retaining ranges based on what
/// sequence (i.e. chromosome) they are on.
pub trait GeneralRangeRecordIterator<R: GenericRange>:
    Iterator<Item = Result<R, GRangesError>> + Sized
{
    fn retain_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, R>;
    fn exclude_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, R>;
}

/// [`GenomicRangeRecordUnwrappable`] defines functionality for unwrapping
/// values from some sort of iterator over [`Result`] (to handle e.g. parsing errors)
/// containing [`GenomicRangeRecord<Option<String>>`], turning them into
/// [`GenomicRangeRecord<String>`] (i.e. "unwrapping" them).
pub trait GenomicRangeRecordUnwrappable:
    Iterator<Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>> + Sized
{
    /// Try unwrapping a [`GenomicRangeRecord<Option<String>>`] into a
    /// [`GenomicRangeRecord<String>`], raising a [`GRangesError::TryUnwrapDataError`] if
    /// if a `None` is encountered. This is used in lazy parsing when we expect additional
    /// data to be parsed, but there isn't any.
    fn try_unwrap_data(self) -> UnwrappedRanges<Self>;
}

/// The [`IterableRangeContainer`] trait defines common functionality for iterating over
/// the range types in range containers.
pub trait IterableRangeContainer
where
    Self: RangeContainer,
    <Self as IterableRangeContainer>::RangeType: GenericRange,
{
    type RangeType: GenericRange;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = Self::RangeType> + '_>;
}

/// The [`IntoIterableRangesContainer`] trait defines common functionality for *consuming* iterating
/// over the range types in range containers.
pub trait IntoIterableRangesContainer<R> {
    fn into_iter_ranges(self) -> Box<dyn Iterator<Item = R>>;
}

/// The main [`IndexedDataContainer`] type is used to interact
/// with arbitrary *indexed* container of data of type `DataContainer::Item`.
///
/// GRanges provides [`IndexedDataContainer`] trait implementations for:
///
///  - [`ndarray::Array1`]
///  - [`ndarray::Array2`]
///  - [`Vec`]
///  - [`polars::DataFrame`]
///
/// # Panics
/// This will panic if `DataContainer.get_value()` tries to access
/// an index that does not exist. This is an explicit design choice:
/// we could have returned an `Option` or a `Result`. However,
/// when the index is invalid, there is a serious problem in the creation
/// of the `GRanges` object's set of ranges; this should not ever
/// happen in normal operation. Returning `Option`/`Result` just to
/// protect against this edge case would force most `GRanges`
/// functions to return a `Result`. This would clog the API usage, so
/// we opt to just panic.
pub trait IndexedDataContainer<'a>: DataContainer {
    // note: I don't think we can refactor out this lifetime,
    // since it is needed often for associated types
    type Item;
    type Output;
    // this type is needed to reference the core underlying type,
    // eg to handle reference types
    fn is_valid_index(&self, index: usize) -> bool;
    fn get_value(&'a self, index: usize) -> <Self as IndexedDataContainer>::Item;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn invalid_indices(&self, indices: &[usize]) -> Vec<usize> {
        let mut invalid = Vec::new();
        for &index in indices {
            if !self.is_valid_index(index) {
                invalid.push(index);
            }
        }
        invalid
    }

    /// Create a new data container using a set of range indices.
    fn new_from_indices(&self, indices: &[usize]) -> Self::Output;
}

/// Defines how to serialize something to TSV.
pub trait TsvSerialize {
    // Serialize something to a TSV [`String`].
    fn to_tsv(&self) -> String;
}

impl TsvSerialize for &String {
    fn to_tsv(&self) -> String {
        self.to_string()
    }
}

impl TsvSerialize for String {
    fn to_tsv(&self) -> String {
        self.to_string()
    }
}

impl TsvSerialize for Option<String> {
    fn to_tsv(&self) -> String {
        self.as_ref()
            .map_or("".to_string(), |x| format!("\t{}", x.to_tsv()))
    }
}

impl<U: TsvSerialize> TsvSerialize for Vec<U> {
    fn to_tsv(&self) -> String {
        self.iter()
            .map(|x| x.to_tsv())
            .collect::<Vec<_>>()
            .join("\t")
    }
}
