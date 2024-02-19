//! Traits used by the GRanges library.
//!

use std::path::PathBuf;

use indexmap::IndexMap;

use crate::{
    error::GRangesError, io::parsers::FilteredRanges, ranges::GenomicRangeRecord, Position, granges::GRanges, PositionOffset,
};

/// Traits for [`GRanges`] types that can be modified.
pub trait GenomicRangesOperationsModifiable<C: RangeContainer> {
    fn adjust_ranges(self, start_delta: PositionOffset, end_delta: PositionOffset) -> Self;
}

pub trait GenomicRangesTsvSerialize<'a, C: RangeContainer> {
    /// Output the TSV version of this [`GRanges`] object.
    fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError>;
}

/// Traits for [`GRanges`] types that can be built from an iterator.
pub trait GenomicRangesOperationsExtended<C: RangeContainer> {
    type DataContainerType;
    type DataElementType;
    /// Build a new [`GRanges`] object from an iterator.
    fn from_iter<I>(iter: I, seqlens: &IndexMap<String, Position>) -> Result<GRanges<C, Self::DataContainerType>, GRangesError> where I: Iterator<Item=Result<GenomicRangeRecord<Self::DataElementType>, GRangesError>>;
}

pub trait GenomicRangesOperations<C: RangeContainer> {
   /// Get the total number of ranges.
    fn len(&self) -> usize;

    /// Return whether the [`GRanges`] object is empty (contains no ranges).
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the raw range container.
    fn get_ranges(&self, seqname: &str) -> Option<&C>;

    /// Get the sequence names.
    fn seqnames(&self) -> Vec<String>;

    /// Get the sequences lengths.
    fn seqlens(&self) -> IndexMap<String, Position>;
}

///
pub trait GenericRange: Clone {
    fn start(&self) -> Position;
    fn end(&self) -> Position;
    fn index(&self) -> Option<usize>;
    fn set_start(&mut self, start: Position);
    fn set_end(&mut self, end: Position);
    fn width(&self) -> Position {
        return self.end() - self.start() - 1;
    }
    /// Calculate how many basepairs overlap this range and other.
    fn overlap_width<R: GenericRange>(&self, other: &R) -> Position {
        let overlap_start = std::cmp::max(self.start(), other.start());
        let overlap_end = std::cmp::min(self.end(), other.end());
        std::cmp::max(overlap_end - overlap_start - 1, 0)
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
}

/// Defines functionality common to all range containers, e.g. [`VecRanges<R>`] and
/// [`COITrees`].
pub trait RangeContainer {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn sequence_length(&self) -> Position;
}

pub trait DataContainer{}

/// Defines functionality for filtering [`RangeRecord`] entries based on their
/// sequence names. [`RangeRecord`] are predominantly used in reading in data,
/// so these trait methods simplify excluding or retaining ranges based on what
/// sequence (i.e. chromosome) they are on.
pub trait GeneralRangeRecordIterator<U>:
    Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>> + Sized
{
    fn retain_seqnames(self, seqnames: Vec<String>) -> FilteredRanges<Self, U>;
    fn exclude_seqnames(self, seqnames: Vec<String>) -> FilteredRanges<Self, U>;
}

/// The [`RangesIterable`] trait defines common functionality for iterating over
/// the range types in range containers.
pub trait RangesIterable
where
    Self: RangeContainer,
{
    type RangeType: GenericRange;
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = Self::RangeType> + '_>;
}

/// The [`RangesIntoIterable`] trait defines common functionality for *consuming* iterating
/// over the range types in range containers.
pub trait RangesIntoIterable<R> {
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

impl<U: TsvSerialize> TsvSerialize for Vec<U> {
    fn to_tsv(&self) -> String {
        self.iter()
            .map(|x| x.to_tsv())
            .collect::<Vec<_>>()
            .join("\t")
    }
}
