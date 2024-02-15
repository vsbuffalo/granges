//! Traits used by the GRanges library.
//!

use crate::{error::GRangesError, io::parsers::FilteredIntervals, ranges::RangeRecord};


/// Defines functionality common to all range containers, e.g. [`VecRanges<R>`] and
/// [`COITrees`]. 
pub trait RangeContainer {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Defines functionality for filtering [`RangeRecord`] entries based on their
/// sequence names. [`RangeRecord`] are predominantly used in reading in data,
/// so these trait methods simplify excluding or retaining ranges based on what
/// sequence (i.e. chromosome) they are on.
pub trait GeneralRangeRecordIterator<U>:
    Iterator<Item = Result<RangeRecord<U>, GRangesError>> + Sized
{
    fn retain_seqnames(self, seqnames: Vec<String>) -> FilteredIntervals<Self, U>;
    fn exclude_seqnames(self, seqnames: Vec<String>) -> FilteredIntervals<Self, U>;
}

/// The [`RangesIterable`] trait defines common functionality for iterating over
/// the range types in range containers.
pub trait RangesIterable<R> {
    fn iter_ranges(&self) -> Box<dyn Iterator<Item = R> + '_>;
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
pub trait IndexedDataContainer<'a> {
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
