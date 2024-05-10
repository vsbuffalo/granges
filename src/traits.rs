//! Traits used by the GRanges library.
//!

use std::path::PathBuf;

use indexmap::IndexMap;

use crate::{
    data::DatumType,
    error::GRangesError,
    granges::GRanges,
    io::{
        parsers::{FilteredRanges, UnwrappedRanges},
        tsv::TsvConfig,
    },
    join::LeftGroupedJoin,
    prelude::VecRangesIndexed,
    ranges::GenomicRangeRecord,
    Position, PositionOffset,
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

/// The [`LeftOverlaps`] trait provides compile time polymorphic behavior
/// over its associated [`LeftOverlaps::Output`] type and its `Right`
/// generic type.
///
/// # Errors
/// This needs to return [`Result`] because:
///  - [`GRanges::push_range_with_join()`], etc could fail, due to an invalid
///     range, chromosome name, etc.
///  - Taking an empty data container.
pub trait LeftOverlaps<'a, Right> {
    type Output;

    fn left_overlaps(self, right: &'a Right) -> Result<Self::Output, GRangesError>;
}

/// The [`GenomicRangesTsvSerialize`] trait defines how to convert a [`GRanges<R, T>`]
/// object, for some mix of generic types, to a TSV file.
pub trait GenomicRangesTsvSerialize<'a, C: RangeContainer> {
    /// Output the TSV version of this [`GRanges`] object.
    fn write_to_tsv(
        &'a self,
        output: Option<impl Into<PathBuf>>,
        config: &TsvConfig,
    ) -> Result<(), GRangesError>;
}

/// The [`GenericRange`] trait defines common functionality for all range types.
pub trait GenericRange: Clone {
    fn start(&self) -> Position;
    fn end(&self) -> Position;
    fn index(&self) -> Option<usize>;
    fn width(&self) -> Position {
        self.end() - self.start()
    }
    fn midpoint(&self) -> Position {
        (self.start() + self.end()) / 2
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

    /// Calculate how many basepairs overlap this range and other.
    fn has_overlap_with<R: GenericRange>(&self, other: &R) -> bool {
        // TODO/OPTIMIZE inefficient FIXME
        self.overlap_width(other) > 0
    }

    fn distance_or_overlap<R: GenericRange>(&self, other: &R) -> PositionOffset {
        if self.end() > other.start() && self.start() < other.end() {
            // The ranges overlap -- return how much as a negative number
            // TODO/OPTIMIZE: is this slow?
            let overlap: PositionOffset = self.overlap_width(other).try_into().unwrap();
            -overlap
        } else if self.end() == other.start() || self.start() == other.end() {
            // The ranges are bookends (adjacent to each other)
            0
        } else {
            // The ranges do not overlap and are not bookends; calculate the distance
            // If `self` is entirely before `other`, the distance is `other.start() - self.end()`.
            // If `self` is entirely after `other`, the distance is `self.start() - other.end()`.
            // Using `max` to ensure we always get the positive distance regardless of range order
            std::cmp::max(
                other.start().saturating_sub(self.end()),
                self.start().saturating_sub(other.end()),
            ) as PositionOffset
        }
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

/// The [`GenericGenomicRange`] extends sequence name comparison and related
/// functionality to [`GenericRange`].
// TODO: to do this right and have it be useful we need to push
// the seqname indices lower into parsing, e.g. so this works for key types
// like GenericGenomicRange.
pub trait GenericGenomicRange: GenericRange {
    /// Returns the sequence name index associated with the range
    fn seqname_index(&self) -> usize;

    /// Calculate how many basepairs overlap this range and another, considering seqname_index
    fn genomic_overlap_width<R: GenericGenomicRange>(&self, other: &R) -> Position {
        if self.seqname_index() != other.seqname_index() {
            return 0; // No overlap if seqname_index is different
        }

        let overlap_start = std::cmp::max(self.start(), other.start());
        let overlap_end = std::cmp::min(self.end(), other.end());

        if overlap_start >= overlap_end {
            0
        } else {
            overlap_end.saturating_sub(overlap_start)
        }
    }

    /// Return a tuple of the range created by an overlap with another range, considering seqname_index
    fn genomic_overlap_range<R: GenericGenomicRange>(
        &self,
        other: &R,
    ) -> Option<(Position, Position)> {
        if self.seqname_index() != other.seqname_index() {
            return None; // No overlap if seqname_index is different
        }

        let overlap_start = std::cmp::max(self.start(), other.start());
        let overlap_end = std::cmp::min(self.end(), other.end());

        if overlap_start <= overlap_end {
            Some((overlap_start, overlap_end))
        } else {
            None
        }
    }
}

/// The [`JoinDataOperations`] trait unifies common operations
/// over combined join data types ([`CombinedJoinData`],
/// [`CombinedJoinDataBothEmpty`], etc).
///
///
/// [`JoinDataOperations`]: crate::traits::JoinDataOperations
/// [`CombinedJoinData`]: crate::join::CombinedJoinData
/// [`CombinedJoinDataBothEmpty`]: crate::join::CombinedJoinDataBothEmpty
pub trait JoinDataOperations<DL, DR> {
    type LeftDataElementType;
    type RightDataElementType;

    fn join(&self) -> &LeftGroupedJoin;
    fn left_data(&self) -> Option<&Self::LeftDataElementType>;
    fn right_data(&self) -> Option<&Vec<Self::RightDataElementType>>;
}

/// Helper trait to convert selected fields into `DataType`
pub trait IntoDatumType {
    fn into_data_type(self) -> DatumType;
}

// TODO
pub trait Selection {
    fn select_by_name(&self, name: &str) -> DatumType;
    fn select(&self, names: &[String]) -> Vec<DatumType> {
        names.iter().map(|name| self.select_by_name(name)).collect()
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
    type InternalRangeType; // the internal stored range type
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
    type RangeType: GenericRange; // the iterator range type
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
///  - [`Vec`], this is the most common and general data container.
///  - [ndarray](https://docs.rs/ndarray/latest/ndarray/index.html) `Array1` and `Array2`, if `--features=ndarray` is set.
///  - [polars](https://pola.rs) datafames.
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
pub trait IndexedDataContainer: DataContainer {
    // note: I don't think we can refactor out this lifetime,
    // since it is needed often for associated types
    type Item<'a>
    where
        Self: 'a;
    type OwnedItem;
    // type Output;
    // this type is needed to reference the core underlying type,
    // eg to handle reference types
    fn is_valid_index(&self, index: usize) -> bool;
    fn get_value(&self, index: usize) -> <Self as IndexedDataContainer>::Item<'_>;
    fn get_owned(&self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem;
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

    ///// Create a new data container using a set of range indices.
    // fn new_from_indices(&self, indices: &[usize]) -> Self::Output;
}

/// The Sequences trait defines generic functionality for
/// per-basepair data, e.g. nucleotide sequences or some
/// per-basepair numeric score.
pub trait Sequences {
    type Container<'a>: 'a
    where
        Self: 'a;
    type Slice<'a>;

    fn seqnames(&self) -> Vec<String>;
    fn get_sequence(&self, seqname: &str) -> Result<Self::Container<'_>, GRangesError>;
    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError>;

    /// Apply a function on a [`Sequences::Slice`] of a sequence.
    ///
    /// # Arguments
    /// * `func` - a function that takes a `Self::Slice`, `start`, and `end` positions
    ///            and returns a [`Result<V, GRangesError>`].
    /// * `seqname` - sequence name.
    /// * `start` - a [`Position`] start position.
    /// * `end` - a [`Position`] *inclusive* end position.
    ///
    /// If you're implementing this trait method, it *must* being with:
    ///
    /// ```no_run
    /// use granges::prelude::try_range;
    /// // let seq = self.get_sequence(seqname)?; // e.g. this normally
    /// let seq = "GATAGAGAGTAGAGTA";
    /// let range = try_range(0, 10, seq.len().try_into().unwrap()).unwrap();
    /// ```
    ///
    /// to validate the range and to avoid panics.
    ///
    /// Note that the `start` and `end` positions can often be ignored 
    /// by the function processing the `Slice`. However, this information
    /// is useful when functions need to know the slice coordinates to 
    /// e.g. combine with other data in this region.
    ///
    fn region_map<V, F>(
        &self,
        func: &F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: Fn(<Self as Sequences>::Slice<'_>, (&str, Position, Position)) -> V;

    fn seqlens(&self) -> Result<IndexMap<String, Position>, GRangesError> {
        let mut seqlens = IndexMap::new();
        for seqname in self.seqnames() {
            let seqlen = self.get_sequence_length(&seqname)?;
            seqlens.insert(seqname, seqlen);
        }
        Ok(seqlens)
    }

    /// Create a new [`GRanges<C, Vec<V>>`] by apply the function `func` on
    /// the genomic ranges from `granges`.
    ///
    /// # Arguments
    fn region_map_into_granges<'b, C, F, V, T: 'b>(
        &self,
        granges: &'b impl AsGRangesRef<'b, C, T>,
        func: &F,
    ) -> Result<GRanges<VecRangesIndexed, Vec<V>>, GRangesError>
    where
        V: Clone,
        C: IterableRangeContainer + 'b,
        F: Fn(<Self as Sequences>::Slice<'_>, (&str, Position, Position)) -> V,
    {
        let granges_ref = granges.as_granges_ref();
        let seqlens = &granges_ref.seqlens().clone();
        let mut gr: GRanges<VecRangesIndexed, Vec<V>> = GRanges::new_vec(seqlens);
        for (seqname, ranges) in granges_ref.ranges.iter() {
            // unwrap should be safe, since seqname is produced from ranges iterator.
            for range in ranges.iter_ranges() {
                let (start, end) = (range.start(), range.end());
                let value = self.region_map(&func, seqname, start, end)?;
                gr.push_range(seqname, start, end, value.clone())?;
            }
        }
        Ok(gr)
    }
}
