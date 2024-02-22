//! # Design
//!
//!
//! # [`GRanges<R, T>`] Generic Types
//!
//! [`GRanges<R, T>`] types are generic over:
//!
//! 1. Their **range container** (`R`). This is because different operations to be fast, and thus
//!    need different data structures. The [`GRanges`] methods allow for efficient conversion of
//!    one range type to another.
//!
//! 2. The *optional* **data container** (`T`). The data container exists if the ranges have some
//!    associated data. When range *do* have data, there is an index that associates each range
//!    with its data element in the data container.
//!
//! This brings up the most important thing to know about working with the GRanges library: because
//! of the emphasis on on knowing all types at compile-time (which has performance benefits over
//! runtime type polymorphism), it must be known at compile-time whether ranges have data or not.
//!
//! In most applications this is known: a user specifies, for example, a GTF file and the contents
//! are parsed and processed accordingly. However, it's not unfeasible to imagine that one would
//! need runtime "polymorphism" over BED3 input (which would lead to a [`GRanges`] object without
//! data) and BED* (e.g. BED5, BED12, etc) input. (For example, the GRanges command line tool
//! `granges` runs into this problem — see it's implementation for examples.)
//!
//! These two possibilities are handled with two differently typed parsing iterators:
//! [`Bed3Iterator`] and [`BedlikeIterator`]. These yield different parsed range types,
//! [`GenomicRangeEmptyRecord`] and [`GenomicRangeRecord`], respectively.
//!
//!
//!
//! This is an important concept when working with [`GRanges<R, T>`] types:
//!
//! **High-level data types**: A key feature of GRanges design is that a single type of ranges are
//! contained in the range containers. By knowing that every range in a range container either has
//! an index to a data element in the data container or it does not ahead of time simplifies
//! downstream ergonomics tremendously.
//!
//! **Emphasis on compile-time**: For this, let's consider a common problem: a bioinformatics tool
//! needs to read in a BED-like file that has a variable, unknown at compile time, number of
//! columns.
//!
//! In Rust, this could be handled in one of two ways. First, it could be handled at *runtime*, by
//! leveraging Rust's dynamic trait system. For example, imagine loading in one of two possible BED
//! formats:
//!
//!  1. *Data-less BED3*: The appropriate `GRanges` object here would be a `GRanges<VecRangesEmpty,
//!     ()>`.
//!
//!  2. *BED-like with data*: Here, we'd need a `GRanges<VecRangesIndexed, Vec<U>>`, where the
//!     `Vec<U>` is data container containing just-loaded-in data.
//!
//! Suppose your code doesn't know, when you're writing it, which of these two cases it will
//! encounter.
//!
//! Because at compile-time, the types *need* to be known, there are a few options here.
//!

use std::{collections::HashSet, path::PathBuf};

use coitrees::IntervalNode;
use genomap::GenomeMap;
use indexmap::IndexMap;

use crate::{
    ensure_eq,
    io::OutputFile,
    iterators::GRangesIterator,
    join::{JoinData, LeftGroupedJoin},
    prelude::GRangesError,
    ranges::{
        coitrees::{COITrees, COITreesEmpty, COITreesIndexed},
        vec::{VecRanges, VecRangesEmpty, VecRangesIndexed},
        GenomicRangeEmptyRecord, GenomicRangeRecord, RangeEmpty, RangeIndexed,
    },
    traits::{
        AdjustableGenericRange, AsGRangesRef, GenericRange, GenericRangeOperations,
        GenomicRangesTsvSerialize, IndexedDataContainer, IterableRangeContainer, RangeContainer,
        TsvSerialize,
    },
    Position, PositionOffset,
};

#[derive(Clone, Debug)]
pub struct GRanges<C, T> {
    pub(crate) ranges: GenomeMap<C>,
    pub(crate) data: Option<T>,
}

#[derive(Clone, Debug)]
pub struct GRangesEmpty<C>(GRanges<C, ()>);

impl<C> GRangesEmpty<C>
where
    C: RangeContainer,
{
    /// Get the total number of ranges.
    pub fn len(&self) -> usize {
        self.0.ranges.values().map(|ranges| ranges.len()).sum()
    }

    /// Return whether the [`GRanges`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.0.len() == 0
    }

    /// Get the raw range container.
    pub fn get_ranges(&self, seqname: &str) -> Option<&C> {
        self.0.ranges.get(seqname)
    }

    /// Get the sequence names.
    pub fn seqnames(&self) -> Vec<String> {
        self.0.ranges.names()
    }

    /// Get the sequences lengths.
    pub fn seqlens(&self) -> IndexMap<String, Position> {
        let seqlens = self
            .0
            .ranges
            .iter()
            .map(|(seqname, ranges)| (seqname.to_string(), ranges.sequence_length()))
            .collect();
        seqlens
    }
}

impl<C> From<GRangesEmpty<C>> for GRanges<C, ()> {
    fn from(value: GRangesEmpty<C>) -> Self {
        value.0
    }
}

impl<'a, C> AsGRangesRef<'a, C, ()> for GRangesEmpty<C> {
    /// Convert a reference to a [`GRangesEmpty<C>`] to a reference to the
    /// underlying [`GRanges<C, ()>`]. This is to greatly improve the ergonomics
    /// of functions that could take either a [`GRanges`] or [`GRangesEmpty] type.
    fn as_granges_ref(&'a self) -> &'a GRanges<C, ()> {
        &self.0
    }
}

impl<'a, C, T> AsGRangesRef<'a, C, T> for GRanges<C, T> {
    /// Return a reference of a [`GRanges<C, T>`] object. This is essentially
    /// a pass-through method. [`IntoGRangesRef`] is not needed in this case,
    /// but is needed elsewhere (see the implementation for [`GRangesEmpty`]) to
    /// improve the ergonomics of working with [`GRanges`] and [`GRangesEmpty`] types.
    fn as_granges_ref(&'a self) -> &'a GRanges<C, T> {
        self
    }
}

impl<C, T> GRanges<C, T>
where
    C: RangeContainer,
{
    /// Get the total number of ranges.
    pub fn len(&self) -> usize {
        self.ranges.values().map(|ranges| ranges.len()).sum()
    }

    /// Return whether the [`GRanges`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the raw range container.
    pub fn get_ranges(&self, seqname: &str) -> Option<&C> {
        self.ranges.get(seqname)
    }

    /// Get the sequence names.
    pub fn seqnames(&self) -> Vec<String> {
        self.ranges.names()
    }

    /// Get the sequences lengths.
    pub fn seqlens(&self) -> IndexMap<String, Position> {
        let seqlens = self
            .ranges
            .iter()
            .map(|(seqname, ranges)| (seqname.to_string(), ranges.sequence_length()))
            .collect();
        seqlens
    }
}

impl<'a, T> GenomicRangesTsvSerialize<'a, VecRangesIndexed> for GRanges<VecRangesIndexed, T>
where
    T: IndexedDataContainer<'a>,
    T: TsvSerialize,
    <T as IndexedDataContainer<'a>>::Item: TsvSerialize,
{
    /// Write
    fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError> {
        // output stream -- header is None for now (TODO)
        let output = output.map_or(OutputFile::new_stdout(None), |file| {
            OutputFile::new(file, None)
        });
        let mut writer = output.writer()?;

        let data_ref = self.data.as_ref().ok_or(GRangesError::NoDataContainer)?;
        let seqnames = self.seqnames();
        for range in self.iter_ranges() {
            let record = range.to_record(&seqnames, data_ref);
            writeln!(writer, "{}", record.to_tsv())?;
        }
        Ok(())
    }
}

impl<'a, R: IterableRangeContainer> GenomicRangesTsvSerialize<'a, R> for GRangesEmpty<R> {
    /// Output a BED3 file for for this data-less [`GRanges<R, ()>`].
    fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError> {
        // output stream -- header is None for now (TODO)
        let output = output.map_or(OutputFile::new_stdout(None), |file| {
            OutputFile::new(file, None)
        });
        let mut writer = output.writer()?;

        let seqnames = self.seqnames();
        for range in self.0.iter_ranges() {
            let record = range.to_record_empty::<()>(&seqnames);
            writeln!(writer, "{}", record.to_tsv())?;
        }
        Ok(())
    }
}

impl<R: GenericRange, T> GRanges<VecRanges<R>, T> {
    /// Create a new [`GRanges`] object, with vector storage for ranges and data.
    ///
    /// This combination of range and data containers is used when loading data into
    /// a new [`GRanges`] object, and the size cannot be known beforehand. Rust's
    /// [`Vec`] will dynamically grow to accommodate new ranges; use [`GRanges.shrink()`]
    /// call the [`Vec`]'s shrink to size methods on the range and data containers
    /// after data loading to shrink to the minimal necessary size (this can reduce
    /// memory usage).
    pub fn new_vec(seqlens: &IndexMap<String, Position>) -> Self {
        let mut ranges = GenomeMap::new();
        for (seqname, length) in seqlens.iter() {
            // this should never happen because the error is only if
            // insert encounters a seqname that's already been inserted -- that
            // cannot happen here.
            ranges
                .insert(seqname, VecRanges::new(*length))
                .expect("Internal error: please report");
        }
        Self { ranges, data: None }
    }

    /// Consume this [`GRanges`] object and sort the ranges.
    pub fn sort(mut self) -> Self {
        self.ranges.values_mut().for_each(|ranges| ranges.sort());
        self
    }

    pub fn shink(&mut self) {
        todo!()
    }
}

impl<R: AdjustableGenericRange, T> GRanges<VecRanges<R>, T> {
    /// Adjust all the ranges in this [`GRanges`] object in place.
    pub fn adjust_ranges(mut self, start_delta: PositionOffset, end_delta: PositionOffset) -> Self {
        self.ranges
            .values_mut()
            .for_each(|ranges| ranges.adjust_ranges(start_delta, end_delta));
        self
    }
}

impl<T> GRanges<VecRangesIndexed, T>
where
    VecRangesIndexed: IterableRangeContainer,
    <VecRangesIndexed as IterableRangeContainer>::RangeType: GenericRangeOperations,
{
    /// Create a new [`GRanges`] object of the flanking ranges of the specified widths.
    ///
    /// # ⚠️ Warnings
    /// This creates an index-only ranges container, meaning the data container is *not*
    /// cloned into this object. This is because future version of the GRanges library
    /// will have special reference counting based data sharing, to reduce memory overhead.
    /// Until then, note that trying to access data elements with these indices may cause
    /// panics.
    ///
    /// # Panics
    /// Data accessing operations in the [`GRanges`] object returned by this function
    /// will likely cause a panic — see note above.
    pub fn flanking_ranges(
        &self,
        left: Option<Position>,
        right: Option<Position>,
    ) -> Result<Self, GRangesError> {
        let mut gr: GRanges<VecRangesIndexed, T> = GRanges::new_vec(&self.seqlens());
        let seqlens = self.seqlens();
        for (seqname, ranges) in self.ranges.iter() {
            // unwrap should be safe, since seqname is produced from ranges iterator.
            let seqlen = seqlens.get(seqname).unwrap();
            for range in ranges.iter_ranges() {
                let flanking_ranges = range.flanking_ranges::<RangeIndexed>(left, right, *seqlen);
                for flanking_range in flanking_ranges {
                    gr.push_range_with_index(
                        seqname,
                        flanking_range.start,
                        flanking_range.end,
                        flanking_range.index,
                    )?;
                }
            }
        }
        Ok(gr)
    }
}

impl<R: GenericRange> GRangesEmpty<VecRanges<R>> {
    /// Create a new [`GRangesEmpty`] object, with vector storage for ranges and no
    /// data container.
    pub fn new_vec(seqlens: &IndexMap<String, Position>) -> Self {
        GRangesEmpty(GRanges::new_vec(seqlens))
    }

    /// Sort the ranges by position for this [`GRangesEmpty`] object.
    ///
    /// This operation is consuming and returns the sorted [`GRangesEmpty`] object.
    pub fn sort(self) -> Self {
        GRangesEmpty(self.0.sort())
    }

    pub fn shink(&mut self) {
        todo!()
    }
}

impl GRangesEmpty<VecRangesEmpty> {
    /// Make a [`GRangesEmpty<VecRanges>`] with ranges from (possibly overlapping) windows.
    ///
    /// # Arguments
    ///  * `seqlens`: the sequence (e.g. chromosome) lengths.
    ///  * `width`: the window width, in basepairs.
    ///  * `chop`: whether to cut off the last window, if there is a remainder less than the width.
    ///  * `step`: the step length, in basepairs; if None, step is `width`.
    ///
    /// # Examples
    ///
    /// ```
    /// use granges::prelude::*;
    ///
    /// let sl = seqlens!( "chr1" => 11);
    ///
    /// // no step don't chop off remainder
    /// let gr = GRangesEmpty::from_windows(&sl, 5, None, false).unwrap();
    ///
    /// let mut range_iter = gr.iter_ranges();
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (0,  5, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (5,  10, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (10,  11, None));
    ///
    /// // no step do chop off remainder
    /// let gr = GRangesEmpty::from_windows(&sl, 5, None, true).unwrap();
    ///
    /// let mut range_iter = gr.iter_ranges();
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (0,  5, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (5,  10, None));
    ///
    /// // with step don't chop off remainder
    /// let gr = GRangesEmpty::from_windows(&sl, 5, Some(2), false).unwrap();
    ///
    /// let mut range_iter = gr.iter_ranges();
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (0,  5, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (2,  7, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (4,  9, None));
    /// assert_eq!(range_iter.next().unwrap().as_tuple(), (6,  11, None));
    ///
    /// ```
    pub fn from_windows(
        seqlens: &IndexMap<String, Position>,
        width: u32,
        step: Option<u32>,
        chop: bool,
    ) -> Result<GRangesEmpty<VecRangesEmpty>, GRangesError> {
        let mut gr: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(seqlens);

        // have we encountered a remainder chunk?

        let mut remainder = false;
        // iterate over each chromosome and create windows
        for (seqname, len) in seqlens {
            let mut start = 0;
            while start < *len {
                let mut end = start + width - 1;

                if end >= *len {
                    if chop {
                        break;
                    } else {
                        end = std::cmp::min(end, len - 1);
                    }
                }
                if end - start + 1 < width {
                    if remainder {
                        break;
                    }
                    remainder = true;
                }
                gr.push_range(seqname, start, end + 1)?;
                start += step.unwrap_or(width);
            }
        }
        Ok(gr)
    }
}

impl<R: AdjustableGenericRange> GRangesEmpty<VecRanges<R>> {
    pub fn adjust_ranges(self, start_delta: PositionOffset, end_delta: PositionOffset) -> Self {
        GRangesEmpty(self.0.adjust_ranges(start_delta, end_delta))
    }
}

impl GRangesEmpty<VecRangesEmpty>
where
    VecRangesEmpty: IterableRangeContainer,
    <VecRangesEmpty as IterableRangeContainer>::RangeType: GenericRangeOperations,
{
    /// Create a new [`GRanges`] object of the flanking ranges of the specified widths.
    ///
    /// # ⚠️ Warnings
    /// This creates an index-only ranges container, meaning the data container is *not*
    /// cloned into this object. This is because future version of the GRanges library
    /// will have special reference counting based data sharing, to reduce memory overhead.
    /// Until then, note that trying to access data elements with these indices may cause
    /// panics.
    ///
    /// # Panics
    /// Data accessing operations in the [`GRanges`] object returned by this function
    /// will likely cause a panic — see note above.
    pub fn flanking_ranges(
        &self,
        left: Option<Position>,
        right: Option<Position>,
    ) -> Result<Self, GRangesError> {
        let mut gr: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(&self.seqlens());
        let seqlens = self.seqlens();
        for (seqname, ranges) in self.0.ranges.iter() {
            // unwrap should be safe, since seqname is produced from ranges iterator.
            let seqlen = seqlens.get(seqname).unwrap();
            for range in ranges.iter_ranges() {
                let flanking_ranges = range.flanking_ranges::<RangeIndexed>(left, right, *seqlen);
                for flanking_range in flanking_ranges {
                    gr.push_range(seqname, flanking_range.start, flanking_range.end)?;
                }
            }
        }
        Ok(gr)
    }
}

impl<U> GRanges<VecRangesIndexed, Vec<U>> {
    /// Push a genomic range with its data to the range and data containers in a [`GRanges] object.
    pub fn push_range(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
        data: U,
    ) -> Result<(), GRangesError> {
        if self.data.is_none() {
            self.data = Some(Vec::new());
        }
        let data_ref = self.data.as_mut().ok_or(GRangesError::NoDataContainer)?;
        // push data to the vec data container, getting the index
        let index: usize = {
            data_ref.push(data);
            data_ref.len() - 1 // new data index
        };
        // push an indexed range
        let range = RangeIndexed::new(start, end, index);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<'a, DL, DR> GRanges<VecRangesIndexed, JoinData<'a, DL, DR>> {
    /// Push a genomic range with its data to the range and data containers in a [`GRanges] object.
    ///
    /// Note that this has slightly different behavior than other [`GRanges.push_range()`]
    /// methods, in that it requires that the [`JoinData`] object be initialized first.
    /// This is because were this not to be the case, each call to `push_ranges()` would
    /// require references to the right and left data containers. This is not ergonomic.
    ///
    /// # Panics
    /// This will panic if the developer did not initialize the new [`JoinData`]
    /// correctly. This is not handled with a user error since it reflects a
    /// developer error.
    pub fn push_range_with_join(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
        data: LeftGroupedJoin,
    ) -> Result<(), GRangesError> {
        if self.data.is_none() {
            // unlike push_range()
            panic!("Internal error: JoinData not initialized.");
        }
        let data_ref = self.data.as_mut().ok_or(GRangesError::NoDataContainer)?;
        // push data to the vec data container, getting the index
        let index: usize = {
            data_ref.push(data);
            data_ref.len() - 1 // new data index
        };
        // push an indexed range
        let range = RangeIndexed::new(start, end, index);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<'a, T> GRanges<VecRanges<RangeIndexed>, T>
where
    T: IndexedDataContainer<'a>,
    T: TsvSerialize,
    <T as IndexedDataContainer<'a>>::Item: TsvSerialize,
{
    /// Write this [`GRanges<VecRanges, T>`] object to a TSV file, using the
    /// [`TsvSerialize`] trait methods defiend for the items in `T`.
    pub fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError> {
        // output stream -- header is None for now (TODO)
        let output = output.map_or(OutputFile::new_stdout(None), |file| {
            OutputFile::new(file, None)
        });
        let mut writer = output.writer()?;

        let seqnames = self.seqnames();
        let data_ref = self.data.as_ref().ok_or(GRangesError::NoDataContainer)?;
        for range in self.iter_ranges() {
            let record = range.to_record(&seqnames, data_ref);
            writeln!(writer, "{}", record.to_tsv())?;
        }
        Ok(())
    }

    /// Conduct a left overlap join, consuming self and returning a new
    /// [`GRanges<VecRangesIndexed, JoinData>>`].
    ///
    /// The [`JoinData`] container contains references to both left and right
    /// data containers and a [`Vec<OverlapJoin>`]. Each [`OverlapJoin`] represents
    /// a summary of an overlap, which downstream operations use to calculate
    /// statistics using the information about overlaps.
    pub fn left_overlaps<M: Clone + 'a, DR: 'a>(
        self,
        right: &'a impl AsGRangesRef<'a, COITrees<M>, DR>,
    ) -> Result<GRanges<VecRanges<RangeIndexed>, JoinData<'a, T, DR>>, GRangesError>
    where
        IntervalNode<M, usize>: GenericRange,
    {
        let mut gr: GRanges<VecRangesIndexed, JoinData<'a, T, DR>> =
            GRanges::new_vec(&self.seqlens());

        let right_ref = right.as_granges_ref();
        gr.data = Some(JoinData::new(self.data, right_ref.data.as_ref()));

        for (seqname, left_ranges) in self.ranges.iter() {
            for left_range in left_ranges.iter_ranges() {
                // Left join: every left range gets a JoinData.
                let mut join_data = LeftGroupedJoin::new(&left_range);
                if let Some(right_ranges) = right_ref.ranges.get(seqname) {
                    right_ranges.query(left_range.start(), left_range.end(), |right_range| {
                        // NOTE: right_range is a coitrees::IntervalNode.
                        join_data.add_right(&left_range, right_range);
                    });
                }
                gr.push_range_with_join(seqname, left_range.start, left_range.end, join_data)?;
            }
        }
        Ok(gr)
    }
}

impl<'a, C, T> GRanges<C, T>
where
    T: IndexedDataContainer<'a>,
{
    /// Get the data in the data container at specified index.
    ///
    /// # Panics
    /// This will panic if there if the index is invalid, or the
    /// data container is `None`. Both of these indicate internal
    /// design errors: please file an issue of you encounter a panic.
    pub fn get_data_value(&'a self, index: usize) -> <T as IndexedDataContainer>::Item {
        self.data
            .as_ref()
            .expect("data container was None")
            .get_value(index)
    }
}

impl GRangesEmpty<VecRanges<RangeEmpty>> {
    /// Push an empty range (no data) to the appropriate [`VecRangesEmpty`] range container.
    pub fn push_range(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<(), GRangesError> {
        // push an unindexed (empty) range
        let range = RangeEmpty::new(start, end);
        let range_container = self
            .0
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<T> GRanges<VecRanges<RangeIndexed>, T> {
    /// Push an empty range (no data) to the [`VecRangesEmpty`] range container.
    pub fn push_range_with_index(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
        index: usize,
    ) -> Result<(), GRangesError> {
        // push an unindexed (empty) range
        let range = RangeIndexed::new(start, end, index);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<'a> GRangesEmpty<VecRangesEmpty> {
    /// Conduct a left overlap join, consuming self and returning a new
    /// [`GRanges<VecRangesIndexed, JoinData>>`].
    ///
    /// The [`JoinData`] container contains references to both left and right
    /// data containers and a [`Vec<OverlapJoin>`]. Each [`OverlapJoin`] represents
    /// a summary of an overlap, which downstream operations use to calculate
    /// statistics using the information about overlaps.
    pub fn left_overlaps<M: Clone + 'a, DR: 'a>(
        self,
        right: &'a impl AsGRangesRef<'a, COITrees<M>, DR>,
    ) -> Result<GRanges<VecRanges<RangeIndexed>, JoinData<'a, (), DR>>, GRangesError>
    where
        IntervalNode<M, usize>: GenericRange,
    {
        let mut gr: GRanges<VecRangesIndexed, JoinData<'a, (), DR>> =
            GRanges::new_vec(&self.seqlens());

        let right_ref = right.as_granges_ref();
        gr.data = Some(JoinData::new(None, right_ref.data.as_ref()));

        for (seqname, left_ranges) in self.0.ranges.iter() {
            for left_range in left_ranges.iter_ranges() {
                // Left join: every left range gets a JoinData.
                let mut join_data = LeftGroupedJoin::new(&left_range);
                if let Some(right_ranges) = right_ref.ranges.get(seqname) {
                    right_ranges.query(left_range.start(), left_range.end(), |right_range| {
                        // NOTE: right_range is a coitrees::IntervalNode.
                        join_data.add_right(&left_range, right_range);
                    });
                }
                gr.push_range_with_join(seqname, left_range.start, left_range.end, join_data)?;
            }
        }
        Ok(gr)
    }
}

impl<U> GRanges<VecRangesIndexed, Vec<U>> {
    /// Create a new [`GRanges<VecRangesIndexed, Vec<U>>`] object from an iterator over
    /// [`GenomicRangeRecord<U>`] records.
    pub fn from_iter<I>(
        iter: I,
        seqlens: &IndexMap<String, Position>,
    ) -> Result<GRanges<VecRangesIndexed, Vec<U>>, GRangesError>
    where
        I: Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    {
        let mut gr = GRanges::new_vec(seqlens);
        for possible_entry in iter {
            let entry = possible_entry?;
            gr.push_range(&entry.seqname, entry.start, entry.end, entry.data)?;
        }
        Ok(gr)
    }
}

impl GRangesEmpty<VecRangesEmpty> {
    /// Create a new [`GRanges<VecRangesEmpty, Vec<U>>`] object from an iterator over
    /// [`GenomicRangeRecord`] records.
    pub fn from_iter<I>(
        iter: I,
        seqlens: &IndexMap<String, Position>,
    ) -> Result<GRangesEmpty<VecRangesEmpty>, GRangesError>
    where
        I: Iterator<Item = Result<GenomicRangeEmptyRecord, GRangesError>>,
    {
        let mut gr = GRangesEmpty::new_vec(seqlens);
        for possible_entry in iter {
            let entry = possible_entry?;
            gr.push_range(&entry.seqname, entry.start, entry.end)?;
        }
        Ok(gr)
    }
}

impl<C> GRangesEmpty<C>
where
    COITrees<()>: From<C>,
{
    /// Convert the [`VecRangesEmpty`] range containers in this [`GRangesEmpty`] to a
    /// cache-oblivious interval tree range container, [`COITreesEmpty`]. This is
    /// done using the [`coitrees`] library by Daniel C. Jones.
    pub fn into_coitrees(self) -> Result<GRangesEmpty<COITreesEmpty>, GRangesError> {
        let old_ranges = self.0.ranges;
        let mut new_ranges = GenomeMap::new();
        for (seqname, vec_ranges) in old_ranges.into_iter() {
            let trees = COITrees::from(vec_ranges);
            new_ranges.insert(&seqname, trees)?;
        }
        Ok(GRangesEmpty(GRanges {
            ranges: new_ranges,
            data: None,
        }))
    }
}

impl<T> GRanges<VecRanges<RangeIndexed>, T> {
    /// Convert the [`VecRangesIndexed`] range containers in this [`GRanges`] to a
    /// cache-oblivious interval tree range container, [`COITreesIndexed`]. This is
    /// done using the [`coitrees`] library by Daniel C. Jones.
    pub fn into_coitrees(self) -> Result<GRanges<COITreesIndexed, T>, GRangesError> {
        let old_ranges = self.ranges;
        let mut new_ranges = GenomeMap::new();
        for (seqname, vec_ranges) in old_ranges.into_iter() {
            let trees = COITrees::from(vec_ranges);
            new_ranges.insert(&seqname, trees)?;
        }
        Ok(GRanges {
            ranges: new_ranges,
            data: self.data,
        })
    }
}

impl<CL: RangeContainer> GRangesEmpty<CL>
where
    CL: IterableRangeContainer,
{
    /// Filter out ranges that do *not* have at least overlap with the `right` ranges.
    ///
    /// In database lingo, this is a type of *filtering join*, in particular a *semi join*.
    /// See Hadley Wickham's excellent [R for Data
    /// Science](https://r4ds.hadley.nz/joins.html#filtering-joins) for more information.
    ///
    /// Note that this consumes the `self` [`GRanges`] object, turning it into a new
    /// [`GRanges<VecRangesEmpty, ()>`].
    ///
    pub fn filter_overlaps<'a, M: Clone + 'a, DR: 'a>(
        self,
        // right: &GRanges<COITrees<M>, DR>,
        right: &'a impl AsGRangesRef<'a, COITrees<M>, DR>,
    ) -> Result<GRangesEmpty<VecRangesEmpty>, GRangesError> {
        let mut gr = GRangesEmpty::new_vec(&self.seqlens());

        let right_ref = right.as_granges_ref();

        for (seqname, left_ranges) in self.0.ranges.iter() {
            for left_range in left_ranges.iter_ranges() {
                if let Some(right_ranges) = right_ref.ranges.get(seqname) {
                    let num_overlaps =
                        right_ranges.count_overlaps(left_range.start(), left_range.end());
                    if num_overlaps == 0 {
                        // no overlaps -- skip
                    } else {
                        gr.push_range(seqname, left_range.start(), left_range.end())?;
                    }
                }
            }
        }
        Ok(gr)
    }
}

impl<CL, U> GRanges<CL, Vec<U>>
where
    CL: IterableRangeContainer,
{
    /// Filter out ranges that do *not* have at least overlap with the `right` ranges.
    ///
    /// In database lingo, this is a type of *filtering join*, in particular a *semi join*.
    /// See Hadley Wickham's excellent [R for Data
    /// Science](https://r4ds.hadley.nz/joins.html#filtering-joins) for more information.
    ///
    /// Note that this consumes the `self` [`GRanges`] object, turning it into a new
    /// [`GRanges<VecRangesIndexed, Vec<U>`]. The data container is rebuilt from indices
    /// into a new [`Vec<U>`] where `U` is the associated type [`IndexedDataContainer::Item`],
    /// which represents the individual data element in the data container.
    pub fn filter_overlaps<'a, M: Clone + 'a, DR: 'a>(
        mut self,
        right: &'a impl AsGRangesRef<'a, COITrees<M>, DR>,
    ) -> Result<GRanges<VecRangesIndexed, Vec<U>>, GRangesError> {
        let mut gr: GRanges<VecRangesIndexed, Vec<U>> = GRanges::new_vec(&self.seqlens());

        let right_ref = right.as_granges_ref();
        let data = std::mem::take(&mut self.data).ok_or(GRangesError::NoDataContainer)?;

        let mut old_indices = HashSet::new(); // the old indices to *keep*
        let mut new_indices = Vec::new();
        let mut current_index = 0;

        for (seqname, left_ranges) in self.ranges.iter() {
            for left_range in left_ranges.iter_ranges() {
                if let Some(right_ranges) = right_ref.ranges.get(seqname) {
                    let num_overlaps =
                        right_ranges.count_overlaps(left_range.start(), left_range.end());
                    if num_overlaps == 0 {
                        // no overlaps -- skip
                    } else {
                        gr.push_range_with_index(
                            seqname,
                            left_range.start(),
                            left_range.end(),
                            current_index,
                        )?;
                        // unwrap should be safe, since this is an indexed GRanges
                        old_indices.insert(left_range.index().unwrap());
                        new_indices.push(current_index);
                        current_index += 1;
                    }
                }
            }
        }

        // Now, we reconstruct the right data. Note that we do not use the
        // standard push_range() method here, which would double the memory
        // usage essentially.
        let mut new_data = Vec::new();
        for (old_index, data_value) in data.into_iter().enumerate() {
            if old_indices.contains(&old_index) {
                new_data.push(data_value)
            }
        }
        ensure_eq!(new_data.len(), new_indices.len());
        gr.data = Some(new_data);
        Ok(gr)
    }
}

impl<R, T> GRanges<R, T>
where
    R: IterableRangeContainer,
{
    /// Create a new [`GRangesIterator`] to iterate through all the ranges in this [`GRanges`] object.
    pub fn iter_ranges(&self) -> GRangesIterator<'_, R> {
        GRangesIterator::new(&self.ranges)
    }
}

impl<R> GRangesEmpty<R>
where
    R: IterableRangeContainer,
{
    /// Create a new [`GRangesIterator`] to iterate through all the ranges in this [`GRangesEmpty`] object.
    pub fn iter_ranges(&self) -> GRangesIterator<'_, R> {
        GRangesIterator::new(&self.0.ranges)
    }
}

/// [`PartialEq`] for [`GRanges`] objects.
///
/// This is a more powerful comparison operator than [`GRanges.is_equal_to()`], since it will first
/// convert one type of ranges to another.
///
/// A [`GRanges`] object is considered equal to another if:
/// 1. Same number or ranges (this is checked first for efficiency).
/// 2. Same metadata.
/// 3. Same ranges. This requires converting ranges from different range container map types to a
///    [`VecRanges`]. **Warning**: this is currently memory intensive.
///
/// # Developer Notes
/// This uses type conversion in the range container types; see
/// [`crate::granges::ranges::comparison.rs`].
impl<CL, CR, DL, DR> PartialEq<GRanges<CR, DR>> for GRanges<CL, DL>
where
    CL: IterableRangeContainer + PartialEq<CL>,
    CR: IterableRangeContainer + PartialEq<CR>,
    DL: PartialEq<DR>,
{
    fn eq(&self, other: &GRanges<CR, DR>) -> bool {
        if self.len() != other.len() {
            return false;
        }

        let data_matches = match (&self.data, &other.data) {
            (Some(self_data), Some(other_data)) => self_data == other_data,
            (None, None) => true,
            _ => false,
        };

        let self_ranges = self.iter_ranges();
        let other_ranges = other.iter_ranges();
        let ranges_eq = self_ranges.zip(other_ranges).all(|(r1, r2)| r1 == r2);
        ranges_eq && data_matches
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        prelude::*,
        test_utilities::{granges_test_case_01, granges_test_case_02, random_vecranges},
        Position,
    };

    #[test]
    fn test_new_vec() {
        let seqlens = seqlens! { "chr1" => 10};
        let mut gr = GRanges::new_vec(&seqlens);
        gr.push_range("chr1", 0, 10, 1.1).unwrap();
        assert_eq!(gr.len(), 1);
    }

    #[test]
    fn test_random_vecranges() {
        let vr = random_vecranges(100);
        assert_eq!(vr.len(), 100)
    }

    #[test]
    fn test_to_coitrees() {
        let gr_vec = granges_test_case_01();
        let gr = gr_vec.clone().into_coitrees().unwrap();
        assert_eq!(gr.len(), 5);
    }

    #[test]
    fn granges_filter_overlaps() {
        let seqlens = seqlens! { "chr1" => 10};

        // test 1: one range that overlaps only the first range
        let gr = granges_test_case_01();
        let mut gr_keep: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(&seqlens);
        gr_keep.push_range("chr1", 0, 1).unwrap();
        let gr_keep = gr_keep.into_coitrees().unwrap();

        let gr_filtered = gr.filter_overlaps(&gr_keep).unwrap();
        assert_eq!(gr_filtered.len(), 1);

        // test 2: one range that overlaps the first two ranges
        let gr = granges_test_case_01();
        let mut gr_keep: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(&seqlens);
        gr_keep.push_range("chr1", 0, 5).unwrap();
        let gr_keep = gr_keep.into_coitrees().unwrap();

        let gr_filtered = gr.filter_overlaps(&gr_keep).unwrap();
        assert_eq!(gr_filtered.len(), 2);

        // test 3: one range that overlaps no ranges
        let gr = granges_test_case_01();
        let mut gr_keep: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(&seqlens);
        gr_keep.push_range("chr1", 8, 9).unwrap();
        let gr_keep = gr_keep.into_coitrees().unwrap();

        let gr_filtered = gr.filter_overlaps(&gr_keep).unwrap();
        assert_eq!(gr_filtered.len(), 0);

        // test 4: ranges on two chromosomes: first overlaps two ranges, second
        // overlaps one
        let gr = granges_test_case_01();
        let seqlens = seqlens! { "chr1" => 10, "chr2" => 10 };
        let mut gr_keep: GRangesEmpty<VecRangesEmpty> = GRangesEmpty::new_vec(&seqlens);
        gr_keep.push_range("chr1", 4, 7).unwrap();
        gr_keep.push_range("chr2", 10, 12).unwrap();
        let gr_keep = gr_keep.into_coitrees().unwrap();

        let gr_filtered = gr.filter_overlaps(&gr_keep).unwrap();
        assert_eq!(gr_filtered.len(), 3);
    }

    #[test]
    fn test_flanking_left() {
        let gr = granges_test_case_02();
        let gr_left = gr.flanking_ranges(Some(10), None).unwrap();

        let mut gr_left_iter = gr_left.iter_ranges();
        let first_range = gr_left_iter.next().unwrap();
        assert_eq!(first_range.start(), 20);
        assert_eq!(first_range.end(), 30);

        let second_range = gr_left_iter.next().unwrap();
        assert_eq!(second_range.start(), 90);
        assert_eq!(second_range.end(), 100);
    }

    #[test]
    fn test_flanking_both() {
        // Now with right flanks too.
        let gr = granges_test_case_02();
        let gr_left = gr.flanking_ranges(Some(10), Some(10)).unwrap();

        // First range is the new left flank.
        let mut gr_left_iter = gr_left.iter_ranges();
        let first_range = gr_left_iter.next().unwrap();
        assert_eq!(first_range.start(), 20);
        assert_eq!(first_range.end(), 30);

        // NOTE: the next flank would have been a zero-width
        // right range for the first range. But this is skipped,
        // because it is zero-width. Subtle edge case!

        // Now the next should be the new right flank.
        let second_range = gr_left_iter.next().unwrap();
        assert_eq!(second_range.start(), 90);
        let second_range = gr_left_iter.next().unwrap();
        assert_eq!(second_range.end(), 210);
    }

    #[test]
    fn test_from_windows() {
        let sl = seqlens!( "chr1" => 35);
        // Test chop
        let gr = GRangesEmpty::from_windows(&sl, 10, None, true).unwrap();
        assert_eq!(gr.len(), 3, "{:?}", gr);
        dbg!(&gr);
        gr.iter_ranges().for_each(|r| assert_eq!(r.width(), 10));

        // Test no chop
        let gr = GRangesEmpty::from_windows(&sl, 10, None, false).unwrap();
        assert_eq!(gr.iter_ranges().last().unwrap().width(), 5);

        // Test sliding with step of 2, with chop
        let sl = seqlens!( "chr1" => 21);
        let gr = GRangesEmpty::from_windows(&sl, 10, Some(2), true).unwrap();
        let mut expected_ranges_chop: Vec<(String, Position, Position)> = vec![
            ("chr1", 0, 10),
            ("chr1", 2, 12),
            ("chr1", 4, 14),
            ("chr1", 6, 16),
            ("chr1", 8, 18),
            ("chr1", 10, 20),
        ]
        .into_iter()
        .map(|(seq, s, e)| (seq.to_string(), s, e))
        .collect();
        let seqnames = sl.keys().map(|x| x.to_string()).collect::<Vec<_>>();
        let actual_ranges: Vec<(String, Position, Position)> = gr
            .iter_ranges()
            .map(|r| (r.seqname(&seqnames), r.start(), r.end()))
            .collect();
        assert_eq!(actual_ranges, expected_ranges_chop);

        // The same as above, without chop -- goes to seqlen with little remainder.
        let gr = GRangesEmpty::from_windows(&sl, 10, Some(2), false).unwrap();
        expected_ranges_chop.push(("chr1".to_string(), 12, 21));
        let actual_ranges: Vec<(String, Position, Position)> = gr
            .iter_ranges()
            .map(|r| (r.seqname(&seqnames), r.start(), r.end()))
            .collect();

        assert_eq!(actual_ranges, expected_ranges_chop);
    }

    #[test]
    fn test_left_overlaps() {
        let sl = seqlens!("chr1" => 50);
        let windows: GRangesEmpty<VecRangesEmpty> =
            GRangesEmpty::from_windows(&sl, 10, None, true).unwrap();

        let windows_len = windows.len();

        let mut right_gr: GRanges<VecRangesIndexed, Vec<f64>> = GRanges::new_vec(&sl);
        right_gr.push_range("chr1", 1, 2, 1.1).unwrap();
        right_gr.push_range("chr1", 5, 7, 2.8).unwrap();
        right_gr.push_range("chr1", 21, 35, 1.2).unwrap();
        right_gr.push_range("chr1", 23, 24, 2.9).unwrap();
        let right_gr = right_gr.into_coitrees().unwrap();

        let joined_results = windows.left_overlaps(&right_gr).unwrap();

        // get join data
        let data = joined_results.data.unwrap();

        // check is left join
        assert_eq!(data.len(), windows_len);

        let mut join_iter = data.iter();
        assert_eq!(join_iter.next().unwrap().num_overlaps(), 2);
        assert_eq!(join_iter.next().unwrap().num_overlaps(), 0);
        assert_eq!(join_iter.next().unwrap().num_overlaps(), 2);
        assert_eq!(join_iter.next().unwrap().num_overlaps(), 1);
        // rest are empty TODO should check
    }

    #[test]
    fn test_partial_eq() {
        // check equality case
        let gr1 = granges_test_case_01();
        let gr2 = granges_test_case_01();
        assert_eq!(gr1, gr2);

        // check minor difference case
        let sl = seqlens!( "chr1" => 30, "chr2" => 100 );
        let mut gr3 = GRanges::<VecRangesIndexed, Vec<f64>>::new_vec(&sl);
        gr3.push_range("chr1", 0, 4, 0.1).unwrap(); // one bp diff end
        gr3.push_range("chr1", 4, 7, 8.1).unwrap();
        gr3.push_range("chr1", 10, 17, 10.1).unwrap();
        gr3.push_range("chr2", 10, 20, 3.7).unwrap();
        gr3.push_range("chr2", 18, 32, 1.1).unwrap();

        assert_ne!(gr1, gr3);

        // check differential order case
        let sl = seqlens!( "chr1" => 30, "chr2" => 100 );
        let mut gr4 = GRanges::<VecRangesIndexed, Vec<f64>>::new_vec(&sl);
        gr4.push_range("chr1", 4, 7, 8.1).unwrap();
        gr4.push_range("chr1", 0, 5, 1.1).unwrap(); // swapped with above
        gr4.push_range("chr1", 10, 17, 10.1).unwrap();
        gr4.push_range("chr2", 10, 20, 3.7).unwrap();
        gr4.push_range("chr2", 18, 32, 1.1).unwrap();

        assert_ne!(gr3, gr4);

        // check manual case (for typos)
        let sl = seqlens!( "chr1" => 30, "chr2" => 100 );
        let mut gr5 = GRanges::<VecRangesIndexed, Vec<f64>>::new_vec(&sl);
        gr5.push_range("chr1", 0, 5, 1.1).unwrap();
        gr5.push_range("chr1", 4, 7, 8.1).unwrap();
        gr5.push_range("chr1", 10, 17, 10.1).unwrap();
        gr5.push_range("chr2", 10, 20, 3.7).unwrap();
        gr5.push_range("chr2", 18, 32, 1.1).unwrap();

        assert_eq!(gr1, gr5)
    }

    #[test]
    fn test_is_equal_to() {
        let vec_orig = granges_test_case_01();
        assert_eq!(vec_orig, vec_orig);

        let mut vec = vec_orig.clone();
        vec.push_range("chr1", 0, 10, 3.4).unwrap();
        assert_ne!(vec_orig, vec);

        let vec = vec.sort();
        let coit = vec.clone().into_coitrees().unwrap();
        assert_eq!(coit, vec);
    }
}
