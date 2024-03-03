//! Filters for parsing iterators.

use crate::error::GRangesError;
use crate::io::TsvRecordIterator;
use crate::ranges::{GenomicRangeRecord, GenomicRangeRecordEmpty};
use crate::traits::{GeneralRangeRecordIterator, GenomicRangeRecordUnwrappable};
use std::collections::HashSet;

use super::bed::{Bed4Addition, Bed4Iterator};
use super::{Bed3Iterator, Bed5Addition, Bed5Iterator, BedlikeIterator};

/// An iterator over a generic "genomic range like " item type `R`, that filters based on sequence name.
///
/// Note that that the exclude filter is prioritized over the retain filter. So,
/// if a pipeline contains both, then if a chromosome is supplied to
/// [`GeneralIntervalIterator.exclude_seqnames()`], then even if it is in retain,
/// it is skipped. This is to prevent validation and returning a [`Result`].
///
/// # Example
///
/// ```
/// use granges::prelude::*;
///
/// let iter = Bed3Iterator::new("tests_data/example.bed")
///            .expect("error reading file")
///            .exclude_seqnames(&vec!["chr1".to_string()]);
///
/// let seqlens = seqlens! { "chr1" => 22, "chr2" => 10, "chr3" => 10, "chr4" => 15 };
/// let gr = GRangesEmpty::from_iter(iter, &seqlens)
///            .expect("parsing error");
/// let mut iter = gr.iter_ranges();
///
/// // the first range should be the third range in the file,
/// // chr2:4-5
/// assert_eq!(iter.next().unwrap().start, 4);
/// assert_eq!(iter.next().unwrap().end, 6);
/// assert_eq!(iter.next().unwrap().end, 15);
/// assert_eq!(iter.next(), None);
/// ```
#[derive(Debug)]
pub struct FilteredRanges<I, R>
where
    I: Iterator<Item = Result<R, GRangesError>>,
{
    inner: I,
    retain_seqnames: Option<HashSet<String>>,
    exclude_seqnames: Option<HashSet<String>>,
}

impl<I, R> FilteredRanges<I, R>
where
    I: Iterator<Item = Result<R, GRangesError>>,
{
    pub fn new(
        inner: I,
        retain_seqnames: Option<&Vec<String>>,
        exclude_seqnames: Option<&Vec<String>>,
    ) -> Self {
        let retain_seqnames = retain_seqnames.cloned().map(HashSet::from_iter);
        let exclude_seqnames = exclude_seqnames.cloned().map(HashSet::from_iter);
        Self {
            inner,
            retain_seqnames,
            exclude_seqnames,
        }
    }
}

/// Range-filtering iterator implementation for [`GenomicRangeRecord<U>`].
impl<I, U> Iterator for FilteredRanges<I, GenomicRangeRecord<U>>
where
    I: Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
{
    type Item = Result<GenomicRangeRecord<U>, GRangesError>;

    /// Get the next filtered entry, prioritizing exclude over retain.
    fn next(&mut self) -> Option<Self::Item> {
        for item in self.inner.by_ref() {
            match &item {
                Ok(entry) => {
                    if self
                        .exclude_seqnames
                        .as_ref()
                        .map_or(false, |ex| ex.contains(&entry.seqname))
                    {
                        continue;
                    }
                    if self
                        .retain_seqnames
                        .as_ref()
                        .map_or(true, |rt| rt.contains(&entry.seqname))
                    {
                        return Some(item);
                    }
                }
                Err(_) => return Some(item),
            }
        }
        None
    }
}

/// Range-filtering iterator implementation for [`GenomicRangeEmptyRecord`].
impl<I> Iterator for FilteredRanges<I, GenomicRangeRecordEmpty>
where
    I: Iterator<Item = Result<GenomicRangeRecordEmpty, GRangesError>>,
{
    type Item = Result<GenomicRangeRecordEmpty, GRangesError>;

    /// Get the next filtered entry, prioritizing exclude over retain.
    fn next(&mut self) -> Option<Self::Item> {
        for item in self.inner.by_ref() {
            match &item {
                Ok(entry) => {
                    if self
                        .exclude_seqnames
                        .as_ref()
                        .map_or(false, |ex| ex.contains(&entry.seqname))
                    {
                        continue;
                    }
                    if self
                        .retain_seqnames
                        .as_ref()
                        .map_or(true, |rt| rt.contains(&entry.seqname))
                    {
                        return Some(item);
                    }
                }
                Err(_) => return Some(item),
            }
        }
        None
    }
}

impl<U> GeneralRangeRecordIterator<GenomicRangeRecord<U>>
    for TsvRecordIterator<GenomicRangeRecord<U>>
where
    U: Clone + for<'de> serde::Deserialize<'de>,
{
    fn retain_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, GenomicRangeRecord<U>> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, GenomicRangeRecord<U>> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

impl GeneralRangeRecordIterator<GenomicRangeRecord<Option<String>>> for BedlikeIterator {
    fn retain_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Option<String>>> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Option<String>>> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

impl GeneralRangeRecordIterator<GenomicRangeRecordEmpty> for Bed3Iterator {
    fn retain_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, GenomicRangeRecordEmpty> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecordEmpty> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

impl GeneralRangeRecordIterator<GenomicRangeRecord<Bed4Addition>> for Bed4Iterator {
    fn retain_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Bed4Addition>> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Bed4Addition>> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

impl GeneralRangeRecordIterator<GenomicRangeRecord<Bed5Addition>> for Bed5Iterator {
    fn retain_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Bed5Addition>> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<Bed5Addition>> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

impl<I> GeneralRangeRecordIterator<GenomicRangeRecord<String>> for UnwrappedRanges<I>
where
    I: Iterator<Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>>,
{
    fn retain_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<String>> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeRecord<String>> {
        FilteredRanges::new(self, None, Some(&seqnames.to_vec()))
    }
}

#[derive(Debug)]
pub struct UnwrappedRanges<I>
where
    I: Iterator<Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>>,
{
    inner: I,
}

impl<I> UnwrappedRanges<I>
where
    I: Iterator<Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>>,
{
    pub fn new(inner: I) -> Self {
        Self { inner }
    }
}

impl<I> Iterator for UnwrappedRanges<I>
where
    I: Iterator<Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>>,
{
    type Item = Result<GenomicRangeRecord<String>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(result) = self.inner.next() {
            return Some(result.and_then(|record| match record.data {
                Some(data) => Ok(GenomicRangeRecord::new(
                    record.seqname,
                    record.start,
                    record.end,
                    data,
                )),
                None => Err(GRangesError::TryUnwrapDataError),
            }));
        }
        None
    }
}

impl GenomicRangeRecordUnwrappable for BedlikeIterator {
    fn try_unwrap_data(self) -> UnwrappedRanges<Self> {
        UnwrappedRanges::new(self)
    }
}
