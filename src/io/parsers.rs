//! Functionality for general parsing, by turning BED-like files into iterators.
//!
//! ## Parsers
//!
//! Parsers in the GRanges library are designed to be:
//!
//!  1. *Iterator-based*: GRanges parsers are iterators, allow entries to be filtered or
//!     manipulated on the fly using Rust's iterator methods like [`Iterator.filter`].
//!
//!  2. *Lazy*: GRanges parsers are lazy. Since all BED file formats are a superset of BED3, we can
//!     parse the range components out of the first three columns and then just store the remaining
//!     unparsed part of the line in a `String`.
//!
//!  3. *Permissive*: All GRanges parsers are built off of the [`TsvRecordIterator`], which reads
//!     plaintext and gzip-compressed files and parses them according to a specified function. This
//!     general parser should be able to accomodate every line-based range bioinformatics format
//!     (BED3, BED5, GTF, GFF, VCF, etc). Downstream users can implement their own specific parsers
//!     for these variants (please feel free to contribute a parser for a commonly-used variant to
//!     GRanges too!).
//!
//! Since GRanges is fundamentally about manipulating genomic range data, all parsers output to the
//! same record type: [`GenomicRangeRecord<Option<String>>`], which is generic over the data. The
//! [`Option<String>`] here is because the lazy parser *at compile time* does not know how many
//! columns it will encounter. If only three columns are encountered (e.g. a BED3 file), the data
//! in this [`GenomicRangeRecord`] are all `None`. Then, the ranges that go into [`GRanges`] object
//! do not have indices, since there is not data container.
//!
//! Otherwise, if there *is* data, this data needs to be pushed to the data container, and the
//! indexed range type ([`RangeIndexed`]) is used in the range containers.
//!
//! While handling of this at compile time leads to very performant code with lower memory
//! overhead, it has the downside that *we must handle both types at compile time*.
//!
//! 
//! # Working Downstream of Parsers
//!
//! Often at runtime the exact file format may not be known. A user could specify a BED3 file,
//! which only contains ranges and no data, or a BED* of BED-like file (see terminology below).
//! Since GRanges aims to handle most situations *at compile time*, it must handle the process of
//! figuring out how to take an iterator of [`GenomicRangeRecord<U>`] entries and 
//!
//!
//!
//! ## Terminology 
//!
//!  - BED3 - BED* - BED-like
//!
//!
//! : it could be ranges-only (i.e. a BED3), or contain data (e.g. BED5). The lazy BED parser will
//! output a [`GenomicRangeRecord<Option<String>>`], where the data would be `None` only in the
//! case that three columns were encountered in the file (which must be a BED3). 
//!
//! In GRanges, there are two types of ranges: ranges with an index to an element in the data
//! container, and ranges without indices (i.e. what we would use in the case of processing a BED3
//! file). Since a [`GRanges`] object needs to have a single, concrete range type in its range
//! containers, it must be known *at compile time* how one should convert the 
//!
//! Downstream pipelines must immediately determine how to handle whether there is additional data,
//! or all of the [`GenomicRangeRecord`] entries are `None`, and 
//!
//!
//! # BED-like File Parser Design 
//!
//! All BED formats (BED3, BED5, etc). are built upon a BED3. Often when working with these types
//! of formats, many operations do not immediately require full parsing of the line, past the
//! *range components*. This is because downstream operations may immediately filter away entries
//! (e.g. based on width), or do an overlap operation, and then filter away entries based on some
//! overlap criteria. Either way, it may be advantageous to work with ranges that just store some
//! generic data.
//!

use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use indexmap::IndexMap;

use crate::error::GRangesError;
use crate::granges::GRanges;
use crate::io::file::InputFile;
use crate::ranges::GenomicRangeRecord;
use crate::Position;
use crate::traits::{RangeContainer, GenomicRangesOperationsExtended};

/// An extensible TSV parser, which uses a supplied parser function to
/// convert a line into a [`RangeRecord<U>`], a range with generic associated
/// data.
pub struct TsvRecordIterator<F, U> {
    reader: BufReader<Box<dyn std::io::Read>>,
    num_columns: usize,
    parser: F,
    phantom: std::marker::PhantomData<U>,
}

impl<F, U> TsvRecordIterator<F, U>
where
    F: Fn(&str) -> Result<GenomicRangeRecord<U>, GRangesError>,
{
    /// Create a new [`TsvRecordIterator`], which parses lines from the supplied
    /// file path into [`RangeRecord<U>`] using the specified parsing function.
    pub fn new(filepath: impl Into<PathBuf>, parser: F) -> Result<Self, GRangesError> {
        let mut input_file = InputFile::new(filepath);
        let _has_metadata = input_file.collect_metadata("#", None);
        let num_columns = input_file.detect_columns("\t")?;
        let reader = input_file.continue_reading()?;

        Ok(Self {
            reader,
            num_columns,
            parser,
            phantom: std::marker::PhantomData,
        })
    }
}

impl<F, U> Iterator for TsvRecordIterator<F, U>
where
    F: Fn(&str) -> Result<GenomicRangeRecord<U>, GRangesError>,
{
    type Item = Result<GenomicRangeRecord<U>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => None,
            Ok(_) => {
                let line = line.trim_end();
                Some((self.parser)(line))
            }
            Err(e) => Some(Err(GRangesError::IOError(e))),
        }
    }
}

/// A lazy parser for BED-like files. This lazily parses only the first three columns, and
/// yields [`GenomicRangeRecord<Option<String>>`] entries. If the file is a BED3 file,
/// the data in the [`GenomicRangeRecord`] will be set to `None`, since there are no remaining
/// string columns to parse.
///
#[allow(clippy::type_complexity)]
pub struct BedlikeIterator {
    iter: TsvRecordIterator<
        fn(&str) -> Result<GenomicRangeRecord<Option<String>>, GRangesError>,
        Option<String>,
    >,
}

pub enum GenomicRangesIteratorVariant {
    WithData(Box<dyn Iterator<Item = Result<GenomicRangeRecord<String>, GRangesError>>>),
    WithoutData(Box<dyn Iterator<Item = Result<GenomicRangeRecord<()>, GRangesError>>>),
}

impl BedlikeIterator {
    /// Create a new lazy-parsing iterator over Bed-like TSV data. This parser
    /// assumes the first three columns are the sequence name, start (0-indexed and inclusive),
    /// and end (0-indeed and exclusive) positions.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let parser: fn(&str) -> Result<GenomicRangeRecord<Option<String>>, GRangesError> =
            parse_bed_lazy;

        let iter = TsvRecordIterator::new(filepath, parser)?;
        Ok(Self { iter })
    }

    /// Detect the number of columns from the first entry.
    ///
    /// Note: this does not guard against the risk of ragged input, i.e. differing
    /// numbers of columns per row.
    pub fn number_columns(&self) -> usize {
        self.iter.num_columns
    }

    pub fn is_bed3(&self) -> bool {
        self.number_columns() == 3
    }

    /// Try to unwrap each [`GenomicRangeRecord<Option<String>>`] iterator item into a
    /// [`GenomicRangeRecord<String>`] iterator item.
    ///
    /// This will raise errors if:
    ///  1. The detected number of columns is < 4 (i.e. there appears to be no data to unwrap).
    ///  2. During iteration, a `None` data element is encountered.
    pub fn try_unwrap_data(
        self,
    ) -> Result<impl Iterator<Item = Result<GenomicRangeRecord<String>, GRangesError>>, GRangesError>
    {
        if self.number_columns() < 4 {
            return Err(GRangesError::TooFewColumns)?;
        }
        Ok(self.iter.map(|result| {
            result.and_then(|record| {
                if let Some(data) = record.data {
                    Ok(GenomicRangeRecord::new(
                        record.seqname,
                        record.start,
                        record.end,
                        data,
                    ))
                } else {
                    Err(GRangesError::TryUnwrapDataError)
                }
            })
        }))
    }

    /// Drop the data in each [`GenomicRangeRecord<Option<String>>`] iterator, converting it to a range-only
    /// [`GenomicRangeRecord<()>`] iterator item.
    pub fn drop_data(self) -> impl Iterator<Item = Result<GenomicRangeRecord<()>, GRangesError>> {
        self.iter.map(|result| {
            result
                .map(|record| {
                    Ok(GenomicRangeRecord::new(
                        record.seqname,
                        record.start,
                        record.end,
                        (),
                    ))
                })
                .unwrap_or_else(|e| Err(e)) // pass through parsing errors
        })
    }

    /// Map the iterator into [`GRanges`] object, by using the specified function to convert the
    /// lazily-parsed [`GenomicRangeRecord<Option<String>>`] into some concrete type.
    pub fn map_into_granges<C: RangeContainer, F, U>(self, func: F, seqlens: &IndexMap<String, Position>) 
        -> Result<GRanges<C, <GRanges<C, Vec<U>> as GenomicRangesOperationsExtended<C>>::DataContainerType>, GRangesError>
        where GRanges<C, Vec<U>>: GenomicRangesOperationsExtended<C>,
        F: Fn(Result<GenomicRangeRecord<Option<String>>, GRangesError>) -> 
            Result<GenomicRangeRecord<<GRanges<C, Vec<U>> as GenomicRangesOperationsExtended<C>>::DataElementType>, GRangesError> {
        GRanges::from_iter(self.iter.map(func), seqlens)
    }


    /// Try to convert the iterator into one of two variants: one with data and one without.
    pub fn into_variant(self) -> Result<GenomicRangesIteratorVariant, GRangesError> {
        let number_columns = self.number_columns();

        if number_columns == 3 {
            let without_data_iterator = self.drop_data().map(|result| {
                result.map(|record| GenomicRangeRecord {
                    seqname: record.seqname,
                    start: record.start,
                    end: record.end,
                    data: (),
                })
            });
            Ok(GenomicRangesIteratorVariant::WithoutData(Box::new(
                without_data_iterator,
            )))
        } else {
            let with_data_iterator = self.try_unwrap_data()?.map(|result| {
                result.map(|record| GenomicRangeRecord {
                    seqname: record.seqname,
                    start: record.start,
                    end: record.end,
                    data: record.data,
                })
            });
            Ok(GenomicRangesIteratorVariant::WithData(Box::new(
                with_data_iterator,
            )))
        }
    }
}

impl Iterator for BedlikeIterator {
    type Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// An iterator over [`IntervalRecord`] items that filters based on sequence name.
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
/// let iter = BedlikeIterator::new("tests_data/example.bed")
///            .expect("error reading file")
///            .exclude_seqnames(vec!["chr1".to_string()]);
///
/// let seqlens = seqlens! { "chr1" => 22, "chr2" => 10, "chr3" => 10, "chr4" => 15 };
/// let gr = GRanges::from_iter(iter, &seqlens)
///            .expect("parsing error");
/// let mut iter = gr.iter_ranges();
///
/// // the first range should be the third range in the file,
/// // chr2:4-5
/// assert_eq!(iter.next().unwrap().start, 4);
/// assert_eq!(iter.next().unwrap().end, 5);
/// ```
pub struct FilteredRanges<I, U>
where
    I: Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
{
    inner: I,
    retain_seqnames: Option<HashSet<String>>,
    exclude_seqnames: Option<HashSet<String>>,
}

impl<I, U> FilteredRanges<I, U>
where
    I: Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
{
    pub fn new(
        inner: I,
        retain_seqnames: Option<Vec<String>>,
        exclude_seqnames: Option<Vec<String>>,
    ) -> Self {
        let retain_seqnames = retain_seqnames.map(HashSet::from_iter);
        let exclude_seqnames = exclude_seqnames.map(HashSet::from_iter);
        Self {
            inner,
            retain_seqnames,
            exclude_seqnames,
        }
    }
}

impl<I, U> Iterator for FilteredRanges<I, U>
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

/// Parses a single column from a string slice into a specified type.
///
/// This function is particularly useful for converting columns in genomic data files
/// from string representation to a specific type like `usize`, `f64`, etc.
///
/// # Type Parameters
///
/// * `T`: The type to which the column should be parsed. It must implement `FromStr`.
///
/// # Arguments
///
/// * `column`: The string slice representing the column to be parsed.
/// * `line`: The entire line from which the column was extracted. Used for error reporting.
///
/// # Errors
///
/// Returns `GRangesError::InvalidColumnType` if the column cannot be parsed into type `T`.
pub fn parse_column<T: std::str::FromStr>(column: &str, line: &str) -> Result<T, GRangesError>
where
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    // NOTE: this is used a lot, and should be benchmarked.
    column
        .parse::<T>()
        .map_err(|_| GRangesError::InvalidColumnType(format!("column '{}' in '{}'", column, line)))
}

/// Parses a BED-like format line into its constituent components.
///
/// BED-like formats contain BED3 columns (chromosome, start, and end) and
/// others. It's an non-standard but extensible (since it's just a TSV), and sensible
/// for a variety of genomic data.
///
/// # Arguments
///
/// * `line`: A string slice representing a line in BED3 format.
///
/// # Errors
///
/// Returns `GRangesError::InvalidColumnType` if the line does not have enough columns
/// or if any column cannot be correctly parsed.
pub fn parse_bedlike(line: &str) -> Result<(String, Position, Position, Vec<&str>), GRangesError> {
    let columns: Vec<&str> = line.split('\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::BedlikeTooFewColumns(line.to_string()));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    // Collect remaining columns
    let additional_columns = columns[3..].to_vec();

    Ok((seqname, start, end, additional_columns))
}

/// Lazily parses a BED* format line into the first three columns defining the range,
/// storing the rest as a `String`.
pub fn parse_bed_lazy(line: &str) -> Result<GenomicRangeRecord<Option<String>>, GRangesError> {
    let columns: Vec<&str> = line.splitn(4, '\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::BedlikeTooFewColumns(line.to_string()));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    let data = if columns.len() > 3 {
        Some(columns[3].to_string())
    } else {
        None
    };

    Ok(GenomicRangeRecord {
        seqname,
        start,
        end,
        data,
    })
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    use super::BedlikeIterator;

    #[test]
    fn test_parser() {
        // based on this example: https://docs.rs/noodles-bed/latest/noodles_bed/struct.Reader.html#method.records
        let iter = BedlikeIterator::new("tests_data/noodles_example.bed").unwrap();

        let seqlens = seqlens! { "sq0" => 10 };
        let gr: GRanges<VecRangesEmpty, _> = GRanges::from_iter(iter.drop_data(), &seqlens).unwrap();

        assert_eq!(gr.len(), 2);

        // let mut gr_iter = gr.iter_ranges();
        // let first_interval = gr_iter.next().unwrap();
        // assert_eq!(first_interval.first, 7);
        // assert_eq!(first_interval.last, 12);
        //
        // let second_interval = gr_iter.next().unwrap();
        // assert_eq!(second_interval.first, 20);
        // assert_eq!(second_interval.last, 33);
    }

    #[test]
    fn test_invalid_bedlike_iterator() {
        let iter = BedlikeIterator::new("tests_data/invalid.bed").unwrap();
        let seqlens = seqlens! { "sq0" => 10 };
        let result: Result<GRanges<VecRangesEmpty, _>, _> = GRanges::from_iter(iter.drop_data(), &seqlens);

        // note: the Rust LSP thinks this isn't used for some reason, so prefaced with _
        // to silence warnings.
        let _msg = "column '-1' in 'chr1\t-1\t20'".to_string();
        assert!(matches!(result, Err(GRangesError::InvalidColumnType(_msg))));
    }
}
