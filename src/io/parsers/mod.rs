//! Various types of parsing iterators for range formats.
//!
//! To work with genomic data in GRanges, one first needs to read it off disk and parse it. This is
//! done with *parsing iterators*. These iterators work on both plaintext and gzip-compressed
//! files (and bgzip and binary format support will hopefully be added). Each row of a file is
//! yielded as a particular parsed range type, which then can be filtered or altered in while
//! in the iterator, using Rust's powerful [`Iterator`] trait methods.
//!
//! Under the hood, this used the [`csv`] crate with [`serde`]. Initial parsers were handrolled,
//! but about 20%-30% less performant.
//!
//! ## Parsing Iterator Item Types
//!
//! These are the generic types that the iterator yields for each call of [`Iterator::next()`],
//! which represent a single row of data.
//!
//! There are two key item types that flow through parsing iterators, each representing a row of a
//! different kind of data file.
//!
//!  1. [`GenomicRangeEmptyRecord`], when a BED3 is loaded. This type indicates the incoming ranges
//!     *do not* have associated data that would go into a data container. It is *empty* of data.
//!
//!  2. [`GenomicRangeRecord<Option<String>>`], which is a range record *possibly* (hence,
//!     the [`Option`]) with remaining *unparsed* `String` data (e.g. the remaining unparsed
//!     columns of a VCF file).
//!
//!  3. [`GenomicRangeRecord<U>`], which is a range record with an *addition* set of data.
//!     For example, the standard BED5 has a two column [`Bed5Addition`], over the standard
//!     first three range columns. Thus, the [`Bed5Iterator`] yields
//!     [`GenomicRangeRecord<Bed5Addition>`] types.
//!
//! ## Why Empty Types?
//!
//! Most genomic file formats can be thought of as a genomic range with some associated data,
//! for example:
//!
//!  - VCF files will have information on the alleles and a vector of genotypes.
//!
//!  - GTF/GFF files store details about the annotated feature (e.g. name, protein ID, etc).
//!
//! However, often empty ranges (those that do not carry any data) are useful: for example, ranges
//! that just define genomic windows at regular widths, or masked repeat ranges. In this case, when
//! data is loaded into a [`GRanges<R, T>`] type, the data container type `T` is `()`. This is a
//! degenerate case that GRanges handles with a specific type called [`GRangesEmpty`], which indicates
//! the ranges in the ranges container in [`GRangesEmpty`] do not have any associated data.
//!
//! In most cases, people working with genomic file formats will know its type at compile-time.
//! However, if this isn't the case, see [⚙️  Handling Range Data when the type isn't know at runtime] below.
//!
//! # ⚙️  Handling Range Data when the type isn't know at runtime
//!
//! Because GRanges is a compile-time library, the type of incoming data is not runtime-polymorphic.
//! GRanges implements file detection methods to infer the filetype from the extension and
//! does additional type validation checks. When you want to use GRanges to handle a type
//! not known at compile-time, you will need to use `match` statements against the parsing
//! iterator types defined in [`GenomicRangesParser`]. The best examples of this are the
//! `granges` subcommands implementations.
//!
//! [`GRanges`]: crate::granges::GRanges
//! [`GRanges<R, T>`]: crate::granges::GRanges
//! [`GRangesEmpty`]: crate::granges::GRangesEmpty
//!

pub mod bed;
pub mod tsv;
pub mod utils;

pub use bed::{Bed3Iterator, Bed5Addition, Bed5Iterator, BedlikeIterator};

use std::collections::HashSet;
use std::io::BufRead;
use std::path::PathBuf;

use crate::error::GRangesError;
use crate::io::file::InputStream;
use crate::ranges::{GenomicRangeEmptyRecord, GenomicRangeRecord};
use crate::traits::{GeneralRangeRecordIterator, GenomicRangeRecordUnwrappable};
use crate::Position;

use self::utils::get_base_extension;

use super::TsvRecordIterator;

/// Inspect the first line to check that it looks like a valid BED-like
/// file, i.e. the first column is there (there are no reasonable checks
/// for sequence names other than presence), and the next to columns can
/// be parsed into a [`Position`].
fn valid_bedlike(input_file: &mut InputStream) -> Result<bool, GRangesError> {
    let _metadata = input_file.collect_metadata("#", None)?;
    let mut reader = input_file.continue_reading()?;
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    let columns = first_line
        .splitn(4, '\t')
        .map(String::from)
        .collect::<Vec<String>>();

    if columns.len() < 3 {
        // too few columns to be BED-like
        return Ok(false);
    }

    // Attempt to parse the second and third columns as positions
    let start_result = columns[1].trim().parse::<Position>();
    let end_result = columns[2].trim().parse::<Position>();

    // Check if both positions are valid
    match (start_result, end_result) {
        (Ok(_), Ok(_)) => Ok(true), // both are valid
        _ => Ok(false),             // one or both is not valid
    }
}

// TODO: we can combine this and the above into an enum
fn valid_bed5(input_file: &mut InputStream) -> Result<bool, GRangesError> {
    let _metadata = input_file.collect_metadata("#", None)?;
    let mut reader = input_file.continue_reading()?;
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    let columns = first_line
        .splitn(4, '\t')
        .map(String::from)
        .collect::<Vec<String>>();

    if columns.len() == 5 {
        // Attempt to parse the second and third columns as positions
        let valid_start = columns[1].trim().parse::<Position>().is_ok();
        let valid_end = columns[2].trim().parse::<Position>().is_ok();
        // this should always be true
        let _ = columns[3].trim().parse::<Position>().is_ok();
        let valid_strand = columns[4].trim().parse::<Position>().is_ok();
        return Ok(valid_start && valid_end && valid_strand);
    }
    Ok(false)
}

/// Enum that connects a genomic ranges file type to its specific parser.
#[derive(Debug)]
pub enum GenomicRangesParser {
    Bed3(Bed3Iterator),
    Bed5(Bed5Iterator),
    Bedlike(BedlikeIterator),
    Unsupported,
}

/// Enum that indicates the filetype of some genomic ranges file.
#[derive(Debug, PartialEq)]
pub enum GenomicRangesFile {
    Bed3(PathBuf),
    Bed5(PathBuf),
    Bedlike(PathBuf),
    Unsupported,
}

impl GenomicRangesFile {
    /// Detect the type of range genomic range file type we are working with, and output
    /// the appropriate [`GenomicRangesFile`] enum variant.
    ///
    /// Detection works like this:
    ///  1. Skip comment lines, starting with `#`.
    ///  2. Retrieve extension, removing any additional compression-related
    ///     extensions (`.gz` and `.bgz`) if present.
    ///  3. Read the first line, and try to parse the second and third columns
    ///     into [`Position`] types, to check if this is a BED-like file (i.e.
    ///     the first three columns are sequence name, start, and end positions).
    ///  4. Match against
    ///
    /// Currently this supports:
    ///  1. BED3 files, without range data.
    ///  2. BED-like files, with BED3 first columns and optional additional columns
    ///     after that. If there is no additional columnar data after the first
    ///     three BED3 columns, the type is [`GenomicRangesFile::Bed3`]. If there
    ///     is additional columnar data, it is [`GenomicRangesFile::Bedlike`].
    ///  3. Files with `.tsv` extension but are BED-like (first three
    ///     columns are BED3) will have the type [`GenomicRangesFile::Bed3`]. This
    ///     is because downstream [`GRanges`] operations need to know if any
    ///     additional data is present, which would need to be put in a data container.
    ///  4. BED5 files, which are BED3 + a *feature name* and a *strand* column.
    ///  5. If the file type does not satisfy any of the rules above, it is
    ///     [`GenomicRangesFile::Unsupported`].
    ///
    /// See the `match` statement in the source code for the exact rules. Additional
    /// genomic range file formats like GTF/GFF can easily be added later.
    ///
    /// [`GRanges`]: crate::granges::GRanges
    pub fn detect(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let path: PathBuf = filepath.into();
        let mut input_file = InputStream::new(&path);
        let _metadata = input_file.collect_metadata("#", None)?;
        let number_columns = input_file.detect_columns("\t")?;

        // get the extension, as a hint
        let extension =
            get_base_extension(&path).ok_or(GRangesError::CouldNotDetectRangesFiletype)?;

        // test if the first row can be parsed into a BED-like file format.
        let is_valid_bedlike = valid_bedlike(&mut input_file)?;
        let is_valid_bed5 = valid_bed5(&mut input_file)?; // TODO OPTIMIZE merge with above
        let file_type = match (
            extension.as_str(),
            number_columns,
            is_valid_bedlike,
            is_valid_bed5,
        ) {
            ("bed", 3, true, false) => GenomicRangesFile::Bed3(path),
            ("tsv", 3, true, false) => GenomicRangesFile::Bed3(path),
            ("bed", 5, true, true) => GenomicRangesFile::Bed5(path),
            ("bed", n, true, false) if n > 3 => GenomicRangesFile::Bedlike(path),
            ("tsv", n, true, false) if n > 3 => GenomicRangesFile::Bedlike(path),
            _ => GenomicRangesFile::Unsupported,
        };
        Ok(file_type)
    }

    /// Detect the genomic range filetype and link it to its parsing iterator, or raise an error
    /// if the filetype is not supported.
    ///
    /// This returns a [`GenomicRangesParser`] enum, since the parsing iterator filetype
    /// cannot be known at compile time.
    pub fn parsing_iterator(
        filepath: impl Clone + Into<PathBuf>,
    ) -> Result<GenomicRangesParser, GRangesError> {
        let path = filepath.into();
        match Self::detect(path)? {
            GenomicRangesFile::Bed3(path) => {
                Ok(GenomicRangesParser::Bed3(Bed3Iterator::new(path)?))
            }
            GenomicRangesFile::Bed5(path) => {
                Ok(GenomicRangesParser::Bed5(Bed5Iterator::new(path)?))
            }
            GenomicRangesFile::Bedlike(path) => {
                Ok(GenomicRangesParser::Bedlike(BedlikeIterator::new(path)?))
            }
            GenomicRangesFile::Unsupported => Err(GRangesError::UnsupportedGenomicRangesFileFormat),
        }
    }
}

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
impl<I> Iterator for FilteredRanges<I, GenomicRangeEmptyRecord>
where
    I: Iterator<Item = Result<GenomicRangeEmptyRecord, GRangesError>>,
{
    type Item = Result<GenomicRangeEmptyRecord, GRangesError>;

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

impl GeneralRangeRecordIterator<GenomicRangeEmptyRecord> for Bed3Iterator {
    fn retain_seqnames(self, seqnames: &[String]) -> FilteredRanges<Self, GenomicRangeEmptyRecord> {
        FilteredRanges::new(self, Some(&seqnames.to_vec()), None)
    }
    fn exclude_seqnames(
        self,
        seqnames: &[String],
    ) -> FilteredRanges<Self, GenomicRangeEmptyRecord> {
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
    column
        .parse::<T>()
        .map_err(|_| GRangesError::InvalidColumnType {
            expected_type: std::any::type_name::<T>().to_string(), // Provides the expected type name
            found_value: column.to_string(),
            line: line.to_string(),
        })
}

/// Parses a BED-like TSV format line into its constituent components.
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
        return Err(GRangesError::BedTooFewColumns(
            columns.len(),
            3,
            line.to_string(),
        ));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    // Collect remaining columns
    let additional_columns = columns[3..].to_vec();

    Ok((seqname, start, end, additional_columns))
}

#[cfg(test)]
mod tests {
    use crate::io::{
        parsers::{utils::get_base_extension, valid_bedlike},
        InputStream,
    };

    use super::GenomicRangesFile;

    #[test]
    fn test_get_base_extension() {
        assert_eq!(get_base_extension("test.bed.gz").unwrap(), "bed");
        assert_eq!(get_base_extension("test.bed").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.bed.gz").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.bed").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.gff").unwrap(), "gff");
        assert_eq!(get_base_extension("some/path/test.gff.gz").unwrap(), "gff");
        assert_eq!(get_base_extension("test.gff.gz").unwrap(), "gff");
        assert_eq!(get_base_extension("test.gff").unwrap(), "gff");
        assert_eq!(get_base_extension("test"), None);
        assert_eq!(get_base_extension("foo/test"), None);
    }

    #[test]
    fn test_rangefiletype_detect() {
        let range_filetype = GenomicRangesFile::detect("tests_data/example.bed");
        assert!(matches!(
            range_filetype.unwrap(),
            GenomicRangesFile::Bed3(_)
        ));

        let range_filetype = GenomicRangesFile::detect("tests_data/example_bedlike.tsv");
        assert!(matches!(
            range_filetype.unwrap(),
            GenomicRangesFile::Bedlike(_)
        ));
    }

    #[test]
    fn test_valid_bedlike() {
        assert_eq!(
            valid_bedlike(&mut InputStream::new("tests_data/example.bed")).unwrap(),
            true
        );
        assert_eq!(
            valid_bedlike(&mut InputStream::new("tests_data/invalid_format.bed")).unwrap(),
            false
        );
    }

    //#[test]
    //fn test_parser() {
    //    // based on this example: https://docs.rs/noodles-bed/latest/noodles_bed/struct.Reader.html#method.records
    //    let iter = Bed3Iterator::new("tests_data/noodles_example.bed").unwrap();

    //    let seqlens = seqlens! { "sq0" => 10 };
    //    let gr: GRanges<VecRangesEmpty, _> = GRanges::from_iter(iter, &seqlens).unwrap();

    //    assert_eq!(gr.len(), 2);

    //    // let mut gr_iter = gr.iter_ranges();
    //    // let first_interval = gr_iter.next().unwrap();
    //    // assert_eq!(first_interval.first, 7);
    //    // assert_eq!(first_interval.last, 12);
    //    //
    //    // let second_interval = gr_iter.next().unwrap();
    //    // assert_eq!(second_interval.first, 20);
    //    // assert_eq!(second_interval.last, 33);
    //}

    //#[test]
    //fn test_invalid_bedlike_iterator() {
    //    let iter = BedlikeIterator::new("tests_data/invalid.bed").unwrap();
    //    let seqlens = seqlens! { "sq0" => 10 };
    //    let result: Result<GRanges<VecRangesEmpty, _>, _> = GRanges::from_iter(iter.drop_data(), &seqlens);

    //    // note: the Rust LSP thinks this isn't used for some reason, so prefaced with _
    //    // to silence warnings.
    //    let _msg = "column '-1' in 'chr1\t-1\t20'".to_string();
    //    assert!(matches!(result, Err(GRangesError::InvalidColumnType(_msg))));
    //}
}
