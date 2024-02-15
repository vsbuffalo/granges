//! Functionality for general parsing, by turning BED-like files into iterators.
//!

use noodles::bed::{self};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use crate::Position;
use crate::error::GRangesError;
use crate::io::io::InputFile;
use crate::ranges::RangeRecord;
use crate::traits::GeneralRangeRecordIterator;

use super::noodles::convert_noodles_position;

pub struct Bed3RecordIterator {
    bed_reader: bed::Reader<BufReader<Box<dyn Read>>>,
}

impl Bed3RecordIterator {
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let input_file = InputFile::new(filepath);
        let reader = input_file.continue_reading()?;
        let bed_reader = bed::Reader::new(reader);
        Ok(Self { bed_reader })
    }

    pub fn from_reader(reader: BufReader<Box<dyn std::io::Read>>) -> Self {
        let bed_reader = bed::Reader::new(reader);
        Self { bed_reader }
    }
}

impl Iterator for Bed3RecordIterator {
    type Item = Result<RangeRecord<()>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.bed_reader.records::<3>().next().map(|res| {
            res.map_err(GRangesError::IOError)
                .map(|record| RangeRecord {
                    seqname: record.reference_sequence_name().to_string(),
                    start: convert_noodles_position(record.start_position()),
                    end: convert_noodles_position(record.end_position()),
                    data: (),
                })
        })
    }
}

impl GeneralRangeRecordIterator<()> for Bed3RecordIterator {
    fn retain_seqnames(self, seqnames: Vec<String>) -> FilteredIntervals<Self, ()> {
        FilteredIntervals::new(self, Some(seqnames), None)
    }

    fn exclude_seqnames(self, seqnames: Vec<String>) -> FilteredIntervals<Self, ()> {
        FilteredIntervals::new(self, None, Some(seqnames))
    }
}

pub struct TsvRecordIterator<F, U> {
    reader: BufReader<Box<dyn std::io::Read>>,
    parser: F,
    phantom: std::marker::PhantomData<U>,
}

impl<F, U> TsvRecordIterator<F, U>
where
    F: Fn(&str) -> Result<RangeRecord<U>, GRangesError>,
{
    pub fn new(filepath: impl Into<PathBuf>, parser: F) -> Result<Self, GRangesError> {
        let input_file = InputFile::new(filepath);
        let reader = input_file.continue_reading()?;

        Ok(Self {
            reader,
            parser,
            phantom: std::marker::PhantomData,
        })
    }
}

impl<F, U> Iterator for TsvRecordIterator<F, U>
where
    F: Fn(&str) -> Result<RangeRecord<U>, GRangesError>,
{
    type Item = Result<RangeRecord<U>, GRangesError>;

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

/// A BED-like file parser. This works by lazy-parsing the first three
/// columns, which are standard to all BED files.
pub struct BedlikeIterator {
    iter: TsvRecordIterator<fn(&str) -> Result<RangeRecord<String>, GRangesError>, String>,
}

impl BedlikeIterator {
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        // Wrap the parse_bedlike_to_range_record function to conform with TsvRecordIterator's expectations.
        let parser: fn(&str) -> Result<RangeRecord<String>, GRangesError> = parse_bed_lazy;
        
        let iter = TsvRecordIterator::new(filepath, parser)?;
        Ok(Self { iter })
    }
}

impl Iterator for BedlikeIterator {
    type Item = Result<RangeRecord<String>, GRangesError>;

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
/// use noodles::bed;
///
/// let iter = Bed3RecordIterator::new("tests_data/example.bed")
///            .expect("error reading file")
///            .exclude_seqnames(vec!["chr1".to_string()]);
///
/// let seqlens = seqlens! { "chr1" => 22, "chr2" => 10, "chr3" => 10, "chr4" => 15 };
/// let gr = GRanges::from_iter(iter, seqlens)
///            .expect("parsing error");
/// let mut iter = gr.iter_ranges();
///
/// // the first range should be the third range in the file,
/// // chr2:4-5
/// assert_eq!(iter.next().unwrap().start, 4);
/// assert_eq!(iter.next().unwrap().end, 5);
/// ```
pub struct FilteredIntervals<I, U>
where
    I: Iterator<Item = Result<RangeRecord<U>, GRangesError>>,
{
    inner: I,
    retain_seqnames: Option<HashSet<String>>,
    exclude_seqnames: Option<HashSet<String>>,
}

impl<I, U> FilteredIntervals<I, U>
where
    I: Iterator<Item = Result<RangeRecord<U>, GRangesError>>,
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

impl<I, U> Iterator for FilteredIntervals<I, U>
where
    I: Iterator<Item = Result<RangeRecord<U>, GRangesError>>,
{
    type Item = Result<RangeRecord<U>, GRangesError>;

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
pub fn parse_bed_lazy(line: &str) -> Result<RangeRecord<String>, GRangesError> {
    let columns: Vec<&str> = line.splitn(4, '\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::BedlikeTooFewColumns(line.to_string()));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    let data = columns[3].to_string();

    Ok(RangeRecord {
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
        let iter = Bed3RecordIterator::new("tests_data/noodles_example.bed").unwrap();

        let seqlens = seqlens! { "sq0" => 10 };
        let gr: GRanges<VecRangesEmpty, _> = GRanges::from_iter_empty(iter, seqlens).unwrap();

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
    #[should_panic]
    fn test_invalid_bed_noodles() {
        let iter = Bed3RecordIterator::new("tests_data/invalid.bed").unwrap();

        let seqlens = seqlens! { "sq0" => 10 };
        let _gr: GRanges<VecRangesEmpty, _> = GRanges::from_iter_empty(iter, seqlens).unwrap();
    }

    #[test]
    fn test_invalid_bedlike_iterator() {
        let iter = BedlikeIterator::new("tests_data/invalid.bed").unwrap();
        let seqlens = seqlens! { "sq0" => 10 };
        let result: Result<GRanges<VecRangesIndexed, _>, _> = GRanges::from_iter(iter, seqlens);
        dbg!(&result);
        // note: the Rust LSP thinks this isn't used for some reason, so prefaced with _
        // to silence warnings.
        let _msg = "column '-1' in 'chr1\t-1\t20'".to_string();
        assert!(matches!(
                result,
                Err(GRangesError::InvalidColumnType(_msg))
                ));
    }
}
