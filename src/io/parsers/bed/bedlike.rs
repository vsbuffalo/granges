//! Bed-like file lazy parsers.

use std::{
    io::{BufRead, BufReader},
    path::PathBuf,
};

use crate::{
    io::{
        parsers::{tsv::build_tsv_reader, utils::parse_column},
        InputStream,
    },
    ranges::GenomicRangeRecord,
    GRangesError, Position,
};

// for BedlikeIterator only, TODO needed?
pub const PARSE_CAPACITY: usize = 512;

/// A lazy parser for BED-like files.
/// yields [`GenomicRangeRecord<Option<Vec<String>>>`] entries. If the file is a BED3 file,
/// the data in the [`GenomicRangeRecord`] will be set to `None`, since there are no remaining
/// string columns to parse.
pub struct BedlikeIterator {
    reader: BufReader<Box<dyn std::io::Read>>,
    line_buffer: String,
}

impl std::fmt::Debug for BedlikeIterator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BedlikeIterator").finish_non_exhaustive()
    }
}

impl BedlikeIterator {
    /// Create a new lazy-parsing iterator over Bed-like TSV data. This parser
    /// assumes the first three columns are the sequence name, start (0-indexed and inclusive),
    /// and end (0-indeed and exclusive) positions.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let input_file = InputStream::new(filepath);
        // let _has_metadata = input_file.collect_metadata("#", None);
        // let reader = input_file.continue_reading()?;
        let reader = input_file.reader()?;
        let line_buffer = String::with_capacity(PARSE_CAPACITY);
        Ok(Self {
            reader,
            line_buffer,
        })
    }
}

impl Iterator for BedlikeIterator {
    type Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buffer.clear();

        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => return None,
                Ok(_) => {
                    if !self.line_buffer.starts_with('#') {
                        let line = self.line_buffer.trim_end();
                        return Some(parse_bed_lazy(line));
                    }
                    // skip the metadata/comment character
                }
                Err(e) => return Some(Err(GRangesError::IOError(e))),
            }
        }
    }
}

/// Lazily parses a BED* format line into the first three columns defining the range,
/// storing the rest as a `String`.
///
/// Unlike [`parse_bedlike()`], this function returns a concrete
/// [`GenomicRangeRecord<Option<String>>`] type, not a [`Vec<String>`] of split TSV columns.
pub fn parse_bed_lazy(line: &str) -> Result<GenomicRangeRecord<Option<String>>, GRangesError> {
    let columns: Vec<&str> = line.splitn(4, '\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::Bed3TooFewColumns(
            columns.len(),
            line.to_string(),
        ));
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

/// Inspect the first line to check that it looks like a valid BED-like
/// file, i.e. the first column is there (there are no reasonable checks
/// for sequence names other than presence), and the next to columns can
/// be parsed into a [`Position`].
// TODO/NOTE: this is an older-style parser. We might want to try
// using csv here.

pub fn valid_bedlike(filepath: impl Into<PathBuf>) -> Result<bool, GRangesError> {
    let filepath = filepath.into();
    let mut reader = build_tsv_reader(filepath)?;
    let mut records = reader.records();
    if let Some(result) = records.next() {
        let record = result?;

        if record.len() < 3 {
            // too few columns to be BED-like
            return Ok(false);
        }

        // Attempt to parse the second and third columns as positions
        let start_result = record.get(1).unwrap().trim().parse::<usize>();
        let end_result = record.get(2).unwrap().trim().parse::<usize>();

        // Check if both positions are valid
        match (start_result, end_result) {
            (Ok(_), Ok(_)) => Ok(true), // both are valid
            _ => Ok(false),             // one or both is not valid
        }
    } else {
        // No records found
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::valid_bedlike;

    #[test]
    fn test_valid_bedlike() {
        // even bed files work
        assert_eq!(valid_bedlike("tests_data/example.bed").unwrap(), true);

        // but not everything
        assert_eq!(
            valid_bedlike("tests_data/invalid_format.bed").unwrap(),
            false
        );

        assert_eq!(
            valid_bedlike("tests_data/example_bedlike.tsv").unwrap(),
            true
        );
    }
}
