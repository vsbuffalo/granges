//! The [`GRangesError`] `enum` definition and error messages.
//!
use crate::Position;
use genomap::GenomeMapError;
use std::{
    num::{ParseFloatError, ParseIntError},
    string::FromUtf8Error,
};
use thiserror::Error;

///// The [`GRangesError`] defines the standard set of errors that should
///// be passed to the user.
//#[derive(Debug, Error)]
//pub enum GRangesError {
//    // IO related errors
//    #[error("File reading eror: {0}")]
//    IOError(#[from] std::io::Error),
//
//    // File parsing related errors
//    #[error("Could not detect genomic ranges filetype from extension.")]
//    CouldNotDetectRangesFiletype,
//    #[error("Integer parsing error: {0}")]
//    ParseIntError(#[from] ParseIntError),
//    #[error("Float parsing error: {0}")]
//    ParseFloatError(#[from] ParseFloatError),
//    #[error("Bed-like file has too few columns. The first three columns must be sequence name, and start and end positions.\nLine: {0}")]
//    BedlikeTooFewColumns(String),
//    #[error("File has invalid column type entry: {0}")]
//    InvalidColumnType(String),
//    #[error("Genome file is invalid: {0}")]
//    InvalidGenomeFile(String),
//    #[error("Invalid BED string: must be either '+', '-', or '.'")]
//    InvalidString,
//
//    // BedlikeIterator errors
//    #[error("GenomicRangeRecord encountered with None data in try_unwrap_data()")]
//    TryUnwrapDataError,
//    #[error("GenomicRangeRecord.try_unwrap_data() called on TSV with fewer than 4 columns")]
//    TooFewColumns,
//
//    // Invalid genomic range errors
//    #[error("Range invalid: start ({0}) must be greater than end ({1})")]
//    InvalidGenomicRange(Position, Position),
//    #[error("Range [{0}, {1}] is invalid for sequence of length {2}")]
//    InvalidGenomicRangeForSequence(Position, Position, Position),
//    #[error("Sequence name '{0}' is not in the ranges container")]
//    MissingSequence(String),
//    #[error("Error encountered in genomap::GenomeMap")]
//    GenomeMapError(#[from] GenomeMapError),
//
//    #[error("Invalid GRanges object: no data container.")]
//    NoDataContainer,
//
//    // Sequences related errors
//    #[error("Sequence name '{0}' was not found.")]
//    MissingSequenceName(String),
//
//    // FASTA/noodles related errors
//    #[error("Error encountered in trying to convert bytes to UTF8 string.")]
//    FromUtf8Error(#[from] FromUtf8Error),
//
//    // Command line tool related errors
//    #[error("Unsupported genomic range format")]
//    UnsupportedGenomicRangesFileFormat,
//    #[error("Command line argument error: {0}")]
//    ArgumentError(#[from] clap::error::Error),
//    #[error("No such operation: {0}")]
//    NoSuchOperation(String),
//}

// AICODE/NOTE: these are AI-generated from the above human-written ones.
#[derive(Debug, Error)]
pub enum GRangesError {
    // IO related errors
    #[error("File reading error: {0}. Please check if the file exists and you have permission to read it.")]
    IOError(#[from] std::io::Error),

    #[error("File parsing error: {0}")]
    TsvParsingError(#[from] csv::Error),

    // File parsing related errors
    #[error("Could not determine the file type based on its extension. Ensure the file has a standard genomic data extension (.bed, .gff, etc.).")]
    CouldNotDetectRangesFiletype,

    #[error("Error parsing an integer value: {0}. Check the file for non-integer entries where integers are expected.")]
    ParseIntError(#[from] ParseIntError),

    #[error("Error parsing a floating-point value: {0}. Check the file for non-numeric entries where floating-point numbers are expected.")]
    ParseFloatError(#[from] ParseFloatError),

    #[error(
        "The provided BED file has fewer columns ({0}) than expected ({1}). Problematic line:\n{2}"
    )]
    BedTooFewColumns(usize, usize, String),

    #[error("The provided BED3 file has fewer columns ({0}) than expected (3).\nAt least three columns are needed: sequence name, start, and end positions.\nProblematic line:\n{1}")]
    Bed3TooFewColumns(usize, String),

    #[error(
        "Invalid column type: expected {expected_type} but got '{found_value}' in line: '{line}'."
    )]
    InvalidColumnType {
        expected_type: String,
        found_value: String,
        line: String,
    },

    #[error("The input file had no rows.")]
    NoRows,

    #[error("The genome file is invalid: {0}. Please verify the file's format and contents.")]
    InvalidGenomeFile(String),

    #[error("Invalid BED format detected. Each entry must be '+', '-', or '.' to represent strand information.")]
    InvalidString,

    // BedlikeIterator errors
    #[error("Attempted to unwrap genomic range data, but none was present. This operation requires data to be associated with each genomic range.")]
    TryUnwrapDataError,

    #[error("Attempted to unwrap genomic range data from a TSV with fewer than 4 columns. Ensure your TSV includes the necessary genomic range columns plus additional data.")]
    TooFewColumns,

    // Invalid genomic range errors
    #[error("Invalid genomic range specified: start position ({0}) must be less than or equal to the end position ({1}).")]
    InvalidGenomicRange(Position, Position),

    #[error("The specified genomic range [{0}, {1}] is invalid for a sequence of length {2}. Adjust the range to fit within the sequence length.")]
    InvalidGenomicRangeForSequence(Position, Position, Position),

    #[error("The sequence name '{0}' is not found within the provided ranges container. Check the sequence names for typos or missing entries.")]
    MissingSequence(String),

    #[error("An error was encountered with the underlying genomap::GenomeMap: {0}")]
    GenomeMapError(#[from] GenomeMapError),

    #[error("The GRanges object is invalid because it lacks an associated data container. Ensure data is loaded or associated with the GRanges object before attempting operations that require data.")]
    NoDataContainer,

    // Sequences related errors
    #[error("The sequence name '{0}' was not found in the dataset. Verify the sequence names and try again.")]
    MissingSequenceName(String),

    // FASTA/noodles related errors
    #[error("An error occurred while converting bytes to a UTF-8 string. This often indicates invalid or corrupted data.")]
    FromUtf8Error(#[from] FromUtf8Error),

    // Command line tool related errors
    #[error("The specified genomic range file format is unsupported. Check the documentation for a list of supported formats.")]
    UnsupportedGenomicRangesFileFormat,

    #[error("There was an error with the provided command line argument: {0}. Review the command line arguments for any mistakes.")]
    ArgumentError(#[from] clap::error::Error),

    #[error(
        "The requested operation '{0}' does not exist. Check the list of supported operations."
    )]
    NoSuchOperation(String),

    #[error("No operation was specified. See granges map --help.")]
    NoOperationSpecified,
}
