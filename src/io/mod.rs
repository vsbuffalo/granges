//! Types and methods for reading and parsing input and writing output.

pub mod file;
pub mod parsers;
pub mod tsv;

pub use file::{InputStream, OutputStream};
pub use parsers::{
    Bed3Iterator, Bed5Iterator, BedlikeIterator, GenomicRangesFile, GenomicRangesParser,
    TsvRecordIterator,
};
pub use tsv::BED_TSV;
