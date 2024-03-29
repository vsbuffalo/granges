//! Types and methods for reading and parsing input and writing output.

pub mod file;
pub mod parsers;
pub mod tsv;

pub use file::{InputStream, OutputStream};
pub use parsers::{
    bed::{Bed3Iterator, Bed4Iterator, Bed5Iterator, BedlikeIterator},
    tsv::TsvRecordIterator,
    GenomicRangesFile, GenomicRangesParser,
};
pub use tsv::{TsvConfig, BED_TSV};
