//! Input/Output
//!

pub mod file;
pub mod noodles;
pub mod parsers;

pub use file::{InputFile, OutputFile};
pub use parsers::{BedlikeIterator, TsvRecordIterator};
