//! Input/Output
//!

pub mod io;
pub mod noodles;
pub mod parsers;

pub use io::{InputFile, OutputFile};
pub use parsers::{Bed3RecordIterator, TsvRecordIterator, BedlikeIterator};
