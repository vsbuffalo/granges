// Copyright (2024) Vince Buffalo
#![crate_name = "granges"]
#![doc(html_root_url = "https://docs.rs/granges/")]

//! # GRanges: Generic Genomic Range and Data Containers
//!
//! GRanges is a Rust library for working with genomic ranges and their associated data. GRanges
//! also provides a separate command line tool built using this library, which currently
//! implements a subset of the features of other great bioinformatics utilities like
//! [bedtools](https://bedtools.readthedocs.io/en/latest/).
//!
//! The idea of GRanges is to facilitate building a new generation of robust, reproducible, and
//! easily installable bioinformatics tools in Rust. It aims to lower the barrier to building
//! simple bioinformatics data processing tools in Rust for biologists that are interested
//! in using the language in their work. Even biologists that have no interest in learning Rust
//! should benefit from Rust-based tools since they are fast. Currently GRanges is engineered
//! as a set of compile-time generic types and methods, though future versions may implement
//! analagous runtime dynamic data structures that would allow the library to be wrapped and used
//! by other languages.
//!
//!
//! ## Design Overview
//!
//! To understand how GRanges could be useful in processing genomic data, it's helpful to first
//! understand the overall design of this software library.
//!
//! Nearly all range data is loaded into GRanges using **parsing iterators** ([parsers]).
//! Parsing iterators read input data from some bioinformatics file format (e.g. `.bed`,
//! `.gtf.gz`, etc) and it then goes down one of two paths:
//!
//!  1. Streaming mode: the entire data set is never read into memory all at once, but rather
//!     is processed entry by entry.
//!
//!  2. In-memory mode: the entire data set *is* read into memory all at once, usually as
//!     a [`GRanges<R, T>`] object. The [`GRanges<R, T>`] object is composed of two parts:
//!     a *ranges container* ([ranges]) and a *data container* ([data]).
//!     [`GRanges`] objects then provide high-level in-memory methods for working with this
//!     genomic data.
//!
//! Both modes then write something to output: a TSV of some statistic calculated on the
//! genomic data, the output of 10 million block bootstraps, or the coefficients of
//! linear regressions run every kilobase on sequence data.
//!
//! The GRanges workflow idea is illustrated in this (super cool) ASCII diagram:
//!
//!
//! ```text
//!   +-------+           +-------------------+           +----------------------+          
//!   | Input |  ----->   |     Iterator      |  ----->   |    GRanges object    | 
//!   +------ +           | (streaming mode)  |           |   (in-memory mode)   |          
//!                       +-------------------+           +----------------------+          
//!                                 |                                |            
//!                                 |                                |            
//!                                 v                                v            
//!                        +-----------------+               +-----------------+       
//!                        | Data Processing |               | Data Processing |       
//!                        +-----------------+               +-----------------+       
//!                                 |                                |            
//!                                 v                                v            
//!                             +--------+                       +--------+       
//!                             | Output |                       | Output |       
//!                             +--------+                       +--------+       
//!
//! ```
//!
//! Nearly all common data processing and analytics tasks in science can be thought of as taking
//! some raw input data and refining it down some pipeline, until eventually outputting it as some sort of statistical
//! result. Bioinformatics has been built on composing workflows out of command line tools and
//! pipes. Analogously, tidyverse unifies data tidying, summarizing, and modeling by passing data
//! in a datafame down a similar pipeline-based interface. Both operations are simple to use
//! because they compose the same essential *operations* into specific data processing workflows.
//!
//! GRanges is designed to emulate these powerful data analytics tools. It implements core
//! operations on genomic ranges and their associated data, and provides generic, extensible
//! data structures for storing genomic data.
//!
//! GRanges attempts to implement the same core set of data joining and processing tools as the
//! [plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html) library and
//! R's [tidyverse](https://www.tidyverse.org). It also frames these genomic operations in
//! standard data analytic and database terms. For example, the operation `bedtools intersect -wa -u -a <left>
//! -b <right>` finds all ranges in the *left* set of genomic ranges that overlaps one or more
//! *right* genomic ranges, and reports only the unique left range. This operation is common in
//! many data analytics workflows and database queries, and is known as a *semi join* between a
//! a left table and right table. (see for example this section in Hadley Wickham's book [R for Data Science](https://r4ds.hadley.nz/joins.html#filtering-joins)).
//! In the `granges` command line tool, this is implemented in the `granges filter` subcommand. In
//! the library, this subcommand is just the method [`GRanges::filter_overlaps()`].
//!
//!
//! ## Documentation Guide
//!
//!  - ⚙️ : technical details that might only be useful to developers handling tricky cases.
//!  - ⚠️: indicates a method or type that is unstable (it may change in future implementations),
//!     or is "sharp" (developer error could cause a panic).
//!  - 🚀: optimization tip.
//!
//! [`GRanges`]: crate::granges::GRanges
//! [`GRanges::filter_overlaps()`]: crate::granges::GRanges
//! [`GRanges<R, T>`]: crate::granges::GRanges
//! [parsers]: crate::io::parsers

pub use indexmap;

pub mod data;
pub mod error;
pub mod granges;
pub mod io;
pub mod iterators;
pub mod join;
pub mod ranges;
pub mod test_utilities;
pub mod traits;

// bringing these CLI modules into lib.rs rather than main/ allows for
// use in integration tests and other Rust-side command line work
pub mod commands;
pub mod reporting;

/// The main position type in GRanges.
pub type Position = u32;

/// The main *signed* position type in GRanges, to represent offsets (e.g.
/// for adjust range coordinates, etc).
pub type PositionOffset = i32;

/// The main exports of the GRanges library.
pub mod prelude {
    pub use crate::error::GRangesError;
    pub use crate::granges::{GRanges, GRangesEmpty};
    pub use crate::io::file::read_seqlens;
    pub use crate::io::{
        Bed3Iterator, BedlikeIterator, GenomicRangesFile, GenomicRangesParser, TsvRecordIterator,
    };

    pub use crate::ranges::vec::{VecRangesEmpty, VecRangesIndexed};
    pub use crate::traits::{
        AsGRangesRef, GeneralRangeRecordIterator, GenericRange, GenericRangeOperations,
        GenomicRangeRecordUnwrappable, GenomicRangesTsvSerialize, IndexedDataContainer,
        IntoIterableRangesContainer, IterableRangeContainer, TsvSerialize,
    };

    pub use crate::seqlens;
}

pub const INTERNAL_ERROR_MESSAGE: &str = r#"
An internal error has occurred. Please file a GitHub issue at 
https://github.com/vsbuffalo/granges/issues
"#;

/// Create an [`IndexMap`] of sequence names and their lengths.
///
/// [`IndexMap`]: indexmap::IndexMap
#[macro_export]
macro_rules! seqlens {
    ($($key:expr => $value:expr),* $(,)?) => {
        $crate::indexmap::indexmap!($($key.to_string() => $value),*)
    };
}

/// Create a new [`GRanges<R, T>`] with sequence length information (used primarily for small examples)
///
/// [`GRanges<R, T>`]: crate::granges::GRanges
///
#[macro_export]
macro_rules! create_granges_with_seqlens {
    ($range_type:ty, $data_type:ty, { $($chr:expr => [$(($start:expr, $end:expr, $data:expr)),*]),* }, seqlens: { $($chr_len:expr => $len:expr),* }) => {
        {
            let mut seqlens = ::indexmap::IndexMap::new();
            $(seqlens.insert($chr_len.to_string(), $len);)*

            let mut gr: GRanges<$range_type, $data_type> = GRanges::new_vec(&seqlens);

            $(
                $(
                    gr.push_range(&$chr.to_string(), $start, $end, $data).expect("Failed to push range");
                )*
            )*

            gr
        }
    };
}

/// a version of [assert_eq!] that is more user-friendly.
#[macro_export]
macro_rules! ensure_eq {
    ($left:expr, $right:expr $(,)?) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                if !(*left_val == *right_val) {
                    panic!("{}\nExpected `{}` but found `{}`.", $crate::INTERNAL_ERROR_MESSAGE, stringify!($left), stringify!($right));

                }
            }
        }
    });
    ($left:expr, $right:expr, $($arg:tt)+) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                if !(*left_val == *right_val) {
                    panic!("{}\n{}\nExpected `{}` but found `{}`.", $crate::INTERNAL_ERROR_MESSAGE, format_args!($($arg)+), stringify!($left), stringify!($right));
                }
            }
        }
    });
}

#[cfg(test)]
mod tests {
    #[test]
    #[should_panic]
    fn test_ensure_eq() {
        ensure_eq!(1, 2);
    }

    #[test]
    #[should_panic]
    fn test_ensure_eq_msg() {
        ensure_eq!(1, 2, "Some description of the problem.");
    }

    #[test]
    fn test_ensure_eq_pass() {
        ensure_eq!(1, 1);
    }
}
