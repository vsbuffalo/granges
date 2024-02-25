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
//! ## Example
//!
//! The best introduction to the GRanges library is an example illustrating
//! how a common genomic data processing task can be expressed in few lines
//! of (relatively simple) code. Below is an annotated implementation of
//! `granges map`, which provides the same functionality of `bedtools map`:
//! arbitrary operations (e.g. mean, median, max, etc) are applied to the numeric
//! scores of all right ranges that overlap a each left range:
//!
//! ```
//! use granges::prelude::*;
//!
//! // Data for example:
//! let genome = seqlens!("chr1" => 100, "chr2" => 100);
//! let left_path = "tests_data/bedtools/map_a.txt";
//! let right_path = "tests_data/bedtools/map_b.txt";
//!
//! // Read in the "genome file" of chromosomes and their lengths.
//! let seqnames: Vec<String> = genome.keys().cloned().collect();
//!
//! // Create parsing iterators to the left and right BED files.
//! let left_iter = Bed3Iterator::new(left_path)
//!                   .expect("error inferring filetype");
//! let right_iter = Bed5Iterator::new(right_path)
//!                   .expect("error inferring filetype");
//!
//! // Filter out any ranges from chromosomes not in our genome file.
//! let left_gr = GRangesEmpty::from_iter(left_iter.retain_seqnames(&seqnames), &genome)
//!                  .expect("error parsing file");
//! let right_gr = GRanges::from_iter(right_iter.retain_seqnames(&seqnames), &genome)
//!                  .expect("error parsing file");
//!
//!
//! // Create the "right" GRanges object, convert the ranges to an
//! // interval trees, and tidy it by selecting out a f64 score.
//! let right_gr = {
//!     right_gr
//!         // Convert to interval trees.
//!         .into_coitrees().expect("error computing interval trees")
//!         // Extract out just the score from the additional BED5 columns.
//!         .map_data(|bed5_cols| {
//!             bed5_cols.score
//!         }).expect("error selecting score")
//! };
//!
//! // Find the overlaps by doing a *left grouped join*.
//! let left_join_gr = left_gr.left_overlaps(&right_gr)
//!                        .expect("error in computing overlaps");
//!
//! // Process all the overlaps.
//! let result_gr = left_join_gr.map_over_joins(|join_data| {
//!     // Get the "right data" -- the BED5 scores.
//!     let overlap_scores = join_data.right_data;
//!     let score_sum: f64 = overlap_scores.iter().sum();
//!     score_sum / (overlap_scores.len() as f64)
//! }).expect("error computing mean score");
//!
//! // Write to a TSV file, using the BED TSV format standards
//! // for missing values, etc.
//! let path = Some("map_results.bed.gz");
//! result_gr.to_tsv(path, &BED_TSV).expect("error writing output");
//! ```
//!
//! ## Design Overview
//!
//! To understand how GRanges could be useful in processing genomic data, it's helpful to first
//! understand the overall design of this software library.
//!
//! Nearly all range data is loaded into GRanges using **parsing iterators** (see the [parsers]
//! module). Parsing iterators read input data from some bioinformatics file format (e.g. `.bed`,
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
//! Both modes then eventually write something to output: a TSV of some statistic calculated on the
//! genomic data, the output of 10 million block bootstraps, or the coefficients of
//! linear regressions run every kilobase on sequence data.
//!
//! The GRanges workflow idea is illustrated in this (super cool) ASCII diagram:
//!
//!
//! ```text
//!   +-------+           +-------------------+           +----------------------+          
//!   | Input |  ----->   |  Parsing Iterator |  ----->   |    GRanges object    |
//!   +------ +           |  (streaming mode) |           |   (in-memory mode)   |          
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
pub mod sequences;
pub mod test_utilities;
pub mod traits;

// bringing these CLI modules into lib.rs rather than main/ allows for
// use in integration tests and other Rust-side command line work
pub mod commands;
pub mod reporting;

/// The main position type in GRanges.
///
/// This type is currently an unwrapped [`u32`]. This should handle
/// chromosome lengths for nearly all species. In fact, the only exception
/// known so far is lungfush (*Neoceratodus forsteri*), which has a chromosomes
/// that reaches 5.4Gb (https://www.nature.com/articles/s41586-021-03198-8l).
/// The [`u32::MAX`] is 4,294,967,295, i.e. 4.29 Gigabases, which means [`u32`] is
/// just barely suitable for even the largest known chromosome. There is a
/// performance and memory-efficiency tradeoff when using [`u64`] over [`u32`],
/// so [`u32`] is used by default since it handles nearly all cases.
///
/// # Feature support for large chromosomes
///
/// If you are working with data from a species with unusually large chromosomes,
/// you can compile GRanges using the `--features=big-position` option, which will set
/// the [`Position`] and [`PositionOffset`] to [`u64`] and [`i64`], respectively.
///
/// [`u32::MAX`]: std::u32::MAX
#[cfg(not(feature = "big-position"))]
pub type Position = u32;
#[cfg(feature = "big-position")]
pub type Position = u64;

/// The main *signed* position type in GRanges, to represent offsets (e.g.
/// for adjust range coordinates, etc).
#[cfg(not(feature = "big-position"))]
pub type PositionOffset = i32;
#[cfg(feature = "big-position")]
pub type PositionOffset = i64;

/// The main exports of the GRanges library.
pub mod prelude {
    pub use crate::error::GRangesError;
    pub use crate::granges::{GRanges, GRangesEmpty};
    pub use crate::io::file::read_seqlens;
    pub use crate::io::tsv::BED_TSV;
    pub use crate::io::{
        Bed3Iterator, Bed5Iterator, BedlikeIterator, GenomicRangesFile, GenomicRangesParser,
        TsvRecordIterator,
    };

    pub use crate::data::DatumType;
    pub use crate::ranges::vec::{VecRangesEmpty, VecRangesIndexed};
    pub use crate::traits::{
        AsGRangesRef, GeneralRangeRecordIterator, GenericRange, GenericRangeOperations,
        GenomicRangeRecordUnwrappable, GenomicRangesTsvSerialize, IndexedDataContainer,
        IntoDatumType, IntoIterableRangesContainer, IterableRangeContainer, JoinDataOperations,
        LeftOverlaps, TsvSerialize,
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
