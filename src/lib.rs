// Copyright (2024) Vince Buffalo
#![crate_name = "granges"]
#![doc(html_root_url = "https://docs.rs/granges/")]

//! # GRanges: Generic Genomic Range and Data Containers
//!
//! GRanges is a Rust library for working with genomic ranges and their associated data. GRanges
//! also implements a [bedtools](https://bedtools.readthedocs.io/en/latest/)-like command line tool
//! as a test and benchmark of the library. In benchmarks, this command line tool is reliably
//! 30-40% faster than bedtools. Internally, GRanges uses the very fast
//! [coitrees](https://github.com/dcjones/coitrees) library written by Daniel C. Jones for overlap
//! operations.
//!
//! The GRanges library aims to simplify the creation of powerful, performant genomics tools in
//! Rust. GRanegs is inspired by the design and ease of use of Bioconductor's
//! [GenomicRanges](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118)
//! and [plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)
//! packages, but tools using GRanges are optimized at compile time and thus are typically much
//! faster.
//!
//! ## GRanges Design
//!
//! Nearly all common data processing and analytics tasks in science can be thought of as taking
//! some raw input data and refining it down some pipeline, until eventually outputting some sort
//! of statistical result. In bioinformatics, genomics data processing workflows are often composed
//! from command line tools and Unix pipes. Analogously, tidyverse unifies data tidying,
//! summarizing, and modeling by passing and manipulating a datafame down a similar pipeline. Unix
//! pipes and tidyverse are powerful because they allow the scientist to build highly specialized
//! analytic tools just by remixing the same core set of data operations. 
//!
//! GRanges implements many of the same core data joining and manipulation operations as
//! [plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html) library and
//! R's [tidyverse](https://www.tidyverse.org). 
//!
//! ## The [`GRanges`] container
//!
//! Central to the GRanges library is the generic [`GRanges<R, T>`] data container. This is generic
//! over the *range container* type `R`, because range data structures have inherent performance
//! tradeoffs for different tasks. For example, when parsing input data into memory, the total size
//! is often not known in advanced. Consequently, the most efficient range container data structure
//! is a dynamically-growing [`Vec<U>`]. However, when overlaps need to be computed between ranges
//! across two containers, the *interval tree* data structures are needed. [`GRanges`] supports
//! conversion between these range container types. 
//!
//! [`GRanges<R, T>`] is also generic over its *data container*, `T`. This allows it to hold *any*
//! type of data, in any data structure: a [`Vec<f64>`] of numeric stores, an
//! [`ndarray::Array2`](https://docs.rs/ndarray/latest/ndarray/index.html) matrix, a
//! [polars](https://pola.rs) dataframe, custom `struct`s, etc. Each range in the ranges container
//! indexes an element in the data container. There is also a special [`GRangesEmpty`] type for
//! data-less [`GRanges`] objects.
//!
//! The [`GRanges`] type implements a small but powerful set of methods that can be used to process
//! and analyze genomic data. GRanges relies heavily on *method-chaining* operations that pass
//! genomic range data down a pipeline, similar to tidyverse and Unix pipes.
//!
//! ## The Overlaps--Map-Combine Pattern
//!
//! A very common data analysis strategy is **Split-Apply-Combine** pattern ([Wickham,
//! 2011](https://vita.had.co.nz/papers/plyr.html)). When working with genomic overlaps
//! and their data, GRanges library relies on its own data analysis pattern known 
//! as **Overlaps--Map-Combine** pattern, which is better suited for common genomic operations and analyses.
//!
//! ### Overlap Joins and the Left Grouped Join
//!
//! In this pattern, the *left* genomic ranges are "joined" by their *overlaps* with another
//! *right* genomic ranges. In other words, the *right ranges* that overlap a single *left range* can 
//! be thought of "joined" to this left range.
//!
//! With genomic data, downstream processing is greatly simplified if the results of this join 
//! are *grouped* by what left range they overlap (see illustration below). This is a special 
//! kind of join known as a **left grouped join**. 
//!
//!
//! ```text
//!     Left range:   ███████████████████████████████████
//!                                                      
//!   Right ranges:     ██    ██████        ██████    ██████
//! ```
//!
//! In GRanges, a left grouped join ([`GRanges::left_overlaps()`]) returns a [`GRanges`] object
//! containing the exact same ranges as the input *left ranges*. In other words, *left grouped
//! joins* are *endomorphic in the range container*. Computationally, this is also extremely
//! efficient, because the left ranges can be passed through to the results directly (a
//! zero-overheard operation). 
//!
//! 
//!
//! The [`GRanges`] object returned by [`GRanges::left_overlaps()`] contains a [`JoinData`] (or related)
//! type as its data container. The [`JoinData`] contains information about each left range's overlaps (e.g. number of 
//! overlapping bases, fraction of overlaps, etc) in a [`Vec<LeftGroupedJoin>`]. This information about 
//! overlapping ranges may then be used by downstream calculations.
//! Additionally [`JoinData`] stores the left ranges' data container, and has a *reference* to 
//! the right ranges' data container. In both cases, the data is *unchanged* by the join. 
//! Downstream processing can also easily access the left ranges's data and all overlapping right ranges'
//! data for calculations.
//! 
//! ### Map-Combine over Joins
//!
//! After a left grouped join, each left range can have zero or more overlapping right ranges.
//! A *map-combine* operation (optionally) maps a function over all right ranges' data and
//! overlaps information, and then combines that into a single data entry per left range. These
//! combined data entries make up a new [`GRanges<R, Vec<V>>`] data container, returned by [`GRanges::map_over_joins()`].
//!
//! Note that like [`GRanges::left_overlaps()`], the [`GRanges::map_over_joins()`] is endomorphic 
//! over its range container. This means it can be passed through without modification, which 
//! is computationally efficient. This results from a Map-Combine operation then can be overlap joined with 
//! other genomic ranges, filtered, have its data arbitrarily manipulated by [`GRanges::map_data()`], etc.
//!
//! ## Example
//!
//! To illustrate these ideas, lets look at an example of how we might use GRanges to do a
//! a commonly-encountered genomic calculation. Suppose you had a set of exon genomic ranges (the *left
//! ranges*) and a multicolumn TSV (*right ranges*) of various genomic features (e.g. assay
//! results, recombination hotspots, variants, etc). Imagine you wanted to form some statistic per
//! each exon, based on data in the overlapping right ranges. This operation is very
//! similar to `bedtools map`, except that suppose the statistic you want to calculate is not one of the 
//! few offered by bedtools. GRanges makes it simple to compose your own, very fast genomic tools to answer these
//! questions.
//!
//! To see how GRanges would make it simple to compose a specialized fast tool to solve this
//! problem, let's fist see how few lines of code it would take to implement `bedtools map`,
//! since our problem is akin to using `bedtools map` with our own function to calculate our statistic.
//! Let's start by getting the mean score by seeing how to get the mean BED5 score across all 
//! overlapping right ranges for each left range (i.e. `bedtools map -a <left> -b <right> -c 5 mean`). 
//! Here is the Rust code to do this using GRanges:
//!
//! ```
//! # use granges::prelude::*;
//! # fn try_main() -> Result<(), granges::error::GRangesError> {
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
//! let left_iter = Bed3Iterator::new(left_path).expect("error inferring filetype");
//! let right_iter = Bed5Iterator::new(right_path).expect("error inferring filetype");
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
//! # Ok(())
//! # }
//! # fn main() { try_main().unwrap(); }
//! ```
//!
//! In relatively few lines of code, we can implement core functionality of bedtools.
//! As mentioned earlier, GRanges code typically runs 30-40% faster than bedtools.
//!
//! ## Loading Data into GRanges
//!
//! Nearly all genomic range data enters GRanges through **parsing iterators** (see the [parsers]
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
//! In the example above, we loaded the items in the parsing iterators directly into 
//! [`GRanges`] objects, since we had to do overlap operations (in the future, GRanges will
//! support streaming overlap joins for position-sorted input data).
//!
//! Both processing modes eventually output something, e.g. a TSV of some statistic calculated 
//! on the genomic data, the output of 10 million block bootstraps, or the coefficients of
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
//! ## Manipulating GRanges objects 
//!
//! 1. *Creation*: [`GRanges::new_vec()`], [`GRanges::from_iter()`], [`GRangesEmpty::from_windows()`].
//!
//! 2. *Range-modifying* functions: [`GRanges::into_coitrees()`], [`GRanges::adjust_ranges()`], [`GRanges::sort()`],
//!        [`GRanges::flanking_ranges()`], [`GRanges::filter_overlaps()`].
//!
//! 3. *Data-modifying* functions: [`GRanges::left_overlaps()`], [`GRanges::map_over_joins()`], [`GRanges::map_data()`].
//!
//! 4. *Range and Data modifying* functions: [`GRanges::push_range()`].
//!
//! Note that not all range and data container types support these operations, e.g.
//! [`GRanges::push_range()`] is not available for ranges stored in an interval tree.
//! However, by design, this is known *at compile time*, due Rust's typing system. Thus 
//! there is lower risk of runtime panics, since more potential issues are caught at 
//! compile time.
//!
//! ## The (Active) Future
//!
//! When I was first learning Rust as a developer, one thing that struck me about the language
//! is that it feels *minimal but never constraining*. This is quite a contrast compared to 
//! languages like C++. As Rust developers know, the minimalness was by design; the
//! cultural norm of slow growth produced a better final product in the end. 
//!
//! The GRanges library attempts to follow this design too. Currently, the alpha release is
//! about implementating a subset of essential core functionality to benchmark against 
//! alternative software. Since the design and ergonomics of the API are under active development,
//! please, *please*, file a GitHub issue if:
//!
//!  1. You want a particular feature.
//!
//!  2. You don't know how you'd implement a feature yourself with GRanges.
//!
//!  3. The ["ergonomics"](https://blog.rust-lang.org/2017/03/02/lang-ergonomics.html) don't
//!     feel right.
//!
//! ## Documentation Guide
//!
//!  - ⚙️ : technical details that might only be useful to developers handling tricky cases.
//!  - ⚠️: indicates a method or type that is unstable (it may change in future implementations),
//!     or is "sharp" (developer error could cause a panic).
//!  - 🚀: optimization tip.
//!
//! [`JoinData`]: crate::join::JoinData
//! [`GRanges`]: crate::granges::GRanges
//! [`GRangesEmpty`]: crate::granges::GRangesEmpty
//! [`GRanges::filter_overlaps()`]: granges::GRanges::filter_overlaps
//! [`GRanges<R, T>`]: crate::granges::GRanges
//! [parsers]: crate::io::parsers
//! [`ndarray::Array2`]: ndarray::Array2
//! [`GRanges::left_overlaps()`]: crate::traits::LeftOverlaps::left_overlaps
//! [`GRanges<R, Vec<V>>`]: crate::granges::GRanges
//! [`GRanges::map_over_joins()`]: crate::granges::GRanges::map_over_joins
//! [`GRanges::map_data()`]: crate::granges::GRanges::map_data
//! [`GRanges::new_vec()`]: crate::granges::GRanges::new_vec
//! [`GRanges::into_coitrees()`]: crate::granges::GRanges::into_coitrees
//! [`GRanges::adjust_ranges()`]: crate::granges::GRanges::adjust_ranges
//! [`GRanges::push_range()`]: crate::granges::GRanges::push_range
//! [`GRanges::flanking_ranges()`]: crate::granges::GRanges::flanking_ranges
//! [`GRanges::sort()`]: crate::granges::GRanges::sort
//! [`GRanges::from_iter()`]: crate::granges::GRanges::from_iter
//! [`GRangesEmpty::from_windows()`]: crate::granges::GRangesEmpty::from_windows
// TODO broken links, gr!

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
/// that reaches 5.4Gb (<https://www.nature.com/articles/s41586-021-03198-8l>).
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
