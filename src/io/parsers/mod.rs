//! Various types of parsing iterators for range formats.
//!
//! To work with genomic data in GRanges, one first needs to read it off disk and parse it. This is
//! done with *parsing iterators*. These iterators work on both plaintext and gzip-compressed
//! files (and bgzip and binary format support will hopefully be added). Each row of a file is
//! yielded as a particular parsed range type, which then can be filtered or altered in while
//! in the iterator, using Rust's powerful [`Iterator`] trait methods.
//!
//! Under the hood, this used the [`csv`] crate with [`serde`]. Initial parsers were handrolled,
//! but about 20%-30% less performant.
//!
//! ## Parsing Iterator Item Types
//!
//! These are the generic types that the iterator yields for each call of [`Iterator::next()`],
//! which represent a single row of data.
//!
//! There are two key item types that flow through parsing iterators, each representing a row of a
//! different kind of data file.
//!
//!  1. [`GenomicRangeEmptyRecord`], when a BED3 is loaded. This type indicates the incoming ranges
//!     *do not* have associated data that would go into a data container. It is *empty* of data.
//!
//!  2. [`GenomicRangeRecord<Option<String>>`], which is a range record *possibly* (hence,
//!     the [`Option`]) with remaining *unparsed* `String` data (e.g. the remaining unparsed
//!     columns of a VCF file).
//!
//!  3. [`GenomicRangeRecord<U>`], which is a range record with an *addition* set of data.
//!     For example, the standard BED5 has a two column [`Bed5Addition`], over the standard
//!     first three range columns. Thus, the [`Bed5Iterator`] yields
//!     [`GenomicRangeRecord<Bed5Addition>`] types.
//!
//! ## Why Empty Types?
//!
//! Most genomic file formats can be thought of as a genomic range with some associated data,
//! for example:
//!
//!  - VCF files will have information on the alleles and a vector of genotypes.
//!
//!  - GTF/GFF files store details about the annotated feature (e.g. name, protein ID, etc).
//!
//! However, often empty ranges (those that do not carry any data) are useful: for example, ranges
//! that just define genomic windows at regular widths, or masked repeat ranges. In this case, when
//! data is loaded into a [`GRanges<R, T>`] type, the data container type `T` is `()`. This is a
//! degenerate case that GRanges handles with a specific type called [`GRangesEmpty`], which indicates
//! the ranges in the ranges container in [`GRangesEmpty`] do not have any associated data.
//!
//! In most cases, people working with genomic file formats will know its type at compile-time.
//! However, if this isn't the case, see [⚙️  Handling Range Data when the type isn't know at runtime] below.
//!
//! # ⚙️  Handling Range Data when the type isn't know at runtime
//!
//! Because GRanges is a compile-time library, the type of incoming data is not runtime-polymorphic.
//! GRanges implements file detection methods to infer the filetype from the extension and
//! does additional type validation checks. When you want to use GRanges to handle a type
//! not known at compile-time, you will need to use `match` statements against the parsing
//! iterator types defined in [`GenomicRangesParser`]. The best examples of this are the
//! `granges` subcommands implementations.
//!
//! [`GRanges`]: crate::granges::GRanges
//! [`GRanges<R, T>`]: crate::granges::GRanges
//! [`GRangesEmpty`]: crate::granges::GRangesEmpty
//!

pub mod bed;
pub mod detect;
pub mod filters;
pub mod tsv;
pub mod utils;

pub use bed::{Bed3Iterator, Bed5Addition, Bed5Iterator, BedlikeIterator};
pub use detect::{GenomicRangesFile, GenomicRangesParser};

pub use filters::{FilteredRanges, UnwrappedRanges};
