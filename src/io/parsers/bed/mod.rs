//! BED Types and Functionality
//!
//! The BED (Browser Extensible Format) is a TSV format in bioinformatics.
//! It has a fairly strict [specification](https://samtools.github.io/hts-specs/BEDv1.pdf),
//! but in practice it is quite permissive, and in bioinformatics one encounters lots
//! of "BED-like" files.
//!
//! # Design
//!
//! Since BED files can be thought of BED3 + an optional *addition*, the type
//! returned by these parsing iterators is [`GenomicRangeRecord<U>`] where
//! `U` is some sort of addition type like [`Bed5Addition`].
//!
//! # ⚠️ Stability
//!
//! This module defines core BED types, but is under active development.
//!

pub mod bed3;
pub mod bed4;
pub mod bed5;
pub mod bedlike;

pub use bed3::Bed3Iterator;
pub use bed4::{Bed4Addition, Bed4Iterator};
pub use bed5::{Bed5Addition, Bed5Iterator};
pub use bedlike::{valid_bedlike, BedlikeIterator};

use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer};
use std::str::FromStr;

/// [`serde`] deserializer for a BED column with a possibly missing value. Note that the [BED
/// specification](https://samtools.github.io/hts-specs/BEDv1.pdf) only technically allows `'.'` to
/// be used for missing strands, but in practice it can be found to represent
/// missing scores, etc too.
pub fn bed_missing<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de> + FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    let missing_chars = &["."];
    deserialize_option_generic(deserializer, missing_chars) // Use the generic deserializer with specific placeholders
}

/// Nucleotide strand enum type.
#[derive(Clone, Debug)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Deserializes some value of type `t` with some possible missing
/// character `missing_chars` into [`Option<T>`].
pub fn deserialize_option_generic<'de, D, T>(
    deserializer: D,
    missing_chars: &'de [&'de str],
) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de> + FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    if missing_chars.contains(&s.as_str()) {
        Ok(None)
    } else {
        s.parse::<T>()
            .map(Some)
            .map_err(|e| DeError::custom(format!("parsing error: {}", e)))
    }
}
// TODO
///// Parses a BED6 format line into the three columns defining the range, and additional
///// columns
/////
//pub fn parse_bed6(line: &str) -> Result<GenomicRangeRecord<Bed5Addition>, GRangesError> {
//    let columns: Vec<&str> = line.splitn(4, '\t').collect();
//    if columns.len() < 3 {
//        return Err(GRangesError::BedlikeTooFewColumns(line.to_string()));
//    }
//
//    let seqname = parse_column(columns[0], line)?;
//    let start: Position = parse_column(columns[1], line)?;
//    let end: Position = parse_column(columns[2], line)?;
//
//    let name = parse_column(columns[3], line)?;
//    // let strand: Option<Strand> = parse_strand(parse_column(columns[3], line)?)?;
//
//    let data = Bed5Addition { name, score };
//
//    Ok(GenomicRangeRecord {
//        seqname,
//        start,
//        end,
//        data,
//    })
//}
