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

use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str::FromStr;

use crate::error::GRangesError;
use crate::io::tsv::TsvConfig;
use crate::io::InputStream;
use crate::ranges::{GenomicRangeEmptyRecord, GenomicRangeRecord};
use crate::traits::TsvSerialize;
use crate::Position;

use super::parse_column;
use super::tsv::TsvRecordIterator;

/// [`serde`] deserializer for a BED column with a possibly missing value. Note that the [BED
/// specification](https://samtools.github.io/hts-specs/BEDv1.pdf) only technically allows `'.'` to
/// be used for missing strands, but in practice it can be found to represent
/// missing scores, etc too.
pub fn deserialize_bed_missing<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de> + FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    let missing_chars = &["."];
    deserialize_option_generic(deserializer, missing_chars) // Use the generic deserializer with specific placeholders
}

/// The additional two BED5 columns.
///
/// # Fields
/// * `name`: the feature name.
/// * `score`: a score.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bed5Addition {
    pub name: String,
    #[serde(deserialize_with = "deserialize_bed_missing")]
    pub score: Option<f64>,
}

/// An iterator over BED5 entries, which contain the three
/// range entries (sequence name, start and end positions),
/// a feature name, and a score.
///
/// Note that the [`Bed5Addition`] is *permissive* in the sense
/// it allows for more input types than the official [BED format
/// specification](https://samtools.github.io/hts-specs/BEDv1.pdf).
// TODO strict type?
#[derive(Debug)]
pub struct Bed5Iterator {
    iter: TsvRecordIterator<GenomicRangeRecord<Bed5Addition>>,
}

impl Bed5Iterator {
    /// Creates a parsing iterator over a BED5 file.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let iter = TsvRecordIterator::new(filepath)?;

        Ok(Self { iter })
    }
}

impl Iterator for Bed5Iterator {
    type Item = Result<GenomicRangeRecord<Bed5Addition>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// An iterator over BED3 entries (which just contain ranges no data).
#[derive(Debug)]
pub struct Bed3Iterator {
    iter: TsvRecordIterator<GenomicRangeEmptyRecord>,
}

impl Bed3Iterator {
    /// Creates a parsing iterator over a BED5 file.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let iter = TsvRecordIterator::new(filepath)?;
        Ok(Self { iter })
    }
}

impl Iterator for Bed3Iterator {
    type Item = Result<GenomicRangeEmptyRecord, GRangesError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// Nucleotide strand enum type.
#[derive(Clone, Debug)]
pub enum Strand {
    Forward,
    Reverse,
}

impl TsvSerialize for Option<Strand> {
    #![allow(unused_variables)]
    fn to_tsv(&self, config: &TsvConfig) -> String {
        match self {
            Some(Strand::Forward) => "+".to_string(),
            Some(Strand::Reverse) => "-".to_string(),
            None => config.no_value_string.to_string(),
        }
    }
}

/// Deserializes some value of type `T` with some possible missing
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

//impl Selection for &Bed5Addition {
//    fn select_by_name(&self, name: &str) -> DatumType {
//        match name {
//            "name" => DatumType::String(self.name.clone()),
//            "score" => DatumType::Float64(self.score),
//            _ => panic!("No item named '{}'", name),
//        }
//    }
//}

impl TsvSerialize for &Bed5Addition {
    #![allow(unused_variables)]
    fn to_tsv(&self, config: &TsvConfig) -> String {
        format!(
            "{}\t{}",
            self.name,
            self.score
                .as_ref()
                .map_or(config.no_value_string.clone(), |x| x.to_string())
        )
    }
}

impl TsvSerialize for Bed5Addition {
    #![allow(unused_variables)]
    fn to_tsv(&self, config: &TsvConfig) -> String {
        format!(
            "{}\t{}",
            self.name,
            self.score
                .as_ref()
                .map_or(config.no_value_string.clone(), |x| x.to_string())
        )
    }
}

///// The additional three BED6 columns.
//#[derive(Clone, Debug)]
//pub struct Bed6Addition {
//    pub name: String,
//    pub score: f64,
//    pub strand: Option<Strand>,
//}

/// A lazy parser for BED-like files.
/// yields [`GenomicRangeRecord<Option<Vec<String>>>`] entries. If the file is a BED3 file,
/// the data in the [`GenomicRangeRecord`] will be set to `None`, since there are no remaining
/// string columns to parse.
pub struct BedlikeIterator {
    reader: BufReader<Box<dyn std::io::Read>>,
}

impl std::fmt::Debug for BedlikeIterator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BedlikeIterator").finish_non_exhaustive()
    }
}

impl BedlikeIterator {
    /// Create a new lazy-parsing iterator over Bed-like TSV data. This parser
    /// assumes the first three columns are the sequence name, start (0-indexed and inclusive),
    /// and end (0-indeed and exclusive) positions.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let mut input_file = InputStream::new(filepath);
        let _has_metadata = input_file.collect_metadata("#", None);
        let reader = input_file.continue_reading()?;
        Ok(Self { reader })
    }
}

impl Iterator for BedlikeIterator {
    type Item = Result<GenomicRangeRecord<Option<String>>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => None,
            Ok(_) => {
                let line = line.trim_end();
                Some((parse_bed_lazy)(line))
            }
            Err(e) => Some(Err(GRangesError::IOError(e))),
        }
    }
}

/// Lazily parses a BED* format line into the first three columns defining the range,
/// storing the rest as a `String`.
///
/// Unlike [`parse_bedlike()`], this function returns a concrete
/// [`GenomicRangeRecord<Option<String>>`] type, not a [`Vec<String>`] of split TSV columns.
pub fn parse_bed_lazy(line: &str) -> Result<GenomicRangeRecord<Option<String>>, GRangesError> {
    let columns: Vec<&str> = line.splitn(4, '\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::Bed3TooFewColumns(
            columns.len(),
            line.to_string(),
        ));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    let data = if columns.len() > 3 {
        Some(columns[3].to_string())
    } else {
        None
    };

    Ok(GenomicRangeRecord {
        seqname,
        start,
        end,
        data,
    })
}

#[allow(dead_code)]
fn parse_strand(symbol: char) -> Result<Option<Strand>, GRangesError> {
    match symbol {
        '+' => Ok(Some(Strand::Forward)),
        '-' => Ok(Some(Strand::Reverse)),
        '.' => Ok(None),
        _ => Err(GRangesError::InvalidString),
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
