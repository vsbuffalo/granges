//! Essential TSV parsing functionality, which wraps the blazingly-fast [`csv`] crate's
//! deserialization method using [`serde`].

use csv::{DeserializeRecordsIntoIter, Reader, ReaderBuilder};
use flate2::read::GzDecoder;
use serde::de::Error as DeError;
use serde::{Deserialize, Deserializer};
use std::fs::File;
use std::io::{self, Read};
use std::path::PathBuf;
use std::str::FromStr;

use crate::error::GRangesError;

/// Build a TSV reader which ignores comment lines, works on gzip-compressed
/// files, etc.
///
/// # ⚠️ Stability
///
/// This will likely change into a type with methods.
///
/// # Developers Notes
///
/// Currently headers are ignored, and not properly passed
/// to the output. If you need this feature prioritized, please submit a
/// GitHub issue.
pub fn build_tsv_reader(
    filepath: impl Into<PathBuf>,
) -> Result<Reader<Box<dyn Read>>, GRangesError> {
    let filepath = filepath.into();
    let file = File::open(&filepath)?;
    let is_gzipped = is_gzipped_file(&filepath)?;
    let stream: Box<dyn Read> = if is_gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_reader(stream);
    Ok(reader)
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

/// An extensible TSV parser, which uses a supplied parser function to
/// convert a line into a [`GenomicRangeRecord<U>`], a range with generic associated
/// data.
pub struct TsvRecordIterator<T> {
    inner: DeserializeRecordsIntoIter<Box<dyn std::io::Read>, T>,
}

impl<T> std::fmt::Debug for TsvRecordIterator<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TsvRecordIterator").finish_non_exhaustive()
    }
}

/// Check if a file is a gzipped by looking for the magic numbers
pub fn is_gzipped_file(file_path: impl Into<PathBuf>) -> io::Result<bool> {
    let mut file = File::open(file_path.into())?;
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;

    Ok(buffer == [0x1f, 0x8b])
}

impl<T> TsvRecordIterator<T>
where
    for<'de> T: Deserialize<'de>,
{
    /// Create a new TSV reader. By default, this will skip lines that being with
    /// `'#'`, since a pseudo-standard is that these indicate metadata, or column
    /// headers.
    ///
    /// # Stability
    /// Future versions may parse comment headers or make this an option.
    /// E.g. for VCF, it would need to be parsed.
    pub fn new(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let filepath = filepath.into();
        let reader = build_tsv_reader(filepath)?;
        let inner = reader.into_deserialize();

        Ok(Self { inner })
    }
}

impl<T> Iterator for TsvRecordIterator<T>
where
    for<'de> T: Deserialize<'de>,
{
    type Item = Result<T, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|res| res.map_err(|e| GRangesError::IOError(e.into())))
    }
}
