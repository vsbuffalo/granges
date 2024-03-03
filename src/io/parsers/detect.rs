//! Filetype detection functionality.
//!

use serde::Deserialize;
use std::path::PathBuf;

use super::{
    bed::{valid_bedlike, Bed4Addition, Bed4Iterator},
    tsv::build_tsv_reader,
    utils::get_base_extension,
    Bed3Iterator, Bed5Addition, Bed5Iterator, BedlikeIterator,
};
use crate::{
    ranges::{GenomicRangeRecord, GenomicRangeRecordEmpty},
    GRangesError,
};

/// Enum that connects a genomic ranges file type to its specific parser.
#[derive(Debug)]
pub enum GenomicRangesParser {
    Bed3(Bed3Iterator),
    Bed4(Bed4Iterator),
    Bed5(Bed5Iterator),
    Bedlike(BedlikeIterator),
    Unsupported,
}

/// Enum that indicates the filetype of some genomic ranges file.
#[derive(Debug, PartialEq)]
pub enum GenomicRangesFile {
    Bed3(PathBuf),
    Bed4(PathBuf),
    Bed5(PathBuf),
    Bedlike(PathBuf),
    Unsupported,
}

/// Detect the BED variant by trying to deserialize the first line.
pub fn detect_bed_variant(
    filepath: impl Into<PathBuf>,
) -> Result<Option<GenomicRangesFile>, GRangesError> {
    let filepath = filepath.into();

    if try_deserialize::<GenomicRangeRecord<Bed5Addition>>(&filepath)? {
        Ok(Some(GenomicRangesFile::Bed5(filepath)))
    } else if try_deserialize::<GenomicRangeRecord<Bed4Addition>>(&filepath)? {
        Ok(Some(GenomicRangesFile::Bed4(filepath)))
    } else if try_deserialize::<GenomicRangeRecordEmpty>(&filepath)? {
        Ok(Some(GenomicRangesFile::Bed3(filepath)))
    } else {
        Ok(None)
    }
}

/// Try to deserialize into a generic type `T`.
fn try_deserialize<T: for<'de> Deserialize<'de> + std::fmt::Debug>(
    filepath: impl Into<PathBuf>,
) -> Result<bool, GRangesError> {
    let filepath = filepath.into();
    let reader = build_tsv_reader(&filepath)?;
    let mut iter = reader.into_deserialize::<T>();
    let next_item = iter.next();
    if let Some(result) = next_item {
        Ok(result.is_ok())
    } else {
        Err(GRangesError::EmptyFile(
            filepath.to_string_lossy().to_string(),
        ))
    }
}

impl GenomicRangesFile {
    /// Detect the type of range genomic range file type we are working with, and output
    /// the appropriate [`GenomicRangesFile`] enum variant.
    ///
    /// Detection works like this:
    ///  1. Skip comment lines, starting with `#`.
    ///  2. Retrieve extension, removing any additional compression-related
    ///     extensions (`.gz` and `.bgz`) if present.
    ///  3. Read the first line, and try to parse the second and third columns
    ///     into [`Position`] types, to check if this is a BED-like file (i.e.
    ///     the first three columns are sequence name, start, and end positions).
    ///  4. Match against
    ///
    /// Currently this supports:
    ///  1. BED3 files, without range data.
    ///  2. BED-like files, with BED3 first columns and optional additional columns
    ///     after that. If there is no additional columnar data after the first
    ///     three BED3 columns, the type is [`GenomicRangesFile::Bed3`]. If there
    ///     is additional columnar data, it is [`GenomicRangesFile::Bedlike`].
    ///  3. Files with `.tsv` extension but are BED-like (first three
    ///     columns are BED3) will have the type [`GenomicRangesFile::Bed3`]. This
    ///     is because downstream [`GRanges`] operations need to know if any
    ///     additional data is present, which would need to be put in a data container.
    ///  4. BED5 files, which are BED3 + a *feature name* and a *strand* column.
    ///  5. If the file type does not satisfy any of the rules above, it is
    ///     [`GenomicRangesFile::Unsupported`].
    ///
    /// See the `match` statement in the source code for the exact rules. Additional
    /// genomic range file formats like GTF/GFF can easily be added later.
    ///
    /// [`GRanges`]: crate::granges::GRanges
    pub fn detect(filepath: impl Into<PathBuf>) -> Result<Self, GRangesError> {
        let filepath: PathBuf = filepath.into();

        let is_valid_bedlike = valid_bedlike(&filepath)?;

        // get the extension, as a hint
        let extension =
            get_base_extension(&filepath).ok_or(GRangesError::CouldNotDetectRangesFiletype)?;

        // If it's got a .tsv extension, take this as a hint it *isn't a BED*,
        // thus, this goes to the BedlikeIterator parser.
        if extension.ends_with("tsv") && is_valid_bedlike {
            return Ok(GenomicRangesFile::Bedlike(filepath));
        }

        // Let's try the strict serde-based deserialization approach first.
        if let Some(bed_filetype) = detect_bed_variant(&filepath)? {
            return Ok(bed_filetype);
        }

        if is_valid_bedlike {
            return Ok(GenomicRangesFile::Bedlike(filepath));
        }

        Ok(GenomicRangesFile::Unsupported)
    }

    /// Detect the genomic range filetype and link it to its parsing iterator, or raise an error
    /// if the filetype is not supported.
    ///
    /// This returns a [`GenomicRangesParser`] enum, since the parsing iterator filetype
    /// cannot be known at compile time.
    pub fn parsing_iterator(
        filepath: impl Clone + Into<PathBuf>,
    ) -> Result<GenomicRangesParser, GRangesError> {
        let path = filepath.into();
        match Self::detect(path)? {
            GenomicRangesFile::Bed3(path) => {
                Ok(GenomicRangesParser::Bed3(Bed3Iterator::new(path)?))
            }
            GenomicRangesFile::Bed4(path) => {
                Ok(GenomicRangesParser::Bed4(Bed4Iterator::new(path)?))
            }
            GenomicRangesFile::Bed5(path) => {
                Ok(GenomicRangesParser::Bed5(Bed5Iterator::new(path)?))
            }
            GenomicRangesFile::Bedlike(path) => {
                Ok(GenomicRangesParser::Bedlike(BedlikeIterator::new(path)?))
            }
            GenomicRangesFile::Unsupported => Err(GRangesError::UnsupportedGenomicRangesFileFormat),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GenomicRangesFile;

    #[test]
    fn test_rangefiletype_detect() {
        let range_filetype = GenomicRangesFile::detect("tests_data/example.bed");
        assert!(matches!(
            range_filetype.unwrap(),
            GenomicRangesFile::Bed3(_)
        ));

        let range_filetype = GenomicRangesFile::detect("tests_data/example_bedlike.tsv");
        assert!(matches!(
            range_filetype.unwrap(),
            GenomicRangesFile::Bedlike(_)
        ));

        let range_filetype = GenomicRangesFile::detect("tests_data/test_case_03.bed");
        assert!(matches!(
            range_filetype.unwrap(),
            GenomicRangesFile::Bed5(_)
        ));
    }
}
