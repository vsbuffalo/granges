use std::path::Path;

use crate::{GRangesError, Position};

/// Parses a single column from a string slice into a specified type.
///
/// This function is particularly useful for converting columns in genomic data files
/// from string representation to a specific type like `usize`, `f64`, etc.
///
/// # Type Parameters
///
/// * `T`: The type to which the column should be parsed. It must implement `FromStr`.
///
/// # Arguments
///
/// * `column`: The string slice representing the column to be parsed.
/// * `line`: The entire line from which the column was extracted. Used for error reporting.
///
/// # Errors
///
/// Returns `GRangesError::InvalidColumnType` if the column cannot be parsed into type `T`.
pub fn parse_column<T: std::str::FromStr>(column: &str, line: &str) -> Result<T, GRangesError>
where
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    column
        .parse::<T>()
        .map_err(|_| GRangesError::InvalidColumnType {
            expected_type: std::any::type_name::<T>().to_string(), // Provides the expected type name
            found_value: column.to_string(),
            line: line.to_string(),
        })
}

/// Parses a BED-like TSV format line into its constituent components.
///
/// BED-like formats contain BED3 columns (chromosome, start, and end) and
/// others. It's an non-standard but extensible (since it's just a TSV), and sensible
/// for a variety of genomic data.
///
/// # Arguments
///
/// * `line`: A string slice representing a line in BED3 format.
///
/// # Errors
///
/// Returns `GRangesError::InvalidColumnType` if the line does not have enough columns
/// or if any column cannot be correctly parsed.
pub fn parse_bedlike(line: &str) -> Result<(String, Position, Position, Vec<&str>), GRangesError> {
    let columns: Vec<&str> = line.split('\t').collect();
    if columns.len() < 3 {
        return Err(GRangesError::BedTooFewColumns(
            columns.len(),
            3,
            line.to_string(),
        ));
    }

    let seqname = parse_column(columns[0], line)?;
    let start: Position = parse_column(columns[1], line)?;
    let end: Position = parse_column(columns[2], line)?;

    // Collect remaining columns
    let additional_columns = columns[3..].to_vec();

    Ok((seqname, start, end, additional_columns))
}

/// Get the *base* extension to help infer filetype, which ignores compression-related
/// extensions (`.gz` and `.bgz`).
pub fn get_base_extension<P: AsRef<Path>>(filepath: P) -> Option<String> {
    let path = filepath.as_ref();

    // get the filename and split by '.'
    let parts: Vec<&str> = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("")
        .split('.')
        .collect();

    let ignore_extensions = ["gz", "bgz"];

    let has_ignore_extension = parts
        .last()
        .map_or(false, |ext| ignore_extensions.contains(ext));

    if parts.len() > 2 && has_ignore_extension {
        // if it's .gz, we return the second to last token,
        // e.g. path/foo.bed.gz would return bed
        Some(parts[parts.len() - 2].to_string())
    } else if parts.len() > 1 {
        // there is no .gz - return the last token.
        Some(parts[parts.len() - 1].to_string())
    } else {
        // no extension found
        None
    }
}

#[cfg(test)]
mod tests {
    use super::get_base_extension;

    #[test]
    fn test_get_base_extension() {
        assert_eq!(get_base_extension("test.bed.gz").unwrap(), "bed");
        assert_eq!(get_base_extension("test.bed").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.bed.gz").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.bed").unwrap(), "bed");
        assert_eq!(get_base_extension("some/path/test.gff").unwrap(), "gff");
        assert_eq!(get_base_extension("some/path/test.gff.gz").unwrap(), "gff");
        assert_eq!(get_base_extension("test.gff.gz").unwrap(), "gff");
        assert_eq!(get_base_extension("test.gff").unwrap(), "gff");
        assert_eq!(get_base_extension("test"), None);
        assert_eq!(get_base_extension("foo/test"), None);
    }
}
