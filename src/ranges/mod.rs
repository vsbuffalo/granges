use std::ops::Range;

use crate::{Position, error::GRangesError};

pub mod coitrees;
pub mod vec;


#[derive(Clone, Default)]
pub struct RangeEmpty {
    range: Range<Position>,
}

impl RangeEmpty {
    /// Create a new 0-indexed right-exclusive range.
    pub fn new(start: Position, end: Position) -> Self {
        let start = start.try_into().unwrap();
        let end = end.try_into().unwrap();
        Self {
            range: start..end,
        }
    }

    pub fn start(&self) -> Position {
        self.range.start
    }

    pub fn end(&self) -> Position {
        self.range.end
    }
}

#[derive(Clone, Debug, Default)]
pub struct RangeIndexed {
    range: Range<Position>,
    index: usize,
}

impl RangeIndexed {
    /// Create a new 0-indexed right-exclusive range.
    pub fn new(start: Position, end: Position, index: usize) -> Self {
        let start = start.try_into().unwrap();
        let end = end.try_into().unwrap();
        Self {
            range: start..end,
            index
        }
    }

    pub fn start(&self) -> Position {
        self.range.start
    }

    pub fn end(&self) -> Position {
        self.range.end
    }

    // Note: this returning a reference is required to 
    // implement coitrees's GenericInterval trait.
    pub fn index(&self) -> &usize {
        &self.index
    }
}

/// Validates whether a given range is valid for accessing a sequence of a given `length`.
///
/// # Arguments
///
/// * `range` - The range to validate.
/// * `length` - The length of the sequence.
///
/// # Returns
///
/// * `bool` - `true` if the range is valid for the sequence; otherwise, `false`.
pub fn validate_range(range: &std::ops::Range<Position>, length: Position) -> 
Result<(), GRangesError> {
    let start = range.start;
    let end = range.start;
    dbg!(&start);
    dbg!(&end);
    if start > end {
        GRangesError::InvalidGenomicRange(start, end);
    }

    if end >= length {
        GRangesError::InvalidGenomicRangeForSequence(start, end, length);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use super::validate_range;

    #[test]
    fn test_invalid_range_start_end() {
        let range = 10..1;
        let result = validate_range(&range, 10);
        dbg!(&range);
        assert!(matches!(result, Err(GRangesError::InvalidGenomicRange(10, 0))));

    }
}
