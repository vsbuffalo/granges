//! Range Operations.
//!
//! - [`adjust()`]: Adjust range start and end positions.

use crate::{PositionOffset, Position};

use super::RangeIndexed;

/// Adjusts the coordinates of a range, ensuring the adjusted range is within [0, length]
/// and returning `None` if the range has zero width after adjustment.
pub fn adjust(range: RangeIndexed, start_delta: PositionOffset, end_delta: PositionOffset, length: Position) -> Option<RangeIndexed> {
    let start: PositionOffset = range.start.try_into().unwrap();
    let end: PositionOffset = range.end.try_into().unwrap();
    let length: PositionOffset = length.try_into().unwrap();

    // ensure within [0, length]
    let new_start = (start + start_delta).max(0).min(length);
    // ensure new_end >= new_start and within [0, length]
    let new_end = (end + end_delta).max(new_start).min(length); 

    // check for zero-width range
    if new_end <= new_start {
        // return None if the range has zero width
        None
    } else {
        Some(RangeIndexed::new(new_start.try_into().unwrap(), new_end.try_into().unwrap(), range.index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_adjustment() {
        let range = RangeIndexed::new(5, 10, 1);
        let adjusted = adjust(range, -2, 3, 15).unwrap();
        assert_eq!(adjusted, RangeIndexed::new(3, 13, 1));
    }

    #[test]
    fn test_out_of_bounds_adjustment() {
        let range = RangeIndexed::new(10, 12, 2);
        let adjusted = adjust(range, -5, 20, 15).unwrap();
        assert_eq!(adjusted, RangeIndexed::new(5, 15, 2));
    }
    
    #[test]
    fn test_zero_width_result() {
        let range = RangeIndexed::new(5, 10, 3);
        assert!(adjust(range, 5, -5, 15).is_none());
    }
}

