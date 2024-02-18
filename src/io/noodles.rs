//!

use crate::Position;
use noodles::core::Position as NoodlesPosition;

/// Convert from [`noodles::core::Position`], which is [1-based
/// indexing](https://docs.rs/noodles-core/latest/noodles_core/position/struct.Position.html)
/// and right-to GRanges internal 0-based indexing.
///
/// # Developers Note
/// See: https://github.com/zaeleus/noodles/issues/226
pub fn convert_noodles_range(start: NoodlesPosition, end: NoodlesPosition) -> (Position, Position) {
    let start_one_based: Position = start.get().try_into().unwrap();
    let end_one_based: Position = end.get().try_into().unwrap();
    (start_one_based - 1, end_one_based)
}
