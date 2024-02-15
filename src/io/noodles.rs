//!

use crate::Position;
use noodles::core::Position as NoodlesPosition;

/// Convert from [`noodles::core::Position`], which is [1-based
/// indexing](https://docs.rs/noodles-core/latest/noodles_core/position/struct.Position.html)
/// to GRanges internal 0-based indexing.
pub fn convert_noodles_position(value: NoodlesPosition) -> Position {
    let one_based_position: Position = value.get().try_into().unwrap();
    one_based_position - 1
}
