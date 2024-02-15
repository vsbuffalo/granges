use thiserror::Error;

use crate::Position;

#[derive(Debug, Error)]
pub enum GRangesError {
    #[error("Range invalid: start ({0}) must be greater than end ({1})")]
    InvalidGenomicRange(Position, Position),

    #[error("Range [{0}, {1}] is invalid for sequence of length {2}")]
    InvalidGenomicRangeForSequence(Position, Position, Position),

}
