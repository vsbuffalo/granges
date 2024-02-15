
pub mod traits;
pub mod ranges;
pub mod granges;
pub mod error;
pub mod test_utilities;

pub type Position = u32;

pub mod prelude {
    pub use crate::granges::GRanges;
    pub use crate::error::GRangesError;
}
