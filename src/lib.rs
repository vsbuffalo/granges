// Copyright (2024) Vince Buffalo
#![crate_name = "granges"]
#![doc(html_root_url = "https://docs.rs/granges/")]

pub use indexmap;

pub mod data;
pub mod error;
pub mod granges;
pub mod io;
pub mod iterators;
pub mod join;
pub mod ranges;
pub mod test_utilities;
pub mod traits;

// bringing these CLI modules into lib.rs rather than main/ allows for
// use in integration tests and other Rust-side command line work
pub mod commands;
pub mod reporting;

pub type Position = u32;
pub type PositionOffset = i32; // signed variant

pub mod prelude {
    pub use crate::error::GRangesError;
    pub use crate::granges::{GRanges, GRangesEmpty};
    pub use crate::io::file::read_seqlens;
    pub use crate::io::{
        Bed3Iterator, BedlikeIterator, GenomicRangesFile, GenomicRangesParser, TsvRecordIterator,
    };

    pub use crate::ranges::vec::{VecRangesEmpty, VecRangesIndexed};
    pub use crate::traits::{
        AsGRangesRef, GeneralRangeRecordIterator, GenericRange, GenericRangeOperations,
        GenomicRangeRecordUnwrappable, GenomicRangesTsvSerialize, IndexedDataContainer,
        IntoIterableRangesContainer, IterableRangeContainer, TsvSerialize,
    };

    pub use crate::seqlens;
}

pub const INTERNAL_ERROR_MESSAGE: &str = r#"
An internal error has occurred. Please file a GitHub issue at 
https://github.com/vsbuffalo/granges/issues
"#;

/// Create and [`indexmap::IndexMap`] of sequence names and their lengths.
#[macro_export]
macro_rules! seqlens {
    ($($key:expr => $value:expr),* $(,)?) => {
        $crate::indexmap::indexmap!($($key.to_string() => $value),*)
    };
}

/// Create a new `GRanges<R, T>` with sequence length information (used primarily for small examples)
#[macro_export]
macro_rules! create_granges_with_seqlens {
    ($range_type:ty, $data_type:ty, { $($chr:expr => [$(($start:expr, $end:expr, $data:expr)),*]),* }, seqlens: { $($chr_len:expr => $len:expr),* }) => {
        {
            let mut seqlens = ::indexmap::IndexMap::new();
            $(seqlens.insert($chr_len.to_string(), $len);)*

            let mut gr: GRanges<$range_type, $data_type> = GRanges::new_vec(&seqlens);

            $(
                $(
                    gr.push_range(&$chr.to_string(), $start, $end, $data).expect("Failed to push range");
                )*
            )*

            gr
        }
    };
}

/// a version of assert_eq! that is more user-friendly.
#[macro_export]
macro_rules! ensure_eq {
    ($left:expr, $right:expr $(,)?) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                if !(*left_val == *right_val) {
                    panic!("{}\nExpected `{}` but found `{}`.", crate::INTERNAL_ERROR_MESSAGE, stringify!($left), stringify!($right));

                }
            }
        }
    });
    ($left:expr, $right:expr, $($arg:tt)+) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                if !(*left_val == *right_val) {
                    panic!("{}\n{}\nExpected `{}` but found `{}`.", crate::INTERNAL_ERROR_MESSAGE, format_args!($($arg)+), stringify!($left), stringify!($right));
                }
            }
        }
    });
}

#[cfg(test)]
mod tests {
    #[test]
    #[should_panic]
    fn test_ensure_eq() {
        ensure_eq!(1, 2);
    }

    #[test]
    #[should_panic]
    fn test_ensure_eq_msg() {
        ensure_eq!(1, 2, "Some description of the problem.");
    }

    #[test]
    fn test_ensure_eq_pass() {
        ensure_eq!(1, 1);
    }
}
