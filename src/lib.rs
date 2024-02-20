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
    pub use crate::io::{Bed3Iterator, BedlikeIterator, GenomicRangesFile, TsvRecordIterator, GenomicRangesParser};

    pub use crate::ranges::vec::{VecRangesEmpty, VecRangesIndexed};
    pub use crate::traits::{
        GeneralRangeRecordIterator, GenericRange, GenomicRangesTsvSerialize, IndexedDataContainer,
        IntoIterableRangesContainer, IterableRangeContainer, TsvSerialize,
    };

    pub use crate::seqlens;
}

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
