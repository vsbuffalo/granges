#![allow(clippy::all)]

use crate::{Position, granges::GRanges, iterators::GRangesIterator, traits::RangesIterable};

pub struct JoinData {
    /// The data index for the left range.
    left: usize,

    /// A `Vec` of the indices for the overlapping right ranges.
    rights: Vec<usize>,

    /// The length of the left range.
    left_length: Position,

    /// The lengths of the right ranges.
    right_lengths: Vec<Position>,

    /// The lengths of the overlaps between the left and right ranges.
    overlaps: Vec<Position>,

    // TODO: we may want some simple summary of whether an overlapping range is
    // up or downstream. I think the cleanest summary is a signed integer
    // representing the side and degree of non-overlap. E.g. a range
    // that overlaps another but overhangs the 3' side of the focal left
    // range by 10bp is +10; if it were 5', it would be -10.
}

pub struct JoinIterator<'a, C, DL, DR> 
where C: RangesIterable {
    left_iter: GRangesIterator<'a, C>,
    left_data: Option<DL>,
    right_data: Option<DR>,
}

impl<'a, C, DL, DR> JoinIterator<'a, C, DL, DR>
where C: RangesIterable {
    pub fn new<CL, CR>(left: GRanges<CL, DL>, right: &'a GRanges<CL, DR>) {
    }
}

impl<'a, C, DL, DR> Iterator for JoinIterator<'a, C, DL, DR> 
where C: RangesIterable {
    type Item = JoinData;
    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}


//
//
// pub fn left_join<CL, CR, DL, DR>(left: GRanges<CL, DL>, right: GRanges<CR, DR>) -> JoinIterator<DL, DR> {
//
// }
