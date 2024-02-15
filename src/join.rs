use std::rc::Rc;

use crate::Position;

pub struct JoinData<DL, DR> {
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

    // TODO: we may want some simple summary of whether something is
    // up or downstream. I think the cleanest summary is a signed integer
    // representing the side and degree of non-overlap. E.g. a range
    // that overlaps another but overhangs the 3' side of the focal left
    // range by 10bp is +10; if it were 5', it would be -10.
    /// A possible reference to the left data of type `T`.
    left_data: Option<Rc<DL>>,

    /// A possible reference to the right data of type `T`.
    right_data: Option<Rc<DR>>,
}

pub struct JoinIterator<DL, DR> {
    left_data: Option<Rc<DL>>,
    right_data: Option<Rc<DR>>,
}
//
// impl<DL, DR> JoinIterator<DL, DR> {
//
// }
//
//
// pub fn left_join<CL, CR, DL, DR>(left: GRanges<CL, DL>, right: GRanges<CR, DR>) -> JoinIterator<DL, DR> {
//
// }
