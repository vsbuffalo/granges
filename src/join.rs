#![allow(clippy::all)]

use genomap::GenomeMap;

use crate::{Position, granges::GRanges, iterators::GRangesIterator, traits::{RangesIterable, GenericRange}, ranges::coitrees::COITreesIndexed};

pub struct JoinData {
    /// The data index for the left range.
    left: Option<usize>,

    /// A `Vec` of the indices for the overlapping right ranges.
    rights: Vec<Option<usize>>,

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

impl JoinData {
    pub fn new<R: GenericRange>(left_range: &R) -> Self {
        Self {
            left: left_range.index(),
            rights: Vec::new(),
            left_length: left_range.width(),
            right_lengths: Vec::new(),
            overlaps: Vec::new(),
        }
    }
    pub fn add_right<R: GenericRange, Q: GenericRange>(&mut self, left: &R, right: &Q) {
        self.rights.push(right.index());
        self.right_lengths.push(right.width());
        self.overlaps.push(left.overlap_width(right));
    }

}

pub struct JoinIterator<'a, CL, DL, DR> 
where CL: RangesIterable {
          seqnames: Vec<String>,
          left_iter: GRangesIterator<'a, CL>,
          right_genomic_ranges: &'a GenomeMap<COITreesIndexed>,
          left_data: Option<&'a DL>,
          right_data: Option<&'a DR>,
      }

impl<'a, CL, DL: 'a, DR> JoinIterator<'a, CL, DL, DR>
where CL: RangesIterable {
          pub fn new(left: &'a GRanges<CL, DL>, right: &'a GRanges<COITreesIndexed, DR>) -> Self {
              let seqnames = left.seqnames();
              let left_iter = left.iter_ranges();
              let left_data = left.data.as_ref();
              let right_data = right.data.as_ref();
              let right_ranges = &right.ranges;
              Self {
                  seqnames,
                  left_iter,
                  right_genomic_ranges: right_ranges,
                  left_data,
                  right_data,
              }
          }

      }

impl<'a, CL, DL, DR> Iterator for JoinIterator<'a, CL, DL, DR> 
where CL: RangesIterable {
          type Item = JoinData;
          fn next(&mut self) -> Option<Self::Item> {
              if let Some(left_range) = self.left_iter.next() {
                  let seqname = &self.seqnames[left_range.seqname_index];
                  let mut join_data = JoinData::new(&left_range);
                  if let Some(right_ranges) = self.right_genomic_ranges.get(&seqname) {
                      right_ranges.query(left_range.start, left_range.end, |interval| {
                          join_data.add_right(&left_range, interval);
                      });
                  }
                  Some(join_data)
              } else {
                  None
              }
          }
      }


//
//
// pub fn left_join<CL, CR, DL, DR>(left: GRanges<CL, DL>, right: GRanges<CR, DR>) -> JoinIterator<DL, DR> {
//
// }
