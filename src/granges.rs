use genomap::GenomeMap;

use crate::{traits::RangeContainer, ranges::{vec::{VecRanges, VecRangesIndexed, VecRangesEmpty}, RangeIndexed, RangeEmpty}, Position};


pub struct GRanges<C, T> {
    ranges: GenomeMap<C>,
    data: Option<T>,
}


impl<C, T> GRanges<C, T>
where C: RangeContainer {
    
    /// Get the total number of ranges.
    pub fn len(&self) -> usize {
        self.ranges.values().map(|ranges| ranges.len()).sum()
    }

    /// Return whether the [`GRanges`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<U> GRanges<VecRangesIndexed, Vec<U>> {

    /// Create a new [`GRanges`] object, with vector storage for ranges and data.
    ///
    /// This combination of range and data containers is used when loading data into
    /// a new [`GRanges`] object, and the size cannot be known beforehand. Rust's 
    /// [`Vec`] will dynamically grow to accommodate new ranges; use [`GRanges.shrink()`] 
    /// call the [`Vec`]'s shrink to size methods on the range and data containers
    /// after data loading to shrink to the minimal necessary size (this can reduce
    /// memory usage).
    pub fn new_vec() -> Self {
        let ranges = GenomeMap::new();
        Self {
            ranges,
            data: None,
        }
    }


    pub fn push_range_with_data(&mut self, seqname: &str, start: Position, end: Position, data: U) {
        // push data to the vec data container, getting the index
        let index: usize = {
            let data_container = self.data.get_or_insert_with(Vec::new);
            data_container.push(data);
            data_container.len() - 1 // new data index
        };
        // push an indexed range
        let range = RangeIndexed::new(start, end, index);
        self.ranges.entry_or_default(seqname).ranges.push(range);
    }
}

impl GRanges<VecRangesEmpty, ()> {

    /// Create a new [`GRanges`] object, with vector storage for ranges and no data container.
    pub fn new_vec_empty() -> Self {
        let ranges = GenomeMap::new();
        Self {
            ranges,
            data: None,
        }
    }

    /// Push an empty range (no data) to the [`VecRangesEmpty`] range container.
    pub fn push_range_only(&mut self, seqname: &str, start: Position, end: Position) {
        // push an unindexed (empty) range
        let range = RangeEmpty::new(start, end);
        self.ranges.entry_or_default(seqname).ranges.push(range);
    }
}



#[cfg(test)]
mod tests {
    use crate::{prelude::*, test_utilities::random_vecranges};

    #[test]
    fn test_new_vec() {
        let mut gr = GRanges::new_vec();
        gr.push_range_with_data("chr1", 0, 10, 1.1);
        assert_eq!(gr.len(), 1);
    }

    #[test]
    fn test_random_vecranges() {
        let vr = random_vecranges(100);
        assert_eq!(vr.len(), 100)
    }

}
