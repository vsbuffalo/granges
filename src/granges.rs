use std::path::PathBuf;

use genomap::GenomeMap;
use indexmap::IndexMap;

use crate::{
    io::OutputFile,
    iterators::GRangesIterator,
    join::JoinIterator,
    prelude::GRangesError,
    ranges::{
        coitrees::{COITrees, COITreesIndexed},
        vec::{VecRanges, VecRangesEmpty, VecRangesIndexed},
        GenomicRangeRecord, RangeEmpty, RangeIndexed,
    },
    traits::{GenericRange, IndexedDataContainer, RangeContainer, RangesIterable, TsvSerialize},
    Position,
};

#[derive(Clone, Debug)]
pub struct GRanges<C, T> {
    pub(crate) ranges: GenomeMap<C>,
    pub(crate) data: Option<T>,
}

impl<C, T> GRanges<C, T>
where
    C: RangeContainer,
{
    /// Get the total number of ranges.
    pub fn len(&self) -> usize {
        self.ranges.values().map(|ranges| ranges.len()).sum()
    }

    /// Return whether the [`GRanges`] object is empty (contains no ranges).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the raw range container.
    pub fn get_ranges(&self, seqname: &str) -> Option<&C> {
        self.ranges.get(seqname)
    }

    /// Get the sequence names.
    pub fn seqnames(&self) -> Vec<String> {
        self.ranges.names()
    }

    /// Get the sequences lengths.
    pub fn seqlens(&self) -> IndexMap<String, Position> {
        let seqlens = self
            .ranges
            .iter()
            .map(|(seqname, ranges)| (seqname.to_string(), ranges.sequence_length()))
            .collect();
        seqlens
    }
}

impl<R: GenericRange, T> GRanges<VecRanges<R>, T> {
    /// Create a new [`GRanges`] object, with vector storage for ranges and data.
    ///
    /// This combination of range and data containers is used when loading data into
    /// a new [`GRanges`] object, and the size cannot be known beforehand. Rust's
    /// [`Vec`] will dynamically grow to accommodate new ranges; use [`GRanges.shrink()`]
    /// call the [`Vec`]'s shrink to size methods on the range and data containers
    /// after data loading to shrink to the minimal necessary size (this can reduce
    /// memory usage).
    pub fn new_vec(seqlens: &IndexMap<String, Position>) -> Self {
        let mut ranges = GenomeMap::new();
        for (seqname, length) in seqlens.iter() {
            // this should never happen because the error is only if
            // insert encounters a seqname that's already been inserted -- that
            // cannot happen here.
            ranges
                .insert(seqname, VecRanges::new(*length))
                .expect("Internal error: please report");
        }
        Self { ranges, data: None }
    }

    /// Consume this [`GRanges`] object and sort the ranges.
    pub fn sort(mut self) -> Self {
        self.ranges.values_mut().for_each(|ranges| ranges.sort());
        self
    }
}

impl<U> GRanges<VecRangesIndexed, Vec<U>> {
    /// Push a genomic range with its data to the range and data containers in a [`GRanges] object.
    pub fn push_range_with_data(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
        data: U,
    ) -> Result<(), GRangesError> {
        // push data to the vec data container, getting the index
        let index: usize = {
            let data_container = self.data.get_or_insert_with(Vec::new);
            data_container.push(data);
            data_container.len() - 1 // new data index
        };
        // push an indexed range
        let range = RangeIndexed::new(start, end, index);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<'a, T> GRanges<VecRanges<RangeIndexed>, T>
where
    T: IndexedDataContainer<'a>,
    T: TsvSerialize,
    <T as IndexedDataContainer<'a>>::Item: TsvSerialize,
{
    ///
    pub fn to_tsv(&'a self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError> {
        // output stream -- header is None for now (TODO)
        let output = output.map_or(OutputFile::new_stdout(None), |file| {
            OutputFile::new(file, None)
        });
        let mut writer = output.writer()?;

        let seqnames = self.seqnames();
        for range in self.iter_ranges() {
            let record = range.to_record(&seqnames, self.data.as_ref());
            writeln!(writer, "{}", record.to_tsv())?;
        }
        Ok(())
    }
}

impl<'a, C, T> GRanges<C, T>
where
    T: IndexedDataContainer<'a>,
{
    /// Get the data in the data container at specified index.
    ///
    /// # Panics
    /// This will panic if there if the index is invalid, or the
    /// data container is `None`. Both of these indicate internal
    /// design errors: please file an issue of you encounter a panic.
    pub fn get_data_value(&'a self, index: usize) -> <T as IndexedDataContainer>::Item {
        self.data
            .as_ref()
            .expect("data container was None")
            .get_value(index)
    }
}

impl<T> GRanges<VecRanges<RangeEmpty>, T> {
    /// Push an empty range (no data) to the [`VecRangesEmpty`] range container.
    pub fn push_range(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<(), GRangesError> {
        // push an unindexed (empty) range
        let range = RangeEmpty::new(start, end);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<T> GRanges<VecRanges<RangeIndexed>, T> {
    /// Push an empty range (no data) to the [`VecRangesEmpty`] range container.
    pub fn push_range_with_index(
        &mut self,
        seqname: &str,
        start: Position,
        end: Position,
        index: usize,
    ) -> Result<(), GRangesError> {
        // push an unindexed (empty) range
        let range = RangeIndexed::new(start, end, index);
        let range_container = self
            .ranges
            .get_mut(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        range_container.push_range(range);
        Ok(())
    }
}

impl<U> GRanges<VecRangesIndexed, Vec<U>> {
    pub fn from_iter<I>(
        iter: I,
        seqlens: &IndexMap<String, Position>,
    ) -> Result<GRanges<VecRangesIndexed, Vec<U>>, GRangesError>
    where
        I: Iterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    {
        let mut gr = GRanges::new_vec(&seqlens);
        for possible_entry in iter {
            let entry = possible_entry?;
            gr.push_range_with_data(&entry.seqname, entry.start, entry.end, entry.data)?;
        }
        Ok(gr)
    }
}

impl GRanges<VecRangesEmpty, ()> {
    pub fn from_iter_empty<I>(
        iter: I,
        seqlens: IndexMap<String, Position>,
    ) -> Result<GRanges<VecRangesEmpty, ()>, GRangesError>
    where
        I: Iterator<Item = Result<GenomicRangeRecord<()>, GRangesError>>,
    {
        let mut gr = GRanges::new_vec(&seqlens);
        for possible_entry in iter {
            let entry = possible_entry?;
            gr.push_range(&entry.seqname, entry.start, entry.end)?;
        }
        Ok(gr)
    }
}

impl<R> GRanges<R, ()>
where
    R: RangeContainer + RangesIterable,
{
    // TODO: candidate for a trait
    pub fn to_bed3(&self, output: Option<impl Into<PathBuf>>) -> Result<(), GRangesError> {
        // output stream -- header is None for now (TODO)
        let output = output.map_or(OutputFile::new_stdout(None), |file| {
            OutputFile::new(file, None)
        });
        let mut writer = output.writer()?;

        let seqnames = self.seqnames();
        for range in self.iter_ranges() {
            let record = range.to_record_empty::<()>(&seqnames);
            writeln!(writer, "{}", record.to_tsv())?;
        }
        Ok(())
    }
}

impl<T> GRanges<VecRangesIndexed, T> {
    /// Convert this [`VecRangesIndexed`] range container to a cache-oblivious interval tree  
    /// range container, [`COITreesIndexed`]. This is done using the [`coitrees`] library
    /// by Daniel C. Jones.
    pub fn to_coitrees(self) -> Result<GRanges<COITreesIndexed, T>, GRangesError> {
        let old_ranges = self.ranges;
        let mut new_ranges = GenomeMap::new();
        for (seqname, vec_ranges) in old_ranges.into_iter() {
            let trees = COITrees::from(vec_ranges);
            new_ranges.insert(&seqname, trees)?;
        }
        Ok(GRanges {
            ranges: new_ranges,
            data: self.data,
        })
    }
}

impl<'a, CL, U> GRanges<CL, Vec<U>>
where
    CL: RangesIterable,
    <CL as RangesIterable>::RangeType: GenericRange,
    U: Clone,
{
    //pub fn left_overlaps<DR>(self, right: &'a GRanges<COITreesIndexed, DR>)
    //-> GRanges<CL, JoinIterator<'a, CL, Vec<U>, DR>> {
    //    //let mut obj = GRanges {
    //    //    ranges: self.ranges,
    //    //    data: None,
    //    //};
    //    //obj.data = Some(JoinIterator::new(&obj, &right));
    //    //obj
    //    todo!()
    //}

    /// Filter out ranges that do *not* have at least overlap with the `right` ranges.
    ///
    /// In database lingo, this is a type of *filtering join*, in particular a *semi join*.
    /// See Hadley Wickham's excellent [R for Data
    /// Science](https://r4ds.hadley.nz/joins.html#filtering-joins) for more information.
    ///
    /// Note that this consumes the `self` [`GRanges`] object, turning it into a new
    /// [`GRanges<VecRangesIndexed, Vec<U>`]. The data container is rebuilt from indices
    /// into a new [`Vec<U>`] where `U` is the associated type [`IndexedDataContainer::Item`],
    /// which represents the individual data element in the data container.
    pub fn filter_overlaps<DR: IndexedDataContainer<'a>>(
        self,
        right: &GRanges<COITreesIndexed, DR>,
    ) -> Result<GRanges<VecRangesIndexed, Vec<U>>, GRangesError> {
        let mut gr: GRanges<VecRangesIndexed, Vec<U>> = GRanges::new_vec(&self.seqlens());

        for (seqname, left_ranges) in self.ranges.iter() {
            for left_range in left_ranges.iter_ranges() {
                if let Some(right_ranges) = right.ranges.get(&seqname) {
                    let num_overlaps =
                        right_ranges.count_overlaps(left_range.start(), left_range.end());
                    if num_overlaps == 0 {
                        // no overlaps -- skip
                    } else {
                        gr.push_range_with_index(
                            seqname,
                            left_range.start(),
                            left_range.end(),
                            left_range.index().unwrap(),
                        )?;
                    }
                }
            }
        }
        Ok(gr)
    }

    pub fn filter_overlaps_anti<DR: IndexedDataContainer<'a>>(&self, right: &GRanges<COITreesIndexed, DR>)
        -> Result<GRanges<VecRangesIndexed, Vec<U>>, GRangesError> {
            todo!()
        }
}

impl<R, T> GRanges<R, T>
where
    R: RangesIterable,
{
    /// Create a new [`GRangesIterator`] to iterate through all the ranges in this [`GRanges`] object.
    pub fn iter_ranges(&self) -> GRangesIterator<'_, R> {
        GRangesIterator::new(&self.ranges)
    }
}

#[cfg(test)]
mod tests {
    use indexmap::indexmap;

    use crate::{
        prelude::*,
        test_utilities::{granges_test_case_01, random_vecranges},
    };

    #[test]
    fn test_new_vec() {
        let seqlens = indexmap! { "chr1".to_string() => 10};
        let mut gr = GRanges::new_vec(&seqlens);
        gr.push_range_with_data("chr1", 0, 10, 1.1).unwrap();
        assert_eq!(gr.len(), 1);
    }

    #[test]
    fn test_random_vecranges() {
        let vr = random_vecranges(100);
        assert_eq!(vr.len(), 100)
    }

    #[test]
    fn test_to_coitrees() {
        let gr_vec = granges_test_case_01();
        let gr = gr_vec.clone().to_coitrees().unwrap();
        assert_eq!(gr.len(), 5);
    }
}
