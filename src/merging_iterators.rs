//! Merging iterators
//!
// TODO: these probably could use better names?
// TODO: we should support weighting by overlaps.
use std::cmp::max;

use crate::{
    error::GRangesError,
    ranges::{GenomicRangeRecord, GenomicRangeRecordEmpty},
    traits::GenericRange,
    PositionOffset,
};

// Create a new [`MergingEmptyResultIterator`], which work over data-less
// "empty" ranges ([`GenomicRangeRecordEmpty`] and merge them based on their
// distance or degree of overlap.
pub struct MergingEmptyResultIterator<I>
where
    I: IntoIterator<Item = Result<GenomicRangeRecordEmpty, GRangesError>>,
{
    last_range: Option<GenomicRangeRecordEmpty>,
    inner: <I as IntoIterator>::IntoIter,
    minimum_distance: PositionOffset,
}

impl<I> MergingEmptyResultIterator<I>
where
    I: IntoIterator<Item = Result<GenomicRangeRecordEmpty, GRangesError>>,
{
    pub fn new(inner: I, minimum_distance: PositionOffset) -> Self {
        Self {
            last_range: None,
            inner: inner.into_iter(),
            minimum_distance,
        }
    }
}

impl<I> Iterator for MergingEmptyResultIterator<I>
where
    I: IntoIterator<Item = Result<GenomicRangeRecordEmpty, GRangesError>>,
{
    type Item = Result<GenomicRangeRecordEmpty, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        for result in self.inner.by_ref() {
            let next_range = match result {
                Ok(range) => range,
                Err(e) => return Some(Err(e)),
            };

            if let Some(last_range) = &mut self.last_range {
                let on_same_chrom = last_range.seqname == next_range.seqname;
                if on_same_chrom
                    && last_range.distance_or_overlap(&next_range) <= self.minimum_distance
                {
                    last_range.end = max(last_range.end, next_range.end);
                } else {
                    let return_range = last_range.clone();
                    self.last_range = Some(next_range);
                    return Some(Ok(return_range));
                }
            } else {
                self.last_range = Some(next_range);
                // If this is the first range, we continue looking for more to potentially merge with.
                continue;
            }
        }

        // If we get here, the inner iterator is exhausted.
        // We need to return the last_range if it exists and make sure it's cleared for subsequent calls.
        self.last_range.take().map(Ok)
    }
}

/// An iterator over [`Result<GenomicRangeRecord<U>, GRangesError>`] that
/// merges ranges that are less than some specified `minimum_distance` apart.
/// If `minimum_distance` is negative, it is taken as a minimum overlap width
/// to merge at.
pub struct MergingResultIterator<I, U, V, F>
where
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
{
    last_range: Option<GenomicRangeRecord<U>>,
    inner: <I as IntoIterator>::IntoIter,
    minimum_distance: PositionOffset,
    func: F,
    accumulated_data: Vec<U>,
}

impl<I, U, V, F> MergingResultIterator<I, U, V, F>
where
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
{
    pub fn new(inner: I, minimum_distance: PositionOffset, func: F) -> Self {
        Self {
            last_range: None,
            inner: inner.into_iter(),
            minimum_distance,
            func,
            accumulated_data: Vec::new(),
        }
    }
}

impl<I, U, V, F> Iterator for MergingResultIterator<I, U, V, F>
where
    GenomicRangeRecord<U>: GenericRange,
    GenomicRangeRecord<V>: GenericRange,
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
{
    type Item = Result<GenomicRangeRecord<V>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        for next_result in self.inner.by_ref() {
            match next_result {
                Ok(next_range) => {
                    if let Some(ref mut last_range) = self.last_range {
                        let on_same_chrom = last_range.seqname == next_range.seqname;
                        if on_same_chrom
                            && last_range.distance_or_overlap(&next_range) <= self.minimum_distance
                        {
                            // this range overlaps the last range, so we keep accumulating data
                            last_range.end = max(last_range.end, next_range.end);
                            self.accumulated_data.push(next_range.data);
                            continue;
                        } else {
                            // New range does not overlap the last range, so we have to finalize.
                            // First, run function on accumulated taken data.
                            let final_data =
                                (self.func)(std::mem::take(&mut self.accumulated_data));
                            let return_range = GenomicRangeRecord {
                                seqname: last_range.seqname.clone(),
                                start: last_range.start,
                                end: last_range.end,
                                data: final_data,
                            };
                            // Next push this new range's data to the data stack (empty after take)
                            self.accumulated_data.push(next_range.data.clone());
                            self.last_range = Some(next_range);
                            // Push the new possibly-extend range with accumulated data.
                            return Some(Ok(return_range));
                        }
                    } else {
                        // There is no last range -- this must be the first.
                        self.last_range = Some(GenomicRangeRecord {
                            seqname: next_range.seqname,
                            start: next_range.start,
                            end: next_range.end,
                            data: next_range.data.clone(),
                        });
                        self.accumulated_data.push(next_range.data);
                        continue;
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        if let Some(last_range) = self.last_range.take() {
            if !self.accumulated_data.is_empty() {
                // Finalize any accumulated data
                let final_data = (self.func)(std::mem::take(&mut self.accumulated_data));
                return Some(Ok(GenomicRangeRecord {
                    seqname: last_range.seqname,
                    start: last_range.start,
                    end: last_range.end,
                    data: final_data,
                }));
            }
        }

        None // Return None when all items have been processed
    }
}

/// A merging iterator (over [`Result<GenomicRangeRecord<U>, GRangesError`] items,
/// i.e. from a parsing iterator) that has an additional merge condition tested
/// by function `G` based on the items.
pub struct ConditionalMergingResultIterator<I, U, V, F, G>
where
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
    G: Fn(&GenomicRangeRecord<U>, &GenomicRangeRecord<U>) -> bool,
{
    last_range: Option<GenomicRangeRecord<U>>,
    inner: <I as IntoIterator>::IntoIter,
    minimum_distance: PositionOffset,
    func: F,
    group: G,
    accumulated_data: Vec<U>,
}

impl<I, U, V, F, G> ConditionalMergingResultIterator<I, U, V, F, G>
where
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
    G: Fn(&GenomicRangeRecord<U>, &GenomicRangeRecord<U>) -> bool,
{
    pub fn new(inner: I, minimum_distance: PositionOffset, func: F, group: G) -> Self {
        Self {
            last_range: None,
            inner: inner.into_iter(),
            minimum_distance,
            func,
            group,
            accumulated_data: Vec::new(),
        }
    }
}

impl<I, U, V, F, G> Iterator for ConditionalMergingResultIterator<I, U, V, F, G>
where
    GenomicRangeRecord<U>: GenericRange,
    GenomicRangeRecord<V>: GenericRange,
    I: IntoIterator<Item = Result<GenomicRangeRecord<U>, GRangesError>>,
    U: Clone,
    V: Clone,
    F: Fn(Vec<U>) -> V,
    G: Fn(&GenomicRangeRecord<U>, &GenomicRangeRecord<U>) -> bool,
{
    type Item = Result<GenomicRangeRecord<V>, GRangesError>;

    fn next(&mut self) -> Option<Self::Item> {
        for next_result in self.inner.by_ref() {
            match next_result {
                Ok(next_range) => {
                    if let Some(ref mut last_range) = self.last_range {
                        // If we have an additional group-by merging condition,
                        // check that. If not set this is always true.
                        let satifies_groupby = (self.group)(last_range, &next_range);

                        let on_same_chrom = last_range.seqname == next_range.seqname;
                        if on_same_chrom
                            && satifies_groupby
                            && last_range.distance_or_overlap(&next_range) <= self.minimum_distance
                        {
                            // this range overlaps the last range, so we keep accumulating data
                            last_range.end = max(last_range.end, next_range.end);
                            self.accumulated_data.push(next_range.data);
                            continue;
                        } else {
                            // New range does not overlap the last range, so we have to finalize.
                            // First, run function on accumulated taken data.
                            let final_data =
                                (self.func)(std::mem::take(&mut self.accumulated_data));
                            let return_range = GenomicRangeRecord {
                                seqname: last_range.seqname.clone(),
                                start: last_range.start,
                                end: last_range.end,
                                data: final_data,
                            };
                            // Next push this new range's data to the data stack (empty after take)
                            self.accumulated_data.push(next_range.data.clone());
                            self.last_range = Some(next_range);
                            // Push the new possibly-extend range with accumulated data.
                            return Some(Ok(return_range));
                        }
                    } else {
                        // There is no last range -- this must be the first.
                        self.last_range = Some(GenomicRangeRecord {
                            seqname: next_range.seqname,
                            start: next_range.start,
                            end: next_range.end,
                            data: next_range.data.clone(),
                        });
                        self.accumulated_data.push(next_range.data);
                        continue;
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        if let Some(last_range) = self.last_range.take() {
            if !self.accumulated_data.is_empty() {
                // Finalize any accumulated data
                let final_data = (self.func)(std::mem::take(&mut self.accumulated_data));
                return Some(Ok(GenomicRangeRecord {
                    seqname: last_range.seqname,
                    start: last_range.start,
                    end: last_range.end,
                    data: final_data,
                }));
            }
        }

        None // Return None when all items have been processed
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        io::{parsers::Bed5Addition, Bed3Iterator, Bed5Iterator},
        merging_iterators::{ConditionalMergingResultIterator, MergingEmptyResultIterator},
        ranges::{GenomicRangeRecord, GenomicRangeRecordEmpty},
    };

    use super::MergingResultIterator;

    #[test]
    fn test_merging_iterators() {
        // NOTE/TODO: the file-based test is because we don't have a good
        // way yet to go from `granges_test_case_01()` into a parsing-like
        // iterator.
        let iter = Bed5Iterator::new("tests_data/test_case_01.bed").unwrap();

        let sum_scores =
            |data: Vec<Bed5Addition>| data.iter().map(|bed5| bed5.score.unwrap()).sum::<f64>();

        let merged_iter = MergingResultIterator::new(iter, 0, sum_scores);

        let results: Vec<_> = Result::from_iter(merged_iter).unwrap();
        // dbg!(&results);

        assert_eq!(
            results[0],
            GenomicRangeRecord::new("chr1".to_string(), 0, 7, 9.2,)
        );

        assert_eq!(
            results[1],
            GenomicRangeRecord::new("chr1".to_string(), 10, 17, 10.1,)
        );

        assert_eq!(
            results[2],
            GenomicRangeRecord::new(
                "chr2".to_string(),
                10,
                32,
                4.800000000000001, // simplifies float comparison
            )
        );
    }

    #[test]
    fn test_merging_empty_iterators() {
        let iter = Bed3Iterator::new("tests_data/test_case_03.bed").unwrap();

        let merged_iter = MergingEmptyResultIterator::new(iter, 0);

        let results: Vec<_> = Result::from_iter(merged_iter).unwrap();

        assert_eq!(
            results[0],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 0, 7,)
        );

        assert_eq!(
            results[1],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 10, 17,)
        );

        assert_eq!(
            results[2],
            GenomicRangeRecordEmpty::new("chr2".to_string(), 10, 32,)
        );
    }

    #[test]
    fn test_merging_empty_iterators_overlap_neg2() {
        let iter = Bed3Iterator::new("tests_data/test_case_03.bed").unwrap();

        // with -2, we require *at least two* overlapping basepairs.
        // chr1 ranges: [0, 5), [4, 7) - these overlap by one, not merged; no others
        // chr2 ranges: [10, 20), [18, 32) - these overlap by two, so merged
        let merged_iter = MergingEmptyResultIterator::new(iter, -2);

        let results: Vec<_> = Result::from_iter(merged_iter).unwrap();

        assert_eq!(
            results[0],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 0, 5,)
        );

        assert_eq!(
            results[1],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 4, 7,)
        );

        assert_eq!(
            results[2],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 10, 17,)
        );

        assert_eq!(
            results[3],
            GenomicRangeRecordEmpty::new("chr2".to_string(), 10, 32,)
        );
    }

    #[test]
    fn test_merging_empty_iterators_overlap_pos10() {
        let iter = Bed3Iterator::new("tests_data/test_case_03.bed").unwrap();

        // with 10, we should just have two range ranges.
        let merged_iter = MergingEmptyResultIterator::new(iter, 10);

        let results: Vec<_> = Result::from_iter(merged_iter).unwrap();

        assert_eq!(
            results[0],
            GenomicRangeRecordEmpty::new("chr1".to_string(), 0, 17)
        );

        assert_eq!(
            results[1],
            GenomicRangeRecordEmpty::new("chr2".to_string(), 10, 32,)
        );
    }

    #[test]
    fn test_conditional_merging_iterators() {
        let iter = Bed5Iterator::new("tests_data/test_case_03.bed").unwrap();

        let sum_scores =
            |data: Vec<Bed5Addition>| data.iter().map(|bed5| bed5.score.unwrap()).sum::<f64>();

        fn group_features(
            range_1: &GenomicRangeRecord<Bed5Addition>,
            range_2: &GenomicRangeRecord<Bed5Addition>,
        ) -> bool {
            range_1.data.name == range_2.data.name
        }

        let merged_iter =
            ConditionalMergingResultIterator::new(iter, 0, sum_scores, group_features);

        let results: Vec<_> = Result::from_iter(merged_iter).unwrap();

        // unlike unconditional merging iterators, the first ranges here,
        // though book-ended, are *not* overlapped because they have different feature
        // names
        assert_eq!(
            results[0],
            GenomicRangeRecord::new("chr1".to_string(), 0, 5, 1.1,)
        );

        assert_eq!(
            results[1],
            GenomicRangeRecord::new("chr1".to_string(), 4, 7, 8.1,)
        );

        assert_eq!(
            results[2],
            GenomicRangeRecord::new("chr1".to_string(), 10, 17, 10.1,)
        );

        // these *are* merged because they have the same feature name
        assert_eq!(
            results[3],
            GenomicRangeRecord::new(
                "chr2".to_string(),
                10,
                32,
                4.800000000000001, // simplifies float comparison
            )
        );
    }
}
