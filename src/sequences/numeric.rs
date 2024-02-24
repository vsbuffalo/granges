//! Per-basepair numeric sequence data.
//!
//! If the crate `ndarray` feature is off, this module will not be loaded, since it requires
//! [ndarray](https://crates.io/crates/ndarray).
//!
use std::{cell::Ref, path::Path};

use genomap::GenomeMap;
use indexmap::IndexMap;
#[cfg(feature = "ndarray")]
use ndarray::{s, Array1, Array2, ArrayView1, ArrayView2};
use ndarray_npy::ReadableElement;
use std::path::PathBuf;

use super::lazy::LazyLoader;
use crate::error::GRangesError;
use crate::traits::Sequences;
use crate::Position;

/// A one-dimensional generic numeric per-basepair container.
///
/// This struct holds one-dimensional sequence data where each basepair position
/// is associated with a numeric value.
///
/// # Examples
///
/// ```
/// use granges::sequences::numeric::NumericSequences1;
/// use granges::test_utilities::random_array1_sequences;
///
/// let mut data = random_array1_sequences(100);
/// let numeric_seq = NumericSequences1::new(data);
/// ```
pub struct NumericSequences1<T> {
    data: GenomeMap<Array1<T>>,
}

impl<T> NumericSequences1<T> {
    /// Create a new in-memory per-basepair `Array1<T>` sequence data container.
    ///
    /// # Arguments
    ///
    /// * `data` - An [`SequenceMap`] mapping sequence names to `Array1<T>` data.
    ///
    /// # Examples
    ///
    /// ```
    /// use ndarray::Array1;
    /// use granges::sequences::numeric::NumericSequences1;
    /// use granges::test_utilities::random_array1_sequences;
    ///
    /// let mut data = random_array1_sequences(100);
    /// let numeric_seq = NumericSequences1::new(data);
    /// ```
    pub fn new(data: GenomeMap<Array1<T>>) -> Self {
        Self { data }
    }
}

impl<'a, T: 'a> Sequences<'a> for NumericSequences1<T> {
    type Container = &'a Array1<T>;
    type Slice = ArrayView1<'a, T>;

    fn seqnames(&self) -> Vec<String> {
        self.data.names()
    }

    fn get_sequence(&'a self, seqname: &str) -> Result<Self::Container, GRangesError> {
        let seq = self
            .data
            .get(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        Ok(seq)
    }

    /// Apply a function to `ndarray` views of the `Array1<T>` sequence data at the specified
    /// range.
    ///
    /// # Arguments
    fn region_apply<V, F>(
        &'a self,
        func: F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: for<'b> Fn(Self::Slice) -> V,
    {
        let seq = self.get_sequence(seqname)?;
        let start_usize: usize = start.try_into().unwrap();
        let end_usize: usize = end.try_into().unwrap();
        let view = seq.slice(s![start_usize..end_usize]);
        Ok(func(view))
    }

    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError> {
        let seq = self.get_sequence(seqname)?;
        Ok(seq.len().try_into().unwrap())
    }
}

/// A two-dimensional generic numeric per-basepair container.
pub struct NumericSequences2<T> {
    data: GenomeMap<Array2<T>>,
}

impl<T> NumericSequences2<T> {
    /// Create a new in-memory per-basepair [`Array2<T>`] sequence data container.
    ///
    /// # Arguments
    ///
    /// * `data` - An [`SequenceMap`] mapping sequence names to `Array2<T>` data.
    ///
    /// # Examples
    ///
    /// ```
    /// use granges::sequences::numeric::NumericSequences2;
    /// use granges::test_utilities::random_array2_sequences;
    ///
    /// let data = random_array2_sequences(20);
    /// let numeric_seq = NumericSequences2::new(data);
    /// ```
    /// [`Array2<T>`]: https://docs.rs/ndarray/latest/ndarray/type.Array2.html
    pub fn new(data: GenomeMap<Array2<T>>) -> Self {
        Self { data }
    }
}

impl<'a, T: 'a> Sequences<'a> for NumericSequences2<T> {
    type Container = &'a Array2<T>;
    type Slice = ArrayView2<'a, T>;

    /// Retrieves the names of all sequences stored in the container.
    ///
    /// # Returns
    ///
    /// A vector of strings, each representing a sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::granges::traits::Sequences;
    /// use granges::sequences::numeric::NumericSequences1;
    /// use granges::test_utilities::random_array1_sequences;
    ///
    /// let data = random_array1_sequences(20);
    /// let numeric_seq = NumericSequences1::new(data);
    /// let seq_names = numeric_seq.seqnames();
    /// ```
    fn seqnames(&self) -> Vec<String> {
        self.data.names()
    }

    /// Gets a reference to the sequence data associated with a given sequence name.
    ///
    /// # Arguments
    ///
    /// * `seqname` - The name of the sequence to retrieve.
    ///
    /// # Returns
    ///
    /// A result containing a reference to the sequence data (`Array1<T>`) on success,
    /// or a `GRangesError` on failure.
    ///
    /// # Examples
    ///
    /// ```
    /// use granges::sequences::numeric::NumericSequences1;
    /// use crate::granges::traits::Sequences;
    /// use granges::test_utilities::random_array1_sequences;
    ///
    /// let mut data = random_array1_sequences(100);
    /// let numeric_seq = NumericSequences1::new(data);
    /// let sequence = numeric_seq.get_sequence("chr1").unwrap();
    /// ```
    fn get_sequence(&'a self, seqname: &str) -> Result<Self::Container, GRangesError> {
        let seq = self
            .data
            .get(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;
        Ok(seq)
    }

    /// Applies a given function to a slice of the sequence data for a specified range.
    ///
    /// # Arguments
    ///
    /// * `func` - A closure or function pointer that takes an `ArrayView2<T>` slice of sequence data and returns a value.
    /// * `seqname` - The name of the sequence.
    /// * `start` - The start position (inclusive) of the range.
    /// * `end` - The end position (inclusive) of the range.
    ///
    /// # Returns
    ///
    /// A result containing the value returned by the function `func`, or a `GRangesError` on failure.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::granges::traits::Sequences;
    /// use granges::sequences::numeric::NumericSequences1;
    /// use granges::test_utilities::random_array1_sequences;
    ///
    ///
    /// let mut data = random_array1_sequences(100);
    /// let numeric_seq = NumericSequences1::new(data);
    /// let result = numeric_seq.region_apply(
    ///     |view| view.sum(),
    ///     "chr1",
    ///     0,
    ///     10
    /// ).unwrap();
    /// ```
    fn region_apply<V, F>(
        &'a self,
        func: F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: for<'b> Fn(Self::Slice) -> V,
    {
        let seq = self.get_sequence(seqname)?;
        let start_usize: usize = start.try_into().unwrap();
        let end_usize: usize = end.try_into().unwrap();
        let view = seq.slice(s![start_usize..end_usize, ..]);
        Ok(func(view))
    }

    /// Retrieves the length of a specific sequence.
    ///
    /// # Arguments
    ///
    /// * `seqname` - The name of the sequence.
    ///
    /// # Returns
    ///
    /// A result containing the length of the sequence as an [`Position`], or a
    /// [`GRangesError`] on failure.
    ///
    /// # Examples
    ///
    /// ```
    /// use granges::sequences::numeric::NumericSequences1;
    /// use crate::granges::traits::Sequences;
    /// use granges::test_utilities::random_array1_sequences;
    ///
    /// let data = random_array1_sequences(20);
    /// let numeric_seq = NumericSequences1::new(data);
    ///
    /// let seqnames = numeric_seq.seqnames();
    /// let first_seqname = seqnames.first().unwrap();
    /// let length = numeric_seq.get_sequence_length(first_seqname).unwrap();
    /// ```
    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError> {
        let seq = self.get_sequence(seqname)?;
        Ok(seq.len().try_into().unwrap())
    }
}

/// A lazy-loaded two-dimensional generic numeric per-basepair container.
pub struct LazyNumericSequences2<T>
where
    T: std::fmt::Debug,
{
    seqlens: IndexMap<String, Position>,
    lazy: LazyLoader<Option<()>, Array2<T>, String>,
}

impl<T> LazyNumericSequences2<T>
where
    T: ReadableElement + std::fmt::Debug,
{
    /// Create a link to an out-of-memory per-basepair [`Array2<T>`] sequence data container.
    ///
    /// For `.npy` files, use [`ndarray_npy::read_npy`] for the `loader`.
    ///
    /// # Arguments
    ///
    /// * `dir` - the directory containing the files.
    /// * `pattern` - the file pattern name with a `seqname` wildcard inside `dir`, e.g.
    ///                "numeric_{}.npy" where `{}` would be a name in `seqlen`.
    ///
    pub fn new<F>(dir: &str, pattern: &str, loader: F, seqlens: IndexMap<String, Position>) -> Self
    where
        F: Fn(PathBuf) -> Result<Array2<T>, GRangesError> + 'static,
    {
        let dir_path = Path::new(dir).to_owned();
        let pattern_owned = pattern.to_owned();

        let lazy = LazyLoader::new(None, move |_reader, seqname: &String| {
            let filename = pattern_owned.replace("{}", seqname);
            let filepath = dir_path.join(filename);
            let data: Array2<T> = loader(filepath)?;
            Ok(data)
        });

        Self { seqlens, lazy }
    }
}

impl<'a, T: 'a> Sequences<'a> for LazyNumericSequences2<T>
where
    T: std::fmt::Debug,
{
    type Container = Ref<'a, Array2<T>>;
    type Slice = ArrayView2<'a, T>;

    /// Retrieves the names of all sequences stored in the container.
    ///
    /// # Returns
    ///
    /// A vector of strings, each representing a sequence name.
    ///
    fn seqnames(&self) -> Vec<String> {
        self.seqlens.keys().cloned().collect()
    }

    /// Gets a reference to the sequence data associated with a given sequence name.
    ///
    /// # Arguments
    ///
    /// * `seqname` - The name of the sequence to retrieve.
    ///
    /// # Returns
    ///
    /// A result containing a reference to the sequence data (`Array1<T>`) on success,
    /// or a `GRangesError` on failure.
    ///
    fn get_sequence(&'a self, seqname: &str) -> Result<Self::Container, GRangesError> {
        self.lazy.get_data(&seqname.to_string())
    }

    /// Applies a given function to a slice of the sequence data for a specified range.
    ///
    /// # Arguments
    ///
    /// * `func` - A closure or function pointer that takes an `ArrayView2<T>` slice of sequence data and returns a value.
    /// * `seqname` - The name of the sequence.
    /// * `start` - The start position (inclusive) of the range.
    /// * `end` - The end position (inclusive) of the range.
    ///
    /// # Returns
    ///
    /// A result containing the value returned by the function `func`, or a `GRangesError` on failure.
    ///
    fn region_apply<V, F>(
        &'a self,
        func: F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: for<'b> Fn(ArrayView2<'b, T>) -> V,
    {
        let seq = self.get_sequence(seqname)?;
        let start_usize: usize = start.try_into().unwrap();
        let end_usize: usize = end.try_into().unwrap();
        let view = seq.slice(s![start_usize..end_usize, ..]);
        let value = func(view);
        Ok(value)
    }

    /// Retrieves the length of a specific sequence.
    ///
    /// # Arguments
    ///
    /// * `seqname` - The name of the sequence.
    ///
    /// # Returns
    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError> {
        let seq = self.get_sequence(seqname)?;
        let seqlen = seq.len().try_into().unwrap();
        Ok(seqlen)
    }
}

#[cfg(test)]
mod tests {}
