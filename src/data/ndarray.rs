//! Data container implementations for [`ndarray::Array1`] and [`ndarray::Array2`].

use crate::error::GRangesError;
use crate::granges::GRanges;
use crate::traits::{DataContainer, IndexedDataContainer, RangeContainer};
use ndarray::{Array1, Array2, ArrayView1};

impl<T> DataContainer for Array1<T> {}
impl<T> DataContainer for Array2<T> {}

impl<C, T> GRanges<C, T>
where
    C: RangeContainer,
{
    /// Take the container data into an iterator, apply a function, a put the results in an
    /// [`Array1`].
    ///
    /// This is a convenience function for `GRanges::map_data().into_array1()`.
    ///
    /// [`Array1`]: ndarray::Array1
    pub fn map_into_array1<F, U>(mut self, func: F) -> Result<GRanges<C, Array1<U>>, GRangesError>
    where
        F: FnMut(<T as IntoIterator>::Item) -> U,
        T: IntoIterator,
    {
        let total_len = self.len();
        let (ranges, data) = self.take_both()?;

        let transformed_data = data.into_iter().map(func).collect::<Array1<U>>();

        assert_eq!(
            transformed_data.len(),
            total_len,
            "Transformed data has different length than GRanges!"
        );

        Ok(GRanges {
            ranges,
            data: Some(transformed_data),
        })
    }

    /// Take the container data into an iterator, apply a function, a put the results in an
    /// [`Array2`].
    ///
    /// This is a convenience function for `GRanges::map_data().into_array2()`.
    ///
    /// [`Array2`]: ndarray::Array2
    pub fn map_into_array2<F, U>(
        mut self,
        ncols: usize,
        func: F,
    ) -> Result<GRanges<C, Array2<U>>, GRangesError>
    where
        F: FnMut(<T as IntoIterator>::Item) -> Vec<U>,
        T: IntoIterator,
    {
        let total_len = self.len();
        let (ranges, data) = self.take_both()?;

        let flat_vec: Vec<U> = data.into_iter().flat_map(func).collect();
        let transformed_data = Array2::from_shape_vec((total_len, ncols), flat_vec)?;

        assert_eq!(
            transformed_data.shape()[0],
            total_len,
            "Transformed data has different length than GRanges!"
        );

        Ok(GRanges {
            ranges,
            data: Some(transformed_data),
        })
    }

    /// Use [`IntoIterator`] to convert the data container into a new [`Array1`] data container.
    ///
    /// This is like [`GRanges::map_into_array1`] without applying a function to each element.
    ///
    /// [`Array1`]: ndarray::Array1
    pub fn into_array1(
        mut self,
    ) -> Result<GRanges<C, Array1<<T as IntoIterator>::Item>>, GRangesError>
    where
        T: IntoIterator,
    {
        let total_len = self.len();
        let (ranges, data) = self.take_both()?;

        let transformed_data = data
            .into_iter()
            .collect::<Array1<<T as IntoIterator>::Item>>();

        assert_eq!(
            transformed_data.len(),
            total_len,
            "Transformed data has different length than GRanges!"
        );

        Ok(GRanges {
            ranges,
            data: Some(transformed_data),
        })
    }

    /// Use [`IntoIterator`] to convert the data container into a new [`Array2`] data container.
    ///
    /// This is like [`GRanges::map_into_array2`] without applying a function to each element.
    ///
    /// [`Array2`]: ndarray::Array2
    pub fn into_array2(
        mut self,
        ncols: usize,
    ) -> Result<GRanges<C, Array2<<<T as IntoIterator>::Item as IntoIterator>::Item>>, GRangesError>
    where
        C: RangeContainer,
        T: IntoIterator,
        T::Item: IntoIterator,
        <<T as IntoIterator>::Item as IntoIterator>::Item: Copy + 'static,
    {
        let total_len = self.len();
        let (ranges, data) = self.take_both()?;

        let flat_vec: Vec<_> = data
            .into_iter()
            .flat_map(|sub_iter| sub_iter.into_iter())
            .collect();
        let transformed_data = Array2::from_shape_vec((total_len, ncols), flat_vec)?;

        assert_eq!(
            transformed_data.shape()[0],
            total_len,
            "Transformed data has different length than GRanges!"
        );

        Ok(GRanges {
            ranges,
            data: Some(transformed_data),
        })
    }
}

impl<U> IndexedDataContainer for Array1<U>
where
    U: Copy + Default + 'static,
{
    type Item<'a> = U;
    type OwnedItem = U;
    type Output = Array1<U>;

    fn get_value(&self, index: usize) -> Self::Item<'_> {
        self[index]
    }

    fn get_owned(&self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem {
        // already owned
        self[index]
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn is_valid_index(&self, index: usize) -> bool {
        index < self.shape()[0]
    }

    fn new_from_indices(&self, indices: &[usize]) -> Self::Output {
        Array1::from_iter(indices.iter().map(|&idx| self.get_value(idx)))
    }
}

impl<U> IndexedDataContainer for Array2<U>
where
    U: Copy + Default + 'static,
{
    type Item<'a> = ArrayView1<'a, U>;
    type OwnedItem = Array1<U>;
    type Output = Array2<U>;

    fn get_value(&self, index: usize) -> Self::Item<'_> {
        self.row(index)
    }

    fn get_owned(&self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem {
        self.row(index).to_owned()
    }

    fn len(&self) -> usize {
        self.shape()[0]
    }

    fn is_valid_index(&self, index: usize) -> bool {
        index < self.shape()[0]
    }

    fn new_from_indices(&self, indices: &[usize]) -> Self::Output {
        let cols = self.shape()[1];

        let rows_data: Vec<U> = indices
            .iter()
            .flat_map(|&idx| self.row(idx).iter().cloned().collect::<Vec<_>>())
            .collect();

        // create a new Array2<U> from the rows
        // shape is (number of indices, number of columns)
        Array2::from_shape_vec((indices.len(), cols), rows_data)
            .expect("Shape and collected data size mismatch")
    }
}

#[cfg(test)]
mod tests {
    use crate::test_utilities::granges_test_case_01;

    #[test]
    fn test_map_into_array1() {
        let gr = granges_test_case_01();
        let new_gr = gr.map_into_array1(|x| x + 1.0).unwrap();
        // note: 5.0, because 5 ranges each 1 added to
        assert!((new_gr.data().unwrap().sum() - (24.1 + 5.0)).abs() < 1e-4);
    }

    #[test]
    fn test_map_into_array2() {
        let gr = granges_test_case_01();
        let new_gr = gr.map_into_array2(2, |x| vec![x + 1.0, 1.0]).unwrap();
        // note: 5.0, because 5 ranges each 1 added to
        let array2 = new_gr.data().unwrap();
        let first_col_sum = array2.sum_axis(ndarray::Axis(0))[0];
        let second_col_sum = array2.sum_axis(ndarray::Axis(0))[1];
        assert!((first_col_sum - (24.1 + 5.0)).abs() < 1e-4);
        assert!((second_col_sum - (5.0)).abs() < 1e-4);
    }

    #[test]
    fn test_into_array1() {
        let gr = granges_test_case_01();
        let new_gr = gr.into_array1().unwrap();
        assert!((new_gr.data().unwrap().sum() - 24.1).abs() < 1e-4);
    }

    #[test]
    fn test_into_array2() {
        // test case: just wipe out the data and replace with (1, 2)
        // rows; check sum
        let gr = granges_test_case_01();
        let new_gr = gr
            .map_data(|_| vec![1.0, 2.0])
            .unwrap()
            .into_array2(2)
            .unwrap();

        let array2 = new_gr.data().unwrap();
        let first_col_sum = array2.sum_axis(ndarray::Axis(0))[0];
        let second_col_sum = array2.sum_axis(ndarray::Axis(0))[1];
        assert!(((first_col_sum as f64) - 5.0).abs() < 1e-4);
        assert!(((second_col_sum as f64) - 10.0).abs() < 1e-4);
    }
}
