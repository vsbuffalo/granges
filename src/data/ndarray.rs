//! Data container implementations for [`ndarray::Array1`] and [`ndarray::Array2`].

use ndarray::{Array1, Array2, ArrayView1};
use crate::traits::IndexedDataContainer;

impl<'a, U> IndexedDataContainer<'a> for Array1<U>
where
    U: Copy + Default + 'a,
{
    type Item = U;
    type OwnedItem = U;
    type Output = Array1<U>;

    fn get_value(&'a self, index: usize) -> Self::Item {
        self[index]
    }

    fn get_owned(&'a self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem {
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

impl<'a, U> IndexedDataContainer<'a> for Array2<U>
where
    U: Copy + Default + 'a,
{
    type Item = ArrayView1<'a, U>;
    type OwnedItem = Array1<U>;
    type Output = Array2<U>;

    fn get_value(&'a self, index: usize) -> Self::Item {
        self.row(index)
    }

    fn get_owned(&'a self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem {
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
