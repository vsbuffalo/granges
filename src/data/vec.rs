//! Data container implementations for [`Vec<U>`].

use crate::traits::IndexedDataContainer;

/// Trait methods for the commonly-used `Vec<U>` data container.
///
/// Note that the associated `Item` type is always a *reference* to the data elements.
impl<'a, U> IndexedDataContainer<'a> for Vec<U>
where
    U: Clone + 'a,
{
    type Item = &'a U;
    type Output = Vec<U>;

    fn get_value(&'a self, index: usize) -> Self::Item {
        self.get(index).unwrap()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn is_valid_index(&self, index: usize) -> bool {
        self.get(index).is_some()
    }

    fn new_from_indices(&self, indices: &[usize]) -> Self::Output {
        Vec::from_iter(indices.iter().map(|&idx| (*self.get_value(idx)).clone()))
    }
}
