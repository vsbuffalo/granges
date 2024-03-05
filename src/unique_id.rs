//! The [`UniqueIdentifier`] type, which stores a mapping between
//! some key type and a `usize` index.
//!
//! Note that these are not thread-safe currently but
//! it is relatively easy to do -- submit a GitHub
//! issue if you need this feature.

use std::{
    cell::RefCell,
    collections::{hash_map::Keys, HashMap},
    hash::Hash,
};

use crate::traits::{DataContainer, IndexedDataContainer};

impl<K: std::cmp::Eq + Hash + Clone> DataContainer for UniqueIdentifier<K> {}
impl<K: std::cmp::Eq + Hash + Clone> IndexedDataContainer for UniqueIdentifier<K> {
    type Item<'a> = K where K: 'a;
    type OwnedItem = K;

    fn get_value(&self, index: usize) -> <Self as IndexedDataContainer>::Item<'_> {
        self.get_key(index).unwrap().clone()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn get_owned(&self, index: usize) -> <Self as IndexedDataContainer>::OwnedItem {
        self.get_value(index).clone()
    }

    fn is_valid_index(&self, index: usize) -> bool {
        self.ensure_reverse_map_created();
        self.reverse_map
            .borrow()
            .as_ref()
            .unwrap()
            .contains_key(&index)
    }
}

#[derive(Clone, Debug)]
pub struct UniqueIdentifier<K> {
    map: HashMap<K, usize>,
    next_key: usize,
    reverse_map: RefCell<Option<HashMap<usize, K>>>,
}

impl<K> Default for UniqueIdentifier<K>
where
    K: std::cmp::Eq + Hash + Clone,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<K> UniqueIdentifier<K>
where
    K: std::cmp::Eq + Hash + Clone,
{
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            next_key: 0,
            reverse_map: RefCell::new(None),
        }
    }

    pub fn get_or_insert(&mut self, s: &K) -> usize {
        match self.map.get(s) {
            Some(&key) => key,
            None => {
                let key = self.next_key;
                self.map.insert(s.clone(), key);
                self.next_key += 1;
                // we need to invalidate the cached reverse map
                let mut reverse_map_borrow = self.reverse_map.borrow_mut();
                *reverse_map_borrow = None;
                key
            }
        }
    }

    fn ensure_reverse_map_created(&self) {
        let mut reverse_map_borrow = self.reverse_map.borrow_mut();
        if reverse_map_borrow.is_none() {
            let reverse_map = self.map.iter().map(|(k, v)| (*v, k.clone())).collect();
            *reverse_map_borrow = Some(reverse_map);
        }
    }

    pub fn get_key(&self, index: usize) -> Option<K> {
        self.ensure_reverse_map_created();
        self.reverse_map
            .borrow()
            .as_ref()
            .unwrap()
            .get(&index)
            .cloned()
    }

    pub fn is_valid_index(&self, index: usize) -> bool {
        self.ensure_reverse_map_created();
        self.reverse_map
            .borrow()
            .as_ref()
            .unwrap()
            .contains_key(&index)
    }

    pub fn keys(&self) -> Keys<'_, K, usize> {
        self.map.keys()
    }

    pub fn indices(&self) -> Vec<usize> {
        self.ensure_reverse_map_created();
        self.reverse_map
            .borrow()
            .as_ref()
            .unwrap()
            .keys()
            .cloned()
            .collect()
    }

    pub fn get_index(&self, s: &K) -> Option<usize> {
        self.map.get(s).cloned()
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert_and_get_key() {
        let mut ui = UniqueIdentifier::new();
        let key = "test_key";
        let index = ui.get_or_insert(&key);

        assert_eq!(index, 0);
        assert_eq!(ui.get_key(index), Some(key));
    }

    #[test]
    fn test_is_valid_index() {
        let mut ui = UniqueIdentifier::new();
        let key = "another_key";
        ui.get_or_insert(&key);

        assert!(ui.is_valid_index(0));
        assert!(!ui.is_valid_index(1));

        let key = "another!";
        let index = ui.get_or_insert(&key);
        assert_eq!(index, 1);
        assert!(ui.is_valid_index(index));
    }

    #[test]
    fn test_len() {
        let mut ui = UniqueIdentifier::new();
        assert_eq!(ui.len(), 0);

        ui.get_or_insert(&"key1");
        assert_eq!(ui.len(), 1);

        ui.get_or_insert(&"key2");
        assert_eq!(ui.len(), 2);

        // Inserting an existing key should not increase the length
        ui.get_or_insert(&"key1");
        assert_eq!(ui.len(), 2);
    }

    #[test]
    fn test_indices() {
        let mut ui = UniqueIdentifier::new();
        ui.get_or_insert(&"key1");
        ui.get_or_insert(&"key2");

        let indices = ui.indices();
        assert_eq!(indices.len(), 2);
        assert!(indices.contains(&0));
        assert!(indices.contains(&1));
    }

    #[test]
    fn test_get_index() {
        let mut ui = UniqueIdentifier::new();
        let key = "key_to_find";
        ui.get_or_insert(&key);

        let index = ui.get_index(&key);
        assert_eq!(index, Some(0));

        let missing_key = "missing_key";
        assert_eq!(ui.get_index(&missing_key), None);
    }
}
