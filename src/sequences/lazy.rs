//! Functionality for lazy-loading sequences off disk into memory.
//!
//! The main functionality is the very generic [`LazyLoader`]. This is generic over the loading
//! function and the key type. It mainly handles loading data into a [`RefCell`].
//!
use std::cell::{Ref, RefCell};

use crate::prelude::GRangesError;

/// A lazy-loader function that takes a reader type `R` and
/// uses it to load in data of type `T`.
type LoaderFunc<R, T, K> = Box<dyn Fn(&mut R, &K) -> Result<T, GRangesError>>;

/// Lazy loader, which uses [`RefCell`] to store mutable reader and data, used for lazy loading and
/// storing one key's worth of data data.
///
/// Generic key types allow for keys to be region tuples, etc.
///
/// # Generics
///  * `R`: the reader type.
///  * `T`: the data type.
///  * `K`: the key type.
///
///
///  * `key`: they key of the loaded data (`None` if nothing loaded).
///  * `reader`: the open reader type (e.g. an indexed FASTA file).
///  * `loader`: a function that takes a key type `K` and an open reader,
///              and retrieves the right data of type `T`.
///  * `data`: the data (`None` if no data loaded).
pub struct LazyLoader<R, T, K>
where
    K: std::fmt::Debug,
    T: std::fmt::Debug,
{
    key: RefCell<Option<K>>,
    reader: RefCell<R>,
    loader: LoaderFunc<R, T, K>,
    data: RefCell<Option<T>>,
}

impl<R, T, K: std::fmt::Debug> std::fmt::Debug for LazyLoader<R, T, K>
where
    K: std::fmt::Debug,
    T: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LazyLoader")
            .field("key:", &self.key)
            .field("data", &self.data)
            .finish_non_exhaustive()
    }
}

impl<R, T, K> LazyLoader<R, T, K>
where
    K: std::fmt::Debug,
    T: std::fmt::Debug,
    K: Clone + PartialEq,
{
    pub fn new<F>(reader: R, loader: F) -> LazyLoader<R, T, K>
    where
        F: Fn(&mut R, &K) -> Result<T, GRangesError> + 'static,
    {
        LazyLoader {
            key: RefCell::new(None),
            reader: RefCell::new(reader),
            loader: Box::new(loader),
            data: RefCell::new(None),
        }
    }

    /// Return a `bool` indicating whether the specified `key` is cached.
    ///
    /// # Arguments
    /// * `key`: a `&str` key passed to the loader function.
    pub fn is_loaded(&self, key: &K) -> bool {
        let loaded_key = self.key.borrow();
        match &*loaded_key {
            None => false,
            Some(existing_key) => *key == *existing_key,
        }
    }

    fn load(&self, key: &K) -> Result<T, GRangesError> {
        let mut reader = self.reader.borrow_mut();
        let new_data = (self.loader)(&mut reader, key)?;
        Ok(new_data)
    }

    /// Clear out the cache.
    pub fn clear(&self) {
        let mut data_borrow = self.data.borrow_mut();
        *data_borrow = None;
        *self.key.borrow_mut() = None;
    }

    /// Check if the cache is empty.
    pub fn is_empty(&self) -> bool {
        let data_borrow = self.data.borrow_mut();
        let is_empty = data_borrow.is_none();
        if is_empty {
            let key_borrow = self.key.borrow_mut();
            assert_eq!(*key_borrow, None, "invalid state: empty cache, but no key");
        }
        is_empty
    }

    /// Load the data corresponding to `key`
    ///
    /// Loads the data associated with a key using the
    /// `loader` function.
    ///
    /// # Arguments
    /// * `key`: a `&str` key passed to the loader function.
    ///
    /// # Returns
    /// Returns a [`Result<T, GRangesError>`].
    pub fn get_data(&self, key: &K) -> Result<Ref<T>, GRangesError> {
        {
            let mut data_borrow = self.data.borrow_mut();
            let is_loaded = self.is_loaded(key);
            if !is_loaded {
                let new_data = self.load(key)?;
                *data_borrow = Some(new_data);
                *self.key.borrow_mut() = Some(key.clone());
            }
        }

        let data_ref = self.data.borrow();
        Ok(Ref::map(data_ref, |opt| {
            opt.as_ref().unwrap_or_else(|| {
                panic!("Data should be loaded at this point, but it's not available.")
            })
        }))
    }
}

pub enum LoadedType {
    Sequence(String),
    Region(String, i32, i32),
}

// pub struct LazyIndexedLoader<R, T> {
//     key: RefCell<LoadedType>,
//     reader: RefCell<R>,
//     loader: LoaderFunc<R, T>,
//     query: LoaderFunc<R, T>,
//     data: RefCell<Option<T>>,
// }
//
// impl<R, T> LazyIndexedLoader<R, T> {
//     pub fn new<F, P>(reader: R, loader: F, query: P) -> LazyLoader<R, T>
//     where
//         F: Fn(&mut R, &str) -> Result<T, GRangesError> + 'static,
//         P: Fn(&mut R, &str, i32, i32) -> Result<T, GRangesError> + 'static,
//     {
//         LazyLoader {
//             key: RefCell::new(String::new()),
//             reader: RefCell::new(reader),
//             loader: Box::new(loader),
//             data: RefCell::new(None),
//         }
//     }
//
//     /// Return a `bool` indicating whether the specified `key` is cached.
//     ///
//     /// # Arguments
//     /// * `key`: a `&str` key passed to the loader function.
//     pub fn is_loaded(&self, key: &str) -> bool {
//         let loaded_key = self.key.borrow();
//         *loaded_key == key
//     }
//
//     fn load(&self, key: &str) -> Result<T, GRangesError> {
//         let mut reader = self.reader.borrow_mut();
//         let new_data = (self.loader)(&mut reader, key)?;
//         Ok(new_data)
//     }
//
//     /// Clear out the cache.
//     pub fn clear(&self) {
//         let mut data_borrow = self.data.borrow_mut();
//         *data_borrow = None;
//         *self.key.borrow_mut() = String::new();
//     }
//
//     /// Check if the cache is empty.
//     pub fn is_empty(&self) -> bool {
//         let data_borrow = self.data.borrow_mut();
//         let is_empty = data_borrow.is_none();
//         if is_empty {
//             let key_borrow = self.key.borrow_mut();
//             assert_eq!(
//                 *key_borrow,
//                 String::new(),
//                 "invalid state: empty cache, but no key"
//             );
//         }
//         is_empty
//     }
//
//     /// Load the data corresponding to `key`
//     ///
//     /// Loads the data associated with a key using the
//     /// `loader` function.
//     ///
//     /// # Arguments
//     /// * `key`: a `&str` key passed to the loader function.
//     ///
//     /// # Returns
//     /// Returns a `Result<T, GRangesError>`.
//     pub fn get_data(&self, key: &str) -> Result<Ref<T>, GRangesError> {
//         {
//             let mut data_borrow = self.data.borrow_mut();
//             let is_loaded = self.is_loaded(key);
//             if !is_loaded {
//                 let new_data = self.load(key)?;
//                 *data_borrow = Some(new_data);
//                 *self.key.borrow_mut() = key.to_string();
//             }
//         }
//
//         let data_ref = self.data.borrow();
//         Ok(Ref::map(data_ref, |opt| {
//             opt.as_ref().unwrap_or_else(|| {
//                 panic!("Data should be loaded at this point, but it's not available.")
//             })
//         }))
//     }
// }
