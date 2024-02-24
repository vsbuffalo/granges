//! Data container implementations.
//!

use crate::traits::DataContainer;

pub mod operations;
pub mod vec;

impl<U> DataContainer for Vec<U> {}
impl DataContainer for () {}
