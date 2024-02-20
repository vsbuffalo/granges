//! Data container implementations.

use crate::traits::DataContainer;

pub mod vec;

impl<U> DataContainer for Vec<U> {}
impl DataContainer for () {}
