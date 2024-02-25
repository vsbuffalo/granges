//! Data container implementations.
//!

use crate::traits::{DataContainer, IntoDatumType};

pub mod operations;
pub mod vec;

/// These are core supported data types stored in an `enum`, to
/// unify the types that come out of standard operations like
/// `select()`.
#[derive(Debug, Clone)]
pub enum DatumType {
    Float32(f32),
    Float64(f64),
    String(String),
    Integer32(i32),
    Integer64(i64),
    Unsigned32(u32),
    Unsigned64(u64),
    NoValue,
}

impl IntoDatumType for f64 {
    fn into_data_type(self) -> DatumType {
        DatumType::Float64(self)
    }
}
impl IntoDatumType for i64 {
    fn into_data_type(self) -> DatumType {
        DatumType::Integer64(self)
    }
}
impl IntoDatumType for String {
    fn into_data_type(self) -> DatumType {
        DatumType::String(self)
    }
}

impl IntoDatumType for f32 {
    fn into_data_type(self) -> DatumType {
        DatumType::Float32(self)
    }
}

impl IntoDatumType for i32 {
    fn into_data_type(self) -> DatumType {
        DatumType::Integer32(self)
    }
}

impl IntoDatumType for u32 {
    fn into_data_type(self) -> DatumType {
        DatumType::Unsigned32(self)
    }
}

impl IntoDatumType for u64 {
    fn into_data_type(self) -> DatumType {
        DatumType::Unsigned64(self)
    }
}

// Conversion from field types to `DatumType`
impl<T: IntoDatumType> From<T> for DatumType {
    fn from(item: T) -> Self {
        item.into_data_type()
    }
}

impl<U> DataContainer for Vec<U> {}
impl DataContainer for () {}
