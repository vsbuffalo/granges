//! Data container implementations.
//!

use crate::{
    io::TsvConfig,
    traits::{DataContainer, IntoDatumType},
};
use serde::ser::Serializer;
use serde::Serialize;

#[cfg(feature = "ndarray")]
pub mod ndarray;
pub mod operations;
pub mod vec;

impl<U> DataContainer for Vec<U> {}
impl DataContainer for () {}

/// These are core supported data types stored in an `enum`, to
/// unify the types that come out of standard operations of
/// heterogeneous output types.
#[derive(Debug, Clone, Serialize)]
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

impl DatumType {
    pub fn into_serializable(self, config: &TsvConfig) -> SerializableDatumType {
        SerializableDatumType {
            datum: self,
            config,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SerializableDatumType<'a> {
    pub datum: DatumType,
    pub config: &'a TsvConfig,
}

impl<'a> Serialize for SerializableDatumType<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match &self.datum {
            DatumType::NoValue => serializer.serialize_str(&self.config.no_value_string),
            DatumType::Float32(value) => serializer.serialize_str(&value.to_string()),
            DatumType::Float64(value) => serializer.serialize_str(&value.to_string()),
            DatumType::String(value) => serializer.serialize_str(value),
            DatumType::Integer32(value) => serializer.serialize_str(&value.to_string()),
            DatumType::Integer64(value) => serializer.serialize_str(&value.to_string()),
            DatumType::Unsigned32(value) => serializer.serialize_str(&value.to_string()),
            DatumType::Unsigned64(value) => serializer.serialize_str(&value.to_string()),
        }
    }
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
