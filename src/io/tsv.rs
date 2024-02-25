//! TSV Serializing helpers, functionality, etc.

use crate::{data::DatumType, traits::TsvSerialize};
use lazy_static::lazy_static;

lazy_static! {
    /// The standard BED format TSV configuration.
    pub static ref BED_TSV: TsvConfig = TsvConfig {
        no_value_string: ".".to_string(),
    };
}

/// This is an extensible type to handle common
/// TSV output configurations, e.g. what to print
/// for `None` or [`DatumType::NoValue`].
pub struct TsvConfig {
    pub no_value_string: String,
}

impl TsvSerialize for &String {
    #![allow(unused_variables)]
    fn to_tsv(&self, config: &TsvConfig) -> String {
        self.to_string()
    }
}

impl TsvSerialize for String {
    #![allow(unused_variables)]
    fn to_tsv(&self, config: &TsvConfig) -> String {
        self.to_string()
    }
}

impl TsvSerialize for Option<String> {
    fn to_tsv(&self, config: &TsvConfig) -> String {
        self.as_ref()
            .map_or("".to_string(), |x| format!("\t{}", x.to_tsv(config)))
    }
}

impl<U: TsvSerialize> TsvSerialize for Vec<U> {
    fn to_tsv(&self, config: &TsvConfig) -> String {
        self.iter()
            .map(|x| x.to_tsv(config))
            .collect::<Vec<_>>()
            .join("\t")
    }
}

impl TsvSerialize for &Vec<DatumType> {
    fn to_tsv(&self, config: &TsvConfig) -> String {
        self.iter()
            .map(|x| x.to_tsv(config))
            .collect::<Vec<_>>()
            .join("\t")
    }
}

impl TsvSerialize for DatumType {
    fn to_tsv(&self, config: &TsvConfig) -> String {
        match self {
            DatumType::Float32(val) => val.to_string(),
            DatumType::Float64(val) => val.to_string(),
            DatumType::String(val) => val.clone(),
            DatumType::Integer32(val) => val.to_string(),
            DatumType::Integer64(val) => val.to_string(),
            DatumType::Unsigned32(val) => val.to_string(),
            DatumType::Unsigned64(val) => val.to_string(),
            DatumType::NoValue => config.no_value_string.clone(),
        }
    }
}
