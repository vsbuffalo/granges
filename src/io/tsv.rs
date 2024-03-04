//! TSV Serializing helpers, functionality, etc.

use lazy_static::lazy_static;

lazy_static! {
    /// The standard BED format TSV configuration.
    pub static ref BED_TSV: TsvConfig = TsvConfig {
        no_value_string: ".".to_string(),
        headers: None,
    };
}

/// This is an extensible type to handle common
/// TSV output configurations, e.g. what to print
/// for `None` or [`DatumType::NoValue`].
#[derive(Debug, Clone)]
pub struct TsvConfig {
    pub no_value_string: String,
    pub headers: Option<Vec<String>>,
}
