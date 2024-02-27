//! Implementations of various operations on data.
//!
//! # 🐌 Performance Note
//!
//! Currently operations are not much faster than bedtools (11%-13%).
//! These methods can be made faster by looping over data once, collecting
//! the quantities that may make up different statistics.

use clap::ValueEnum;
use num_traits::{Float, ToPrimitive};
use std::iter::Sum;
use ndarray;

use super::DatumType;
use crate::traits::IntoDatumType;

/// The (subset of) standard `bedtools map` operations.
#[derive(Clone, Debug, ValueEnum)]
pub enum Operation {
    Sum,
    Min,
    Max,
    Mean,
    Median,
    Collapse,
}

impl Operation {
    pub fn requires_sort(&self) -> bool {
        match self {
            Operation::Median => true,
            Operation::Max => true,
            Operation::Min => true,
            _ => false,
        }
    }

    /// Run an operation on the supplied [`Array1`] of data.
    ///
    /// Note that order statistics ([`Operation::Min`], [`Operation::Max`], 
    /// [`Operation::Median`]) require the input data be **presorted**. This
    /// sorting is *not* validated here.
    #[inline(always)]
    pub fn run<T: IntoDatumType + Copy>(&self, data: &ndarray::Array1<T>) -> DatumType
    where
        T: Float + Sum<T> + ToPrimitive + Clone + ToString,
    {
        match self {
            Operation::Sum => {
                let sum = data.sum();
                sum.into_data_type()
            }
            Operation::Min => {
                if data.is_empty() {
                    return DatumType::NoValue
                }
                let min = data[0];
                min.into_data_type()
            }
            Operation::Max => {
                if data.is_empty() {
                    return DatumType::NoValue
                }
                let max = data[data.len()-1];
                max.into_data_type()
            }
            Operation::Mean => {
                if data.is_empty() {
                    DatumType::NoValue
                } else {
                    let sum = data.sum();
                    let mean = sum / T::from(data.len()).unwrap();
                    mean.into_data_type()
                }
            }
            Operation::Median => {
                let mid = data.len() / 2;
                if data.len() % 2 == 0 {
                    let lower = data[mid - 1];
                    let upper = data[mid];
                    let median = (lower + upper) / T::from(2.0).unwrap();
                    median.into_data_type()
                } else {
                    let median = data[mid];
                    median.into_data_type()
                }
            },
            Operation::Collapse => {
                let collapsed = data
                    .iter()
                    .map(|num| num.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                DatumType::String(collapsed)
            }
        }
    }
}

#[cfg(test)]
mod tests {
}
