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

use super::DatumType;
use crate::traits::IntoDatumType;

/// Calculate the median.
///
/// This will clone and turn `numbers` into a `Vec`.
pub fn median<F: Float + Sum>(numbers: &[F]) -> Option<F> {
    if numbers.is_empty() {
        return None;
    }
    let mut numbers = numbers.to_vec();
    numbers.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        Some((numbers[mid - 1] + numbers[mid]) / F::from(2.0).unwrap())
    } else {
        Some(numbers[mid])
    }
}

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
    pub fn run<T: IntoDatumType + Copy>(&self, data: &[T]) -> DatumType
    where
        T: Float + Sum<T> + ToPrimitive + Clone + ToString,
    {
        match self {
            Operation::Sum => {
                let sum: T = data.iter().copied().sum();
                sum.into_data_type()
            }
            Operation::Min => {
                let min = data
                    .iter()
                    .filter(|x| x.is_finite())
                    .copied()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));
                min.map_or(DatumType::NoValue, |x| x.into_data_type())
            }
            Operation::Max => {
                let max = data
                    .iter()
                    .filter(|x| x.is_finite())
                    .copied()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Less));
                max.map_or(DatumType::NoValue, |x| x.into_data_type())
            }
            Operation::Mean => {
                if data.is_empty() {
                    DatumType::NoValue
                } else {
                    let sum: T = data.iter().copied().sum();
                    let mean = sum / T::from(data.len()).unwrap();
                    mean.into_data_type()
                }
            }
            Operation::Median => median(data).map_or(DatumType::NoValue, |x| x.into_data_type()),
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
