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
pub fn median<F: Float + Sum>(numbers: &mut [F]) -> Option<F> {
    if numbers.is_empty() {
        return None;
    }
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        numbers.select_nth_unstable_by(mid - 1, |a, b| a.partial_cmp(b).unwrap());
        let lower = numbers[mid - 1];
        numbers.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
        let upper = numbers[mid];
        Some((lower + upper) / F::from(2.0).unwrap())
    } else {
        numbers.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
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
    #[inline(always)]
    pub fn run<T: IntoDatumType + Copy>(&self, data: &mut [T]) -> DatumType
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_empty() {
        let mut numbers: Vec<f64> = vec![];
        assert_eq!(median(&mut numbers), None);
    }

    #[test]
    fn test_median_odd() {
        let mut numbers = vec![2.0, 3.0, 1.0];
        assert_eq!(median(&mut numbers), Some(2.0));
    }

    #[test]
    fn test_median_even() {
        let mut numbers = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(median(&mut numbers), Some(2.5));
    }

    #[test]
    fn test_median_negative() {
        let mut numbers = vec![-3.0, -1.0, -2.0];
        assert_eq!(median(&mut numbers), Some(-2.0));
    }
}
