//! Benchmarks comparisons against bedtools.
//!
//! For comparison accuracy, any benchmarks conducted here must also have
//! corresponding integration tests in `tests/bedtools_validation.rs`, to
//! ensure the output is the *exact* same.

use criterion::{criterion_group, criterion_main, Criterion};
use granges::{commands::granges_random_bed, test_utilities::granges_binary_path};
use std::process::Command;
use tempfile::NamedTempFile;

const BED_LENGTH: u32 = 10_000;

/// Create a random BED3 file based on the hg38 sequence lengths, and write to disk.
pub fn random_bedfile() -> NamedTempFile {
    let temp_file = NamedTempFile::new().expect("Failed to create temp file");
    let random_bedfile_path = temp_file.path().to_path_buf();
    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        BED_LENGTH,
        Some(&random_bedfile_path),
        true,
    )
    .expect("could not generate random BED file");
    temp_file
}

fn bench_range_adjustment(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("slop vs adjust");

    // create the test data
    let input_bedfile = random_bedfile();

    // configure the sample size for the group
    // group.sample_size(10);

    // bedtools slop
    group.bench_function("bedtools_slop", |b| {
        b.iter(|| {
            let bedtools_output = Command::new("bedtools")
                .arg("slop")
                .arg("-g")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("-b")
                .arg("10")
                .arg("-i")
                .arg(input_bedfile.path())
                .output()
                .expect("bedtools slop failed");
            assert!(bedtools_output.status.success());
        });
    });

    group.bench_function("granges_adjust", |b| {
        b.iter(|| {
            let granges_output = Command::new(granges_binary_path())
                .arg("adjust")
                .arg("--seqlens")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("--both")
                .arg("10")
                .arg("--sort")
                .arg(input_bedfile.path())
                .output()
                .expect("granges adjust failed");
            assert!(granges_output.status.success());
        });
    });
}
criterion_group!(benches, bench_range_adjustment);
criterion_main!(benches);
