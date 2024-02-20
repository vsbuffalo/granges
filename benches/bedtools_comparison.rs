//! Benchmarks comparisons against bedtools.
//!
//! For comparison accuracy, any benchmarks conducted here must also have
//! corresponding integration tests in `tests/bedtools_validation.rs`, to
//! ensure the output is the *exact* same.

use criterion::{criterion_group, criterion_main, Criterion};
use granges::test_utilities::{granges_binary_path, random_bedfile};
use std::process::Command;

const BED_LENGTH: usize = 10_000;

fn bench_range_adjustment(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("slop vs adjust");

    // create the test data
    let input_bedfile = random_bedfile(BED_LENGTH);

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
                .arg("--genome")
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

fn bench_filter_adjustment(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("intersect vs filter");

    // create the test data
    let random_bedfile_left_tempfile = random_bedfile(BED_LENGTH);
    let random_bedfile_right_tempfile = random_bedfile(BED_LENGTH);
    let random_bedfile_left = random_bedfile_left_tempfile.path();
    let random_bedfile_right = random_bedfile_right_tempfile.path();

    // configure the sample size for the group
    // group.sample_size(10);
    group.bench_function("bedtools_intersect", |b| {
        b.iter(|| {
            let bedtools_output = Command::new("bedtools")
                .arg("intersect")
                .arg("-a")
                .arg(&random_bedfile_left)
                .arg("-b")
                .arg(&random_bedfile_right)
                .arg("-wa")
                .arg("-u")
                .output()
                .expect("bedtools intersect failed");
            assert!(bedtools_output.status.success());
        });
    });

    group.bench_function("granges_filter", |b| {
        b.iter(|| {
            let granges_output = Command::new(granges_binary_path())
                .arg("filter")
                .arg("--genome")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("--left")
                .arg(&random_bedfile_left)
                .arg("--right")
                .arg(&random_bedfile_right)
                .output()
                .expect("granges adjust failed");
            assert!(granges_output.status.success());
        });
    });
}

criterion_group!(benches, bench_filter_adjustment, bench_range_adjustment);
criterion_main!(benches);
