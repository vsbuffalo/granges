//! Benchmarks comparisons against bedtools.
//!
//! For comparison accuracy, any benchmarks conducted here must also have
//! corresponding integration tests in `tests/bedtools_validation.rs`, to
//! ensure the output is the *exact* same.

use criterion::{criterion_group, criterion_main, Criterion};
use granges::test_utilities::{
    granges_binary_path, random_bed3file, random_bed5file, temp_bedfile,
};
use std::{
    fs::File,
    process::{Command, Stdio},
};

const BED_LENGTH: usize = 100_000;

fn bench_range_adjustment(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("adjust");

    // create the test data
    let input_bedfile = random_bed3file(BED_LENGTH);

    // configure the sample size for the group
    group.sample_size(10);

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
    let mut group = c.benchmark_group("filter");

    // create the test data
    let random_bedfile_left_tempfile = random_bed3file(BED_LENGTH);
    let random_bedfile_right_tempfile = random_bed3file(BED_LENGTH);
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

fn bench_flank(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("flank");

    // create the test data
    let random_bedfile_tempfile = random_bed3file(BED_LENGTH);
    let random_bedfile = random_bedfile_tempfile.path();

    // configure the sample size for the group
    // group.sample_size(10);
    group.bench_function("bedtools_flank", |b| {
        b.iter(|| {
            let bedtools_output = Command::new("bedtools")
                .arg("flank")
                .arg("-g")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("-l")
                .arg("20")
                .arg("-r")
                .arg("30")
                .arg("-i")
                .arg(random_bedfile)
                .output()
                .expect("bedtools flank failed");
            assert!(bedtools_output.status.success());
        });
    });

    group.bench_function("granges_filter", |b| {
        b.iter(|| {
            let granges_output = Command::new(granges_binary_path())
                .arg("flank")
                .arg("--genome")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("--left")
                .arg("20")
                .arg("--right")
                .arg("30")
                .arg(random_bedfile)
                .output()
                .expect("granges adjust failed");
            assert!(granges_output.status.success());
        });
    });
}

fn bench_windows(c: &mut Criterion) {
    let width = 1_000_000;
    let step = 10_000;

    // create the benchmark group
    let mut group = c.benchmark_group("windows");

    // configure the sample size for the group
    // group.sample_size(10);
    group.bench_function("bedtools_makewindows", |b| {
        b.iter(|| {
            let bedtools_output = Command::new("bedtools")
                .arg("makewindows")
                .arg("-g")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("-w")
                .arg(width.to_string())
                .arg("-s")
                .arg(step.to_string())
                .output()
                .expect("bedtools makewindows failed");
            assert!(bedtools_output.status.success());
        });
    });

    group.bench_function("granges_windows", |b| {
        b.iter(|| {
            let granges_output = Command::new(granges_binary_path())
                .arg("windows")
                .arg("--genome")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("--width")
                .arg(width.to_string())
                .arg("--step")
                .arg(step.to_string())
                .output()
                .expect("granges windows failed");
            assert!(granges_output.status.success());
        });
    });
}

fn bench_map(c: &mut Criterion) {
    let num_ranges = BED_LENGTH;
    let width = 100_000;
    #[allow(unused_variables)]
    let step = 1_000;

    // make windows
    let windows_file = temp_bedfile();
    let granges_windows_output = Command::new(granges_binary_path())
        .arg("windows")
        .arg("--genome")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--width")
        .arg(width.to_string())
        // .arg("--step")
        // .arg(step.to_string())
        .arg("--output")
        .arg(windows_file.path())
        .output()
        .expect("granges windows failed");
    assert!(
        granges_windows_output.status.success(),
        "{:?}",
        granges_windows_output
    );

    let operations = vec!["sum", "min", "max", "mean", "median"];

    for operation in operations {
        // create the random data BED5
        let bedscores_file = random_bed5file(num_ranges);
        let bedtools_path = temp_bedfile();

        // create the benchmark group
        let mut group = c.benchmark_group(format!("map_{}", operation).to_string());

        // configure the sample size for the group
        group.sample_size(10);
        group.bench_function("bedtools_map", |b| {
            b.iter(|| {
                let bedtools_output_file = File::create(&bedtools_path).unwrap();
                let bedtools_output = Command::new("bedtools")
                    .arg("map")
                    .arg("-a")
                    .arg(windows_file.path())
                    .arg("-b")
                    .arg(&bedscores_file.path())
                    .arg("-c")
                    .arg("5")
                    .arg("-o")
                    .arg(operation)
                    .stdout(Stdio::from(bedtools_output_file))
                    .output()
                    .expect("bedtools map failed");

                assert!(bedtools_output.status.success());
            });
        });

        group.bench_function("granges_map", |b| {
            b.iter(|| {
                let granges_output_file = temp_bedfile();
                let granges_output = Command::new(granges_binary_path())
                    .arg("map")
                    .arg("--genome")
                    .arg("tests_data/hg38_seqlens.tsv")
                    .arg("--left")
                    .arg(windows_file.path())
                    .arg("--right")
                    .arg(bedscores_file.path())
                    .arg("--func")
                    .arg(operation)
                    .arg("--output")
                    .arg(granges_output_file.path())
                    .output()
                    .expect("granges map failed");

                assert!(granges_output.status.success());
            });
        });
    }
}

fn bench_map_all_operations(c: &mut Criterion) {
    let num_ranges = BED_LENGTH;
    let width = 100_000;
    #[allow(unused_variables)]
    let step = 1_000;

    // make windows
    let windows_file = temp_bedfile();
    let granges_windows_output = Command::new(granges_binary_path())
        .arg("windows")
        .arg("--genome")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--width")
        .arg(width.to_string())
        // .arg("--step")
        // .arg(step.to_string())
        .arg("--output")
        .arg(windows_file.path())
        .output()
        .expect("granges windows failed");
    assert!(
        granges_windows_output.status.success(),
        "{:?}",
        granges_windows_output
    );

    // create the random data BED5
    let bedscores_file = random_bed5file(num_ranges);
    let bedtools_path = temp_bedfile();

    // create the benchmark group
    let mut group = c.benchmark_group("map_multiple");

    // configure the sample size for the group
    group.sample_size(10);
    group.bench_function("bedtools_map", |b| {
        b.iter(|| {
            let bedtools_output_file = File::create(&bedtools_path).unwrap();
            let bedtools_output = Command::new("bedtools")
                .arg("map")
                .arg("-a")
                .arg(windows_file.path())
                .arg("-b")
                .arg(&bedscores_file.path())
                .arg("-c")
                .arg("5")
                .arg("-o")
                .arg("min,max,mean,sum,median")
                .stdout(Stdio::from(bedtools_output_file))
                .output()
                .expect("bedtools map failed");

            assert!(bedtools_output.status.success());
        });
    });

    group.bench_function("granges_map", |b| {
        b.iter(|| {
            let granges_output_file = temp_bedfile();
            let granges_output = Command::new(granges_binary_path())
                .arg("map")
                .arg("--genome")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("--left")
                .arg(windows_file.path())
                .arg("--right")
                .arg(bedscores_file.path())
                .arg("--func")
                .arg("min,max,mean,sum,median")
                .arg("--output")
                .arg(granges_output_file.path())
                .output()
                .expect("granges map failed");

            assert!(granges_output.status.success());
        });
    });
}

criterion_group!(
    benches,
    bench_filter_adjustment,
    bench_range_adjustment,
    bench_flank,
    bench_windows,
    bench_map,
    bench_map_all_operations,
);
criterion_main!(benches);
