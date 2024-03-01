//! Validation against bedtools

use granges::{
    commands::granges_random_bed,
    io::{parsers::bed::bed_missing, InputStream},
    prelude::{read_seqlens, BedlikeIterator, GRanges, GenomicRangesFile, TsvRecordIterator},
    ranges::GenomicRangeRecord,
    test_utilities::{granges_binary_path, random_bed3file, random_bed5file, temp_bedfile},
};
use std::{
    fs::File,
    io::BufRead,
    path::PathBuf,
    process::{Command, Stdio},
};

use serde::Deserialize;

// adjustable test sizes -- useful also for making very small
// so string diffs are easier to compare.
#[cfg(not(feature = "bench-big"))]
const BED_LENGTH: usize = 10;
#[cfg(feature = "bench-big")]
const BED_LENGTH: usize = 1_000_000;

// helpers

fn head_file(filepath: impl Into<PathBuf>) {
    let filepath = filepath.into();
    let file = InputStream::new(&filepath);
    let mut reader = file.reader().unwrap();
    for _ in 0..10 {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        dbg!(&line);
    }
}

#[macro_export]
macro_rules! assert_float_tol {
    ($a:expr, $b:expr) => {{
        let tol = 1e-5;
        assert!(
            ($a - $b).abs() < tol,
            "assertion failed: `(left !== right)`\n left: `{:?}`, right: `{:?}`, tolerance: `{:?}`",
            $a,
            $b,
            tol
        );
    }};
}

#[macro_export]
macro_rules! assert_option_float_tol {
    ($a:expr, $b:expr) => {{
        let tol = 1e-5;
        match ($a, $b) {
            (Some(a_val), Some(b_val)) => {
                assert!(
                    (a_val - b_val).abs() < tol,
                    "assertion failed: `(left !== right)`\n left: `{:?}`, right: `{:?}`, tolerance: `{:?}`",
                    a_val,
                    b_val,
                    tol
                    );
            }
            (None, None) => (), // Both are None, so do nothing.
            _ => panic!(
                "assertion failed: one is None and the other is Some\n left: `{:?}`, right: `{:?}`",
                $a, $b
                ),
        }
    }};
}

#[test]
fn test_random_bed3file_filetype_detect() {
    let random_bedfile_path = temp_bedfile().path().to_path_buf();

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        BED_LENGTH,
        Some(&random_bedfile_path),
        true,
        false,
    )
    .expect("could not generate random BED file");

    head_file(&random_bedfile_path);
    match GenomicRangesFile::detect(random_bedfile_path).unwrap() {
        GenomicRangesFile::Bed3(_) => (),
        _ => panic!("<IN-TESTS>: could not detect correct filetype"),
    }
}

#[test]
fn test_against_bedtools_slop() {
    let random_bedfile = temp_bedfile();
    let random_bedfile_path = random_bedfile.path();

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        BED_LENGTH,
        Some(&random_bedfile_path),
        true,
        false,
    )
    .expect("could not generate random BED file");

    let width = 10;

    let bedtools_output = Command::new("bedtools")
        .arg("slop")
        .arg("-g")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("-b")
        .arg(width.to_string())
        .arg("-i")
        .arg(&random_bedfile_path)
        .output()
        .expect("bedtools slop failed");

    let granges_output = Command::new(granges_binary_path())
        .arg("adjust")
        .arg("--genome")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--both")
        .arg(width.to_string())
        .arg("--sort")
        .arg(&random_bedfile_path)
        .output()
        .expect("granges adjust failed");

    assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
    assert!(granges_output.status.success(), "{:?}", granges_output);

    assert_eq!(
        String::from_utf8_lossy(&bedtools_output.stdout),
        String::from_utf8_lossy(&granges_output.stdout)
    );
}

/// Test bedtools intersect -a <left> -b <right> -wa -u
/// against
/// granges filter --genome <genome> --left <left> --right <right>
#[test]
fn test_against_bedtools_intersect_wa() {
    let num_ranges = 1_000_000;

    let random_bedfile_left_tempfile = random_bed3file(num_ranges);
    let random_bedfile_right_tempfile = random_bed3file(num_ranges);
    let random_bedfile_left = random_bedfile_left_tempfile.path();
    let random_bedfile_right = random_bedfile_right_tempfile.path();

    // for testing: uncomment and results are local for inspection
    // let random_bedfile_left = Path::new("test_left.bed");
    // let random_bedfile_right = Path::new("test_right.bed");

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        num_ranges,
        Some(&random_bedfile_right),
        true,
        false,
    )
    .expect("could not generate random BED file");

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

    assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
    assert!(granges_output.status.success(), "{:?}", granges_output);

    assert_eq!(
        String::from_utf8_lossy(&bedtools_output.stdout),
        String::from_utf8_lossy(&granges_output.stdout)
    );
}

/// Test bedtools flank -g <genome> -i <input> -l 10 -r 20
/// against
/// granges filter --genome <genome> --left 10 --right 20 <input>
#[test]
fn test_against_bedtools_flank() {
    let num_ranges = 1_000;

    let random_bedfile_tempfile = random_bed3file(num_ranges);
    let random_bedfile = random_bedfile_tempfile.path();

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        num_ranges,
        Some(&random_bedfile),
        true,
        false,
    )
    .expect("could not generate random BED file");

    let bedtools_output = Command::new("bedtools")
        .arg("flank")
        .arg("-g")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("-l")
        .arg("20")
        .arg("-r")
        .arg("30")
        .arg("-i")
        .arg(&random_bedfile)
        .output()
        .expect("bedtools flank failed");

    let granges_output = Command::new(granges_binary_path())
        .arg("flank")
        .arg("--genome")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--left")
        .arg("20")
        .arg("--right")
        .arg("30")
        .arg(&random_bedfile)
        .output()
        .expect("granges flank failed");

    assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
    assert!(granges_output.status.success(), "{:?}", granges_output);

    let bedtools_str = String::from_utf8_lossy(&bedtools_output.stdout);
    let granges_str = String::from_utf8_lossy(&granges_output.stdout);

    let mut bedtools_ranges: Vec<_> = bedtools_str.split("\n").collect();
    let mut granges_ranges: Vec<_> = granges_str.split("\n").collect();

    bedtools_ranges.sort();
    granges_ranges.sort();

    assert_eq!(bedtools_ranges, granges_ranges);
}

#[test]
fn test_against_bedtools_makewindows() {
    // some weird widths, steps to try to catch remainder issues
    let widths = vec![131123, 1_000_0013];
    let steps = vec![10_001, 10_113];

    for width in widths.iter() {
        for step in steps.iter() {
            let bedtools_output = Command::new("bedtools")
                .arg("makewindows")
                .arg("-g")
                .arg("tests_data/hg38_seqlens.tsv")
                .arg("-w")
                .arg(width.to_string())
                .arg("-s")
                .arg(step.to_string())
                .output()
                .expect("bedtools slop failed");

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

            assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
            assert!(granges_output.status.success(), "{:?}", granges_output);

            assert_eq!(
                String::from_utf8_lossy(&bedtools_output.stdout),
                String::from_utf8_lossy(&granges_output.stdout)
            );
        }
    }
}

#[test]
fn test_against_bedtools_map() {
    let num_ranges = BED_LENGTH;
    let width = 1_000_000;
    #[allow(unused_variables)]
    let step = 10_000; // can uncomment lines below to test this

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

    // we're going to test all of these operations
    // TODO/TEST need to test collapse
    let operations = vec!["sum", "min", "max", "mean", "median"];

    for operation in operations {
        // create the random data BED5
        let bedscores_file = random_bed5file(num_ranges);

        let bedtools_path = temp_bedfile();
        let bedtools_output_file = File::create(&bedtools_path).unwrap();

        // compare map commands
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

        assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
        assert!(granges_output.status.success(), "{:?}", granges_output);

        let genome = read_seqlens("tests_data/hg38_seqlens.tsv").unwrap();

        let bedtools_iter = BedlikeIterator::new(bedtools_path.path()).unwrap();
        let mut bedtools_gr = GRanges::from_iter(bedtools_iter, &genome).unwrap();

        let granges_iter = BedlikeIterator::new(granges_output_file.path().to_path_buf()).unwrap();
        let mut granges_gr = GRanges::from_iter(granges_iter, &genome).unwrap();

        let granges_data = granges_gr.take_data().unwrap();
        let granges_data = granges_data.iter().map(|extra_cols| {
            // explicitly handle missing, then unwrap (rather than ok(), which is less safe)
            if *extra_cols == Some(".".to_string()) {
                return None;
            }
            let score: Result<f64, _> = extra_cols.as_ref().unwrap().parse();
            // dbg!(&extra_cols);
            Some(score.unwrap())
        });

        let bedtools_data = bedtools_gr.take_data().unwrap();
        let bedtools_data = bedtools_data.iter().map(|extra_cols| {
            if *extra_cols == Some(".".to_string()) {
                return None;
            }
            let score: f64 = extra_cols.as_ref().unwrap().parse().unwrap();
            Some(score)
        });
        assert_eq!(granges_data.len(), bedtools_data.len());

        granges_data
            .zip(bedtools_data)
            .for_each(|(gr_val, bd_val)| match (gr_val, bd_val) {
                (Some(gr), Some(bd)) => assert!((gr - bd).abs() < 1e-5),
                // NOTE: for some sum operations with no data,
                // bedtools returns '.' not 0.0. The latter is more correct
                // (the sum of the empty set is not NA, it's 0.0).
                // This is a shim so tests don't stochastically break
                // in this case.
                (Some(n), None) if n == 0.0 => (),
                (None, None) => (),
                _ => panic!("{:?}", (&operation, &gr_val, &bd_val)),
            });
    }
}

#[test]
fn test_against_bedtools_map_multiple() {
    // try multiple operations at once
    let num_ranges = BED_LENGTH;
    let width = 1_000_000;
    #[allow(unused_variables)]
    let step = 10_000; // can uncomment lines below to test this

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
    let bedtools_output_file = File::create(&bedtools_path).unwrap();

    // compare map commands
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
        .arg("min,max,mean,sum-not-empty,median")
        .arg("--output")
        .arg(granges_output_file.path())
        .output()
        .expect("granges map failed");

    assert!(bedtools_output.status.success(), "{:?}", bedtools_output);
    assert!(granges_output.status.success(), "{:?}", granges_output);

    let genome = read_seqlens("tests_data/hg38_seqlens.tsv").unwrap();

    #[derive(Deserialize)]
    struct Stats {
        #[serde(deserialize_with = "bed_missing")]
        min: Option<f64>,
        #[serde(deserialize_with = "bed_missing")]
        max: Option<f64>,
        #[serde(deserialize_with = "bed_missing")]
        mean: Option<f64>,
        #[serde(deserialize_with = "bed_missing")]
        sum: Option<f64>,
        #[serde(deserialize_with = "bed_missing")]
        median: Option<f64>,
    }

    head_file(&bedtools_path.path());
    let bedtools_iter =
        TsvRecordIterator::<GenomicRangeRecord<Stats>>::new(bedtools_path.path()).expect("HERE");
    let mut bedtools_gr = GRanges::from_iter(bedtools_iter, &genome).unwrap();

    let granges_iter = TsvRecordIterator::<GenomicRangeRecord<Stats>>::new(
        granges_output_file.path().to_path_buf(),
    )
    .unwrap();
    let mut granges_gr = GRanges::from_iter(granges_iter, &genome).unwrap();

    let granges_data = granges_gr.take_data().unwrap();
    let bedtools_data = bedtools_gr.take_data().unwrap();

    granges_data
        .iter()
        .zip(bedtools_data.iter())
        .for_each(|(gr, bd)| {
            assert_option_float_tol!(gr.min, bd.min);
            assert_option_float_tol!(gr.max, bd.max);

            // NOTE: this breaks because with bedools sum,
            // zero overlapping ranges = '.' (None), not 0.0.
            // Hence, our sum above is sum-not-empty
            assert_option_float_tol!(gr.sum, bd.sum);

            assert_option_float_tol!(gr.median, bd.median);
            assert_option_float_tol!(gr.mean, bd.mean);
        });
}
