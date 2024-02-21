//! Validation against bedtools

use granges::{
    commands::granges_random_bed,
    prelude::GenomicRangesFile,
    test_utilities::{granges_binary_path, random_bedfile, temp_bedfile},
};
use std::process::Command;

#[test]
fn test_random_bed3file_filetype_detect() {
    let random_bedfile_path = temp_bedfile().path().to_path_buf();

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        100_000,
        Some(&random_bedfile_path),
        true,
    )
    .expect("could not generate random BED file");

    match GenomicRangesFile::detect(random_bedfile_path).unwrap() {
        GenomicRangesFile::Bed3(_) => (),
        _ => panic!("could not detect correct filetype"),
    }
}

#[test]
fn test_against_bedtools_slop() {
    let random_bedfile = temp_bedfile();
    let random_bedfile_path = random_bedfile.path();

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        100_000,
        Some(&random_bedfile_path),
        true,
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
    let num_ranges = 100_000;

    let random_bedfile_left_tempfile = random_bedfile(num_ranges);
    let random_bedfile_right_tempfile = random_bedfile(num_ranges);
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
    let num_ranges = 100_000;

    let random_bedfile_tempfile = random_bedfile(num_ranges);
    let random_bedfile = random_bedfile_tempfile.path();

    // for testing: uncomment and results are local for inspection
    // let random_bedfile_left = Path::new("test_left.bed");
    // let random_bedfile_right = Path::new("test_right.bed");

    granges_random_bed(
        "tests_data/hg38_seqlens.tsv",
        num_ranges,
        Some(&random_bedfile),
        true,
    )
    .expect("could not generate random BED file");

    let bedtools_outfile = temp_bedfile();
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
        .arg(bedtools_outfile.path())
        .output()
        .expect("bedtools flank failed");

    let bedtools_sorted_output = Command::new("bedtools")
        .arg("sort")
        .arg("-i")
        .arg(bedtools_outfile.path())
        .output()
        .expect("bedtools intersect failed");

    let granges_outfile = temp_bedfile();
    let granges_output = Command::new(granges_binary_path())
        .arg("flank")
        .arg("--genome")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--left")
        .arg("20")
        .arg("--right")
        .arg("30")
        .arg(&random_bedfile)
        .arg(granges_outfile.path())
        .output()
        .expect("granges flank failed");

    let granges_sorted_output = Command::new("bedtools")
        .arg("sort")
        .arg("-i")
        .arg(bedtools_outfile.path())
        .output()
        .expect("bedtools intersect failed");

    assert!(
        bedtools_sorted_output.status.success(),
        "{:?}",
        bedtools_sorted_output
    );
    assert!(
        granges_sorted_output.status.success(),
        "{:?}",
        granges_sorted_output
    );

    assert_eq!(
        String::from_utf8_lossy(&bedtools_output.stdout),
        String::from_utf8_lossy(&granges_output.stdout)
    );
}
