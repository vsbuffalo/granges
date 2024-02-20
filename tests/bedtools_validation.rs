//! Validation against bedtools

use granges::{commands::granges_random_bed, test_utilities::granges_binary_path, prelude::{Bed3Iterator, GenomicRangesFile}};
use std::process::Command;
use tempfile::{NamedTempFile, Builder};


fn temp_bedfile() -> NamedTempFile {
    Builder::new()
        .suffix(".bed")
        .tempfile()
        .expect("Failed to create temp file")
}

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

#[test] fn test_against_bedtools_slop() { let random_bedfile = temp_bedfile(); let
    random_bedfile_path = random_bedfile.path();

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
        .arg("--seqlens")
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
