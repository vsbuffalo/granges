//! Validation against bedtools

use granges::{commands::granges_random_bed, test_utilities::granges_binary_path};
use std::process::Command;
use tempfile::NamedTempFile;

#[test]
fn test_against_bedtools_slop() {
    let random_bedfile = NamedTempFile::new().expect("Failed to create temp file");
    let random_bedfile_path = random_bedfile.path().to_path_buf();
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
        .arg(random_bedfile.path())
        .output()
        .expect("bedtools slop failed");

    let granges_output = Command::new(granges_binary_path())
        .arg("adjust")
        .arg("--seqlens")
        .arg("tests_data/hg38_seqlens.tsv")
        .arg("--both")
        .arg(width.to_string())
        .arg("--sort")
        .arg(random_bedfile.path())
        .output()
        .expect("granges adjust failed");

    assert!(bedtools_output.status.success());
    assert!(granges_output.status.success());

    // TODO
    //assert_eq!(
    //    String::from_utf8_lossy(&bedtools_output.stdout),
    //    String::from_utf8_lossy(&granges_output.stdout)
    //);
}
