use granges::prelude::*;

// Our overlap data processing function - a weighted mean.
pub fn weighted_mean_score(join_data: CombinedJoinDataLeftEmpty<Option<f64>>) -> f64 {
    // Get the widths of overlaps.
    let overlap_widths = join_data.join().overlap_widths();

    // Get the "right data" -- the BED5 scores.
    let overlap_scores: Vec<Option<f64>> = join_data
        .right_data
        .into_iter()
        .collect();

    // Calculate the weighted mean
    let mut score_sum = 0.0;
    let mut total_width = 0.0;
    for (score, width) in overlap_scores.iter().zip(overlap_widths) {
        if let Some(value) = score {
            score_sum += value * (width as f64);
            total_width += width as f64;
        } else {
            // This is a missing value, a '.' in BED.
            continue;
        }
    }
    if total_width > 0.0 {
        score_sum / total_width
    } else {
        0.0
    }
}

fn try_main() -> Result<(), granges::error::GRangesError> {
    // Load in the chromosome lengths.
    let genome = read_seqlens("tests_data/hg38_seqlens.tsv")?;

    // Create 1Mb windows, non-sliding, chopping off 
    // remainder.
    let width = 1_000_000;
    let windows = GRangesEmpty::from_windows(&genome, width, None, true)?;

    // Create a new iterator over a BED5 file.
    let scores_iter = Bed5Iterator::new("tests_data/scores.bed.gz")?;

    // Load the ranges from the iterator into a GRanges object.
    let scores = GRanges::from_iter(scores_iter, &genome)?
        // Convert to interval trees.
        .into_coitrees()?
        // Extract out just the score from the additional 
        // BED5 columns.
        .map_data(|bed5_cols| bed5_cols.score)?;

    windows
        // Do a "left overlap join" with scores.
        .left_overlaps(&scores)?
        // Process each overlapping join by applying
        // a weighted mean function.
        .map_joins(weighted_mean_score)?
        // Write to a BED file.
        .write_to_tsv(Some("tests_data/weighted_means.bed"), &BED_TSV)?;
    Ok(())
}

fn main() {
    try_main().unwrap();
}
