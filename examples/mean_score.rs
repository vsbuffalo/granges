use granges::{join::CombinedJoinDataLeftEmpty, prelude::*};

// Our overlap data processing function.
pub fn mean_score(join_data: CombinedJoinDataLeftEmpty<Option<f64>>) -> f64 {
    // Get the "right data" -- the BED5 scores.
    let overlap_scores: Vec<f64> = join_data
        .right_data
        .into_iter()
        // filter out missing values ('.' in BED)
        .filter_map(|x| x)
        .collect();

    // calculate the mean
    let score_sum: f64 = overlap_scores.iter().sum();
    score_sum / (overlap_scores.len() as f64)
}

fn try_main() -> Result<(), granges::error::GRangesError> {
    // Mock sequence lengths (e.g. "genome" file)
    let genome = seqlens!("chr1" => 100, "chr2" => 100);

    // Create parsing iterators to the left and right BED files.
    let left_iter = Bed3Iterator::new("tests_data/bedtools/map_a.txt")?;
    let right_iter = Bed5Iterator::new("tests_data/bedtools/map_b.txt")?;

    // Filter out any ranges from chromosomes not in our genome file.
    let left_gr = GRangesEmpty::from_iter(left_iter, &genome)?;
    let right_gr = GRanges::from_iter(right_iter, &genome)?;

    // Create the "right" GRanges object.
    let right_gr = right_gr
        // Convert to interval trees.
        .into_coitrees()?
        // Extract out just the score from the additional BED5 columns.
        .map_data(|bed5_cols| bed5_cols.score)?;

    // Compute overlaps and combine scores into mean.
    let results_gr = left_gr
        .left_overlaps(&right_gr)?
        .map_joins(mean_score)?;

    results_gr.to_tsv(None::<String>, &BED_TSV)?;
    Ok(())
}

fn main() {
    try_main().unwrap();
}
