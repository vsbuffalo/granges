use std::path::PathBuf;

use granges::io::OutputFile;
use granges::ranges::operations::adjust_range_record;
use granges::test_utilities::random_granges;
use granges::{prelude::*, PositionOffset};

use crate::reporting::{CommandOutput, Report};

/// Adjust the genomic ranges in a bedfile by some specified amount.
// NOTE: we don't do build the full GRanges objects here, for efficiency.
// But it would be a good benchmark to see how much slower that is.
pub fn granges_adjust(
    bedfile: &PathBuf,
    seqlens: &PathBuf,
    both: PositionOffset,
    output: Option<&PathBuf>,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;

    // input iterator
    let bedlike_iterator = BedlikeIterator::new(bedfile)?;

    // output stream -- header is None for now (TODO)
    let output = output.map_or(OutputFile::new_stdout(None), |file| {
        OutputFile::new(file, None)
    });
    let mut writer = output.writer()?;

    // for reporting stuff to the user
    let mut report = Report::new();

    let mut skipped_ranges = 0;
    for record in bedlike_iterator {
        let range = record?;
        let seqname = &range.seqname;
        let length = *genome
            .get(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

        let possibly_adjusted_range = adjust_range_record(range, -both, both, length);
        if let Some(range_adjusted) = possibly_adjusted_range {
            writeln!(writer, "{}", range_adjusted.to_tsv())?;
        } else {
            skipped_ranges += 1;
        }

        if skipped_ranges > 0 {
            report.add_issue(format!(
                "{} ranges were removed because their widths after adjustment were ≤ 0",
                skipped_ranges
            ))
        }
    }
    Ok(CommandOutput::new((), report))
}

/// Generate a random BED-like file with genomic ranges.
pub fn granges_random_bed(
    seqlens: &PathBuf,
    num: u32,
    output: Option<&PathBuf>,
) -> Result<CommandOutput<()>, GRangesError> {
    // get the genome info
    let genome = read_seqlens(seqlens)?;

    let gr = random_granges(&genome, num.try_into().unwrap())?;

    gr.to_bed3(output)?;

    let report = Report::new();
    Ok(CommandOutput::new((), report))
}
