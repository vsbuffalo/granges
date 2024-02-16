use std::path::PathBuf;

use crate::{prelude::*, PositionOffset, reporting::{CommandOutput, Report}, io::OutputFile, ranges::operations::adjust_range, test_utilities::random_granges};


/// Adjust the genomic ranges in a bedfile by some specified amount.
// NOTE: we don't do build the full GRanges objects here, for efficiency.
// But it would be a good benchmark to see how much slower that is.
pub fn granges_adjust(
    bedfile: &PathBuf,
    seqlens: &PathBuf,
    both: PositionOffset,
    output: Option<&PathBuf>,
    sort: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;

    // input iterator
    let bedlike_iterator = BedlikeIterator::new(bedfile)?;

    // output stream -- header is None for now (TODO)
    let output_stream = output.map_or(OutputFile::new_stdout(None), |file| {
        OutputFile::new(file, None)
    });
    let mut writer = output_stream.writer()?;

    // for reporting stuff to the user
    let mut report = Report::new();

    // if we need to sort, we need to accumulate ranges in memory
    enum GRangesVariant {
        Indexed(GRanges<VecRangesIndexed, Vec<String>>),
        Empty(GRanges<VecRangesEmpty, ()>),
    }

    let number_columns = bedlike_iterator.number_columns();
    let mut gr = match number_columns {
        3 => GRangesVariant::Empty(GRanges::new_vec(&genome)),
        n if n > 3 => GRangesVariant::Indexed(GRanges::new_vec(&genome)),
        _ => panic!("Unexpected number of columns"),
    };

    let mut skipped_ranges = 0;
    for record in bedlike_iterator {
        let range = record?;
        let seqname = &range.seqname;
        let length = *genome
            .get(seqname)
            .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

        let possibly_adjusted_range = adjust_range(range, -both, both, length);

        if let Some(range_adjusted) = possibly_adjusted_range {
            if !sort {
                writer.write_all(&range_adjusted.to_tsv().into_bytes())?;
            } else {
                // we need to sort, so we build up the appropriate type of GRanges
                // object, depending on if we need to hold data or not.
                match gr {
                    GRangesVariant::Empty(ref mut obj) => obj.push_range_empty(&range_adjusted.seqname, range_adjusted.start, range_adjusted.end)?,
                    GRangesVariant::Indexed(ref mut obj) => obj.push_range_with_data(&range_adjusted.seqname, range_adjusted.start, range_adjusted.end, range_adjusted.data)?,
                }
            }
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

    // if we need to sort the ranges, do that and then output them
    if sort {
        match gr {
            GRangesVariant::Empty(obj) => obj.sort().to_bed3(output)?,
            GRangesVariant::Indexed(obj) => obj.sort().to_tsv(output)?,
        }

    }
    Ok(CommandOutput::new((), report))
}

/// Generate a random BED-like file with genomic ranges.
pub fn granges_random_bed(
    seqlens: impl Into<PathBuf>,
    num: u32,
    output: Option<impl Into<PathBuf>>,
    sort: bool,
    ) -> Result<CommandOutput<()>, GRangesError> {
    // get the genome info
    let genome = read_seqlens(seqlens)?;

    let mut gr = random_granges(&genome, num.try_into().unwrap())?;

    if sort {
        gr = gr.sort();
    }

    gr.to_bed3(output)?;

    let report = Report::new();
    Ok(CommandOutput::new((), report))
}
