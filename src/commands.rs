//! Command functions that implement each of the `granges` subcommands.

use std::{fs::File, io, path::PathBuf};

use csv::WriterBuilder;

use crate::{
    data::{operations::Operation, SerializableDatumType},
    io::{
        parsers::{Bed5Iterator, GenomicRangesParser},
        tsv::BED_TSV,
    },
    prelude::*,
    ranges::{operations::adjust_range, GenomicRangeRecord, GenomicRangeRecordEmpty},
    reporting::{CommandOutput, Report},
    test_utilities::{random_granges, random_granges_mock_bed5},
    Position, PositionOffset,
};

/// An `enum` to indicate whether an streaming or in-memory algorithm should be used.
#[derive(Clone)]
pub enum ProcessingMode {
    Streaming,
    InMemory,
}

/// Adjusts genomic ranges in a BED file by a specified amount.
///
/// This function modifies the start and end positions of each range in the input BED file based on
/// the provided offsets, ensuring the adjusted ranges do not exceed the sequence lengths specified
/// in the `seqlens` file. If sorting is enabled, the output will be sorted based on the sequence names
/// and start positions.
///
/// # Arguments
///
/// * `bedfile` - A reference to a `PathBuf` for the input BED file.
/// * `seqlens` - A reference to a `PathBuf` for the file containing sequence lengths.
/// * `both` - A [`PositionOffset`] specifying how much to adjust the start and end positions.
/// * `output` - An optional reference to a `PathBuf` where the adjusted ranges will be written. Writes
///   to stdout if `None`.
/// * `sort` - A boolean indicating whether to sort the output.
///
/// # Returns
///
/// A `Result` wrapping [`CommandOutput<()>`] on success, or [`GRangesError`] on failure.
///
/// # Errors
///
/// Returns `GRangesError` if the input BED file or sequence lengths file cannot be read, or if
/// an adjusted range exceeds the sequence boundaries.
pub fn granges_adjust(
    bedfile: &PathBuf,
    seqlens: &PathBuf,
    both: PositionOffset,
    output: Option<&PathBuf>,
    sort: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;

    let writer_boxed: Box<dyn io::Write> = match output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(io::stdout()),
    };

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer_boxed);

    // For reporting stuff to the user.
    let mut report = Report::new();
    let mut skipped_ranges = 0;

    if !sort {
        // Create the parsing iterator, and detect which variant we need based on
        // column number of the first entry.
        let bedlike_iterator = BedlikeIterator::new(bedfile)?;

        // If we don't need to sort, use iterator-based streaming processing.
        for record in bedlike_iterator {
            let range = record?;
            let seqname = &range.seqname;
            let length = *genome
                .get(seqname)
                .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

            let possibly_adjusted_range = adjust_range(range, -both, both, length);

            if let Some(range_adjusted) = possibly_adjusted_range {
                writer.serialize(range_adjusted)?;
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
    } else {
        // If we do need to sort, build up a GRanges variant and adjust ranges through
        // the GRanges interface. Note we need to detect and build a specific iterator
        // for the filetype.

        let ranges_iter = GenomicRangesFile::parsing_iterator(bedfile)?;
        match ranges_iter {
            GenomicRangesParser::Bed3(iter) => {
                let gr = GRangesEmpty::from_iter(iter, &genome)?;
                gr.adjust_ranges(-both, both)
                    .write_to_tsv(output, &BED_TSV)?
            }
            GenomicRangesParser::Bed5(iter) => {
                let gr = GRanges::from_iter(iter, &genome)?;
                gr.adjust_ranges(-both, both)
                    .write_to_tsv(output, &BED_TSV)?
            }
            GenomicRangesParser::Bedlike(iter) => {
                // Note the call to try_unwrap_data() here: this is because
                // we know that the records *do* have data. Unwrapping the Option<String>
                // values means that writing to TSV doesn't have to deal with this (which
                // always creates headaches).
                let gr = GRanges::from_iter(iter.try_unwrap_data(), &genome)?;
                gr.adjust_ranges(-both, both)
                    .write_to_tsv(output, &BED_TSV)?
            }
            GenomicRangesParser::Unsupported => {
                return Err(GRangesError::UnsupportedGenomicRangesFileFormat)
            }
        }
    }
    Ok(CommandOutput::new((), report))
}

/// Filters genomic ranges based on overlaps with another set of ranges.
///
/// Retains only the ranges from the `left_path` file that have at least one overlap with
/// the ranges in the `right_path` file. The function can optionally skip ranges that do not
/// exist in the provided sequence lengths (e.g. a "genome" file in `bedtools` lingo).
///
/// # Arguments
///
/// * `seqlens` - A reference to a `PathBuf` for the file containing sequence lengths.
/// * `left_path` - A reference to a `PathBuf` for the input BED file containing the ranges to filter.
/// * `right_path` - A reference to a `PathBuf` for the BED file containing ranges to check for overlaps.
/// * `output` - An optional reference to a `PathBuf` where the filtered ranges will be written. Writes
///   to stdout if `None`.
/// * `skip_missing` - A boolean indicating whether to skip ranges missing in the sequence lengths file.
///
/// # Returns
///
/// A `Result` wrapping [`CommandOutput<()>`] on success, or [`GRangesError`] on failure.
///
/// # Errors
///
/// Returns [`GRangesError`] if any input file cannot be read, or if there's an issue processing the ranges.
pub fn granges_filter(
    seqlens: &PathBuf,
    left_path: &PathBuf,
    right_path: &PathBuf,
    output: Option<&PathBuf>,
    skip_missing: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    let seqnames: Vec<String> = genome.keys().cloned().collect();

    let left_iter = GenomicRangesFile::parsing_iterator(left_path)?;
    let right_iter = GenomicRangesFile::parsing_iterator(right_path)?;

    // for reporting stuff to the user
    let report = Report::new();

    match (left_iter, right_iter) {
        (GenomicRangesParser::Bed3(left), GenomicRangesParser::Bed3(right)) => {
            let left_gr;
            let right_gr;

            if skip_missing {
                left_gr = GRangesEmpty::from_iter(left.retain_seqnames(&seqnames), &genome)?;
                right_gr = GRangesEmpty::from_iter(right.retain_seqnames(&seqnames), &genome)?;
            } else {
                left_gr = GRangesEmpty::from_iter(left, &genome)?;
                right_gr = GRangesEmpty::from_iter(right, &genome)?;
            }

            let right_gr = right_gr.into_coitrees()?;

            let semijoin = left_gr.filter_overlaps(&right_gr)?;
            semijoin.write_to_tsv(output, &BED_TSV)?;

            Ok(CommandOutput::new((), report))
        }
        (GenomicRangesParser::Bed3(left), GenomicRangesParser::Bedlike(right)) => {
            let left_gr;
            let right_gr;

            if skip_missing {
                left_gr = GRangesEmpty::from_iter(left.retain_seqnames(&seqnames), &genome)?;
                right_gr = GRanges::from_iter(
                    right.try_unwrap_data().retain_seqnames(&seqnames),
                    &genome,
                )?;
            } else {
                left_gr = GRangesEmpty::from_iter(left, &genome)?;
                right_gr = GRanges::from_iter(right.try_unwrap_data(), &genome)?;
            }

            let right_gr = right_gr.into_coitrees()?;

            let semijoin = left_gr.filter_overlaps(&right_gr)?;
            semijoin.write_to_tsv(output, &BED_TSV)?;

            Ok(CommandOutput::new((), report))
        }
        (GenomicRangesParser::Bedlike(left), GenomicRangesParser::Bed3(right)) => {
            let left_gr;
            let right_gr;

            if skip_missing {
                left_gr =
                    GRanges::from_iter(left.try_unwrap_data().retain_seqnames(&seqnames), &genome)?;
                right_gr = GRangesEmpty::from_iter(right.retain_seqnames(&seqnames), &genome)?;
            } else {
                left_gr = GRanges::from_iter(left.try_unwrap_data(), &genome)?;
                right_gr = GRangesEmpty::from_iter(right, &genome)?;
            }

            let right_gr = right_gr.into_coitrees()?;

            let semijoin = left_gr.filter_overlaps(&right_gr)?;
            semijoin.write_to_tsv(output, &BED_TSV)?;

            Ok(CommandOutput::new((), report))
        }
        (GenomicRangesParser::Bedlike(left), GenomicRangesParser::Bedlike(right)) => {
            let left_gr;
            let right_gr;

            if skip_missing {
                left_gr =
                    GRanges::from_iter(left.try_unwrap_data().retain_seqnames(&seqnames), &genome)?;
                right_gr = GRanges::from_iter(
                    right.try_unwrap_data().retain_seqnames(&seqnames),
                    &genome,
                )?;
            } else {
                left_gr = GRanges::from_iter(left.try_unwrap_data(), &genome)?;
                right_gr = GRanges::from_iter(right.try_unwrap_data(), &genome)?;
            }

            let right_gr = right_gr.into_coitrees()?;

            let intersection = left_gr.filter_overlaps(&right_gr)?;
            intersection.write_to_tsv(output, &BED_TSV)?;

            Ok(CommandOutput::new((), report))
        }
        _ => Err(GRangesError::UnsupportedGenomicRangesFileFormat),
    }
}

/// Generates flanking regions for genomic ranges in a BED file.
///
/// For each range in the input BED file, this function computes the flanking regions based on
/// the specified left and right offsets. The flanking regions are constrained by the sequence lengths
/// provided in the `seqlens` file. The function supports both in-memory and streaming modes for processing.
///
/// # Arguments
///
/// * `seqlens` - A reference to a `PathBuf` for the file containing sequence lengths.
/// * `bedfile` - A reference to a `PathBuf` for the input BED file.
/// * `left` - An optional `Position` specifying the left flank size.
/// * `right` - An optional `Position` specifying the right flank size.
/// * `output` - An optional reference to a `PathBuf` for the output file. Writes to stdout if `None`.
/// * `skip_missing` - A boolean indicating whether to skip ranges missing in the sequence lengths file.
/// * `mode` - A [`ProcessingMode`] indicating whether to use in-memory or streaming processing.
///
/// # Returns
///
/// A `Result` wrapping [`CommandOutput<()>`] on success, or [`GRangesError`] on failure.
///
/// # Errors
///
/// Returns [`GRangesError`] if the input BED file or sequence lengths file cannot be read, or if there's
/// an issue generating the flanking regions.
pub fn granges_flank(
    seqlens: &PathBuf,
    bedfile: &PathBuf,
    left: Option<Position>,
    right: Option<Position>,
    output: Option<&PathBuf>,
    skip_missing: bool,
    mode: ProcessingMode,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    let seqnames: Vec<String> = genome.keys().cloned().collect();
    let ranges_iter = GenomicRangesFile::parsing_iterator(bedfile)?;

    let report = Report::new();

    match mode {
        // Note: this is kept for benchmarking, to see how costly building GRanges
        // objects is versus using streaming.
        ProcessingMode::InMemory => match ranges_iter {
            GenomicRangesParser::Bed3(iter) => {
                let gr = if skip_missing {
                    GRangesEmpty::from_iter(iter.retain_seqnames(&seqnames), &genome)?
                } else {
                    GRangesEmpty::from_iter(iter, &genome)?
                };
                gr.flanking_ranges(left, right)?
                    .write_to_tsv(output, &BED_TSV)?
            }
            GenomicRangesParser::Bed5(_iter) => {
                unimplemented!()
            }
            GenomicRangesParser::Bedlike(iter) => {
                let gr = if skip_missing {
                    GRanges::from_iter(iter.try_unwrap_data().retain_seqnames(&seqnames), &genome)?
                } else {
                    GRanges::from_iter(iter.try_unwrap_data(), &genome)?
                };
                gr.flanking_ranges(left, right)?
                    .write_to_tsv(output, &BED_TSV)?
            }
            GenomicRangesParser::Unsupported => {
                return Err(GRangesError::UnsupportedGenomicRangesFileFormat)
            }
        },
        ProcessingMode::Streaming => {
            let writer_boxed: Box<dyn io::Write> = match output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(io::stdout()),
            };

            let mut writer = WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_writer(writer_boxed);

            match ranges_iter {
                // FIXME: code redundancy. But too early now to design traits, etc.
                GenomicRangesParser::Bed3(iter) => {
                    if skip_missing {
                        for record in iter.retain_seqnames(&seqnames) {
                            let range = record?;
                            let seqname = &range.seqname;
                            let length = *genome
                                .get(seqname)
                                .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

                            let flanking_ranges = range
                                .flanking_ranges::<GenomicRangeRecord<String>>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writer.serialize(flanking_range)?;
                            }
                        }
                    } else {
                        for record in iter {
                            let range = record?;
                            let seqname = &range.seqname;
                            let length = *genome
                                .get(seqname)
                                .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

                            let flanking_ranges = range
                                .flanking_ranges::<GenomicRangeRecordEmpty>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writer.serialize(flanking_range)?;
                            }
                        }
                    }
                }
                GenomicRangesParser::Bed5(_iter) => {
                    unimplemented!()
                }
                GenomicRangesParser::Bedlike(iter) => {
                    if skip_missing {
                        for record in iter.retain_seqnames(&seqnames) {
                            let range = record?;
                            let seqname = &range.seqname;
                            let length = *genome
                                .get(seqname)
                                .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

                            let flanking_ranges = range
                                .flanking_ranges::<GenomicRangeRecord<String>>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writer.serialize(flanking_range)?;
                            }
                        }
                    } else {
                        for record in iter {
                            let range = record?;
                            let seqname = &range.seqname;
                            let length = *genome
                                .get(seqname)
                                .ok_or(GRangesError::MissingSequence(seqname.to_string()))?;

                            let flanking_ranges = range
                                .flanking_ranges::<GenomicRangeRecordEmpty>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writer.serialize(flanking_range)?;
                            }
                        }
                    }
                }
                GenomicRangesParser::Unsupported => {
                    return Err(GRangesError::UnsupportedGenomicRangesFileFormat)
                }
            }
        }
    }
    Ok(CommandOutput::new((), report))
}

/// # Developer Notes
/// This function is a great way to see GRange's methods in action.
pub fn granges_map(
    seqlens: impl Into<PathBuf>,
    left_path: &PathBuf,
    right_path: &PathBuf,
    operations: Vec<Operation>,
    output: Option<&PathBuf>,
    skip_missing: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    let seqnames: Vec<String> = genome.keys().cloned().collect();
    let report = Report::new();

    let left_iter = Bed3Iterator::new(left_path)?;
    let right_iter = Bed5Iterator::new(right_path)?;

    let left_gr;
    let right_gr;

    if skip_missing {
        left_gr = GRangesEmpty::from_iter(left_iter.retain_seqnames(&seqnames), &genome)?;
        right_gr = GRanges::from_iter(right_iter.retain_seqnames(&seqnames), &genome)?;
    } else {
        left_gr = GRangesEmpty::from_iter(left_iter, &genome)?;
        right_gr = GRanges::from_iter(right_iter, &genome)?;
    }

    if left_gr.is_empty() {
        return Err(GRangesError::NoRows);
    }
    if right_gr.is_empty() {
        return Err(GRangesError::NoRows);
    }

    let right_gr = {
        // Convert to interval trees for join.
        right_gr
            .into_coitrees()?
            // Select out the score.
            .map_data(|bed5_cols| {
                // Extract out just the score.
                bed5_cols.score
            })?
    };

    // Find the overlaps.
    let left_join_gr = left_gr.left_overlaps(&right_gr)?;

    // Process all the overlaps.
    let result_gr = left_join_gr.map_joins(|join_data| {
        // Get the "right data" -- the BED5 scores
        let mut overlap_scores: Vec<f64> = join_data
            .right_data
            .into_iter()
            // Filter out the `None` values.
            .flatten()
            .collect();

        // Run all operations on the scores.
        operations
            .iter()
            .map(|operation| {
                operation
                    .run(&mut overlap_scores)
                    .into_serializable(&BED_TSV)
            })
            .collect::<Vec<SerializableDatumType>>()
    })?;

    result_gr.write_to_tsv(output, &BED_TSV)?;

    Ok(CommandOutput::new((), report))
}

/// Generate a BED3 file of genomic windows.
pub fn granges_windows(
    seqlens: impl Into<PathBuf>,
    width: Position,
    step: Option<Position>,
    chop: bool,
    output: Option<impl Into<PathBuf>>,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    GRangesEmpty::from_windows(&genome, width, step, chop)?.write_to_tsv(output, &BED_TSV)?;
    let report = Report::new();
    Ok(CommandOutput::new((), report))
}

/// Generate a random BED-like file with genomic ranges.
pub fn granges_random_bed(
    seqlens: impl Into<PathBuf>,
    num: usize,
    output: Option<impl Into<PathBuf>>,
    sort: bool,
    bed5: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    // get the genome info
    let genome = read_seqlens(seqlens)?;

    if bed5 {
        let mut gr = random_granges_mock_bed5(&genome, num)?;
        if sort {
            gr = gr.sort()
        }
        gr.write_to_tsv(output, &BED_TSV)?;
    } else {
        let mut gr = random_granges(&genome, num)?;
        if sort {
            gr = gr.sort();
        }
        gr.write_to_tsv(output, &BED_TSV)?;
    };

    let report = Report::new();
    Ok(CommandOutput::new((), report))
}
