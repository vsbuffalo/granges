//! Command functions that implement each of the `granges` subcommands.

use std::path::PathBuf;

use crate::{
    io::{parsers::GenomicRangesParser, OutputStream},
    prelude::*,
    ranges::{operations::adjust_range, GenomicRangeEmptyRecord, GenomicRangeRecord},
    reporting::{CommandOutput, Report},
    test_utilities::random_granges,
    traits::TsvSerialize,
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

    // Setup Output stream -- header is None for now (TODO).
    let output_stream = output.map_or(OutputStream::new_stdout(None), |file| {
        OutputStream::new(file, None)
    });
    let mut writer = output_stream.writer()?;

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
                writeln!(writer, "{}", &range_adjusted.to_tsv())?;
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
                gr.adjust_ranges(-both, both).to_tsv(output)?
            }
            GenomicRangesParser::Bed5(iter) => {
                unimplemented!()
            }
            GenomicRangesParser::Bedlike(iter) => {
                // Note the call to try_unwrap_data() here: this is because
                // we know that the records *do* have data. Unwrapping the Option<String>
                // values means that writing to TSV doesn't have to deal with this (which
                // always creates headaches).
                let gr = GRanges::from_iter(iter.try_unwrap_data(), &genome)?;
                gr.adjust_ranges(-both, both).to_tsv(output)?
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
    let mut report = Report::new();

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
            semijoin.to_tsv(output)?;

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
            semijoin.to_tsv(output)?;

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
            semijoin.to_tsv(output)?;

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
            intersection.to_tsv(output)?;

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
                gr.flanking_ranges(left, right)?.to_tsv(output)?
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
                gr.flanking_ranges(left, right)?.to_tsv(output)?
            }
            GenomicRangesParser::Unsupported => {
                return Err(GRangesError::UnsupportedGenomicRangesFileFormat)
            }
        },
        ProcessingMode::Streaming => {
            // Setup Output stream -- header is None for now (TODO).
            let output_stream = output.map_or(OutputStream::new_stdout(None), |file| {
                OutputStream::new(file, None)
            });
            let mut writer = output_stream.writer()?;

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
                                writeln!(writer, "{}", &flanking_range.to_tsv())?;
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
                                .flanking_ranges::<GenomicRangeEmptyRecord>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writeln!(writer, "{}", &flanking_range.to_tsv())?;
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
                                writeln!(writer, "{}", &flanking_range.to_tsv())?;
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
                                .flanking_ranges::<GenomicRangeEmptyRecord>(left, right, length);
                            for flanking_range in flanking_ranges {
                                writeln!(writer, "{}", &flanking_range.to_tsv())?;
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

pub fn granges_map(
    seqlens: impl Into<PathBuf>,
    left_path: &PathBuf,
    right_path: &PathBuf,
    output: Option<&PathBuf>,
    skip_missing: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    let seqnames: Vec<String> = genome.keys().cloned().collect();
    let mut report = Report::new();

    let left_iter = GenomicRangesFile::parsing_iterator(left_path)?;
    let right_iter = GenomicRangesFile::parsing_iterator(right_path)?;

    match (left_iter, right_iter) {
        (GenomicRangesParser::Bed5(left), GenomicRangesParser::Bed5(right)) => {
            let left_gr;
            let right_gr;

            if skip_missing {
                left_gr = GRanges::from_iter(left.retain_seqnames(&seqnames), &genome)?;
                right_gr = GRanges::from_iter(right.retain_seqnames(&seqnames), &genome)?;
            } else {
                left_gr = GRanges::from_iter(left, &genome)?;
                right_gr = GRanges::from_iter(right, &genome)?;
            }

            let right_gr = right_gr.into_coitrees()?;

            let left_join = left_gr.left_overlaps(&right_gr)?;

            // TODO -- map function
            // left_join.to_tsv(output)?;

            Ok(CommandOutput::new((), report))
        }
        (GenomicRangesParser::Bed3(left), GenomicRangesParser::Bedlike(right)) => {
            todo!();
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

            let intersection = left_gr.filter_overlaps(&right_gr)?;
            intersection.to_tsv(output)?;

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

            let intersection = left_gr.filter_overlaps(&right_gr)?;
            intersection.to_tsv(output)?;

            Ok(CommandOutput::new((), report))
        }
        (GenomicRangesParser::Bedlike(left), GenomicRangesParser::Bedlike(right)) => {
            todo!();
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
            intersection.to_tsv(output)?;

            Ok(CommandOutput::new((), report))
        }
        _ => Err(GRangesError::UnsupportedGenomicRangesFileFormat),
    }
}

/// Generate a random BED-like file with genomic ranges.
pub fn granges_random_bed(
    seqlens: impl Into<PathBuf>,
    num: usize,
    output: Option<impl Into<PathBuf>>,
    sort: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    // get the genome info
    let genome = read_seqlens(seqlens)?;

    let mut gr = random_granges(&genome, num)?;

    if sort {
        gr = gr.sort();
    }

    gr.to_tsv(output)?;

    let report = Report::new();
    Ok(CommandOutput::new((), report))
}
