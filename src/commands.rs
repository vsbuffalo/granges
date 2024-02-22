use std::path::PathBuf;

use crate::{
    io::{parsers::GenomicRangesParser, OutputFile},
    prelude::*,
    ranges::{operations::adjust_range, GenomicRangeEmptyRecord, GenomicRangeRecord},
    reporting::{CommandOutput, Report},
    test_utilities::random_granges,
    traits::TsvSerialize,
    Position, PositionOffset,
};

#[derive(Clone)]
pub enum ProcessingMode {
    Streaming,
    InMemory,
}

/// Adjust the genomic ranges in a bedfile by some specified amount.
pub fn granges_adjust(
    bedfile: &PathBuf,
    seqlens: &PathBuf,
    both: PositionOffset,
    output: Option<&PathBuf>,
    sort: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;

    // Setup Output stream -- header is None for now (TODO).
    let output_stream = output.map_or(OutputFile::new_stdout(None), |file| {
        OutputFile::new(file, None)
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

/// Retain only the ranges that have at least one overlap with another set of ranges.
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

/// Adjust the genomic ranges in a bedfile by some specified amount.
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
            let output_stream = output.map_or(OutputFile::new_stdout(None), |file| {
                OutputFile::new(file, None)
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
    skip_missing: bool)
    -> Result<CommandOutput<()>, GRangesError> {
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

                left_join.to_tsv(output)?;

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
