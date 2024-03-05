//! Command functions that implement each of the `granges` subcommands.
//!
// TODO: these functions should be methods of the input struct.

use clap::Parser;
use csv::{Writer, WriterBuilder};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, Write},
    path::PathBuf,
};

use crate::{
    data::{operations::FloatOperation, SerializableDatumType},
    io::{
        parsers::{Bed5Iterator, GenomicRangesParser},
        tsv::BED_TSV,
        TsvConfig,
    },
    merging_iterators::{MergingEmptyResultIterator, MergingResultIterator},
    prelude::*,
    ranges::{operations::adjust_range, GenomicRangeRecord, GenomicRangeRecordEmpty},
    reporting::{CommandOutput, Report},
    test_utilities::{random_granges, random_granges_mock_bed5},
    unique_id::UniqueIdentifier,
    Position, PositionOffset,
};

/// Build a new TSV writer
pub fn build_tsv_writer(
    output: Option<impl Into<PathBuf>>,
) -> Result<Writer<Box<dyn Write>>, GRangesError> {
    let output = output.map(|path| path.into());
    let writer_boxed: Box<dyn io::Write> = match &output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(io::stdout()),
    };

    let writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer_boxed);

    Ok(writer)
}

/// Build a new TSV writer with config, i.e. for manual headers, metadata etc.
// TODO: use proper builder pattern here and above?
pub fn build_tsv_writer_with_config(
    output: Option<impl Into<PathBuf>>,
    config: &TsvConfig,
) -> Result<Writer<Box<dyn Write>>, GRangesError> {
    let output = output.map(|path| path.into());
    let mut writer_boxed: Box<dyn io::Write> = match &output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(io::stdout()),
    };

    if let Some(headers) = &config.headers {
        writeln!(writer_boxed, "{}", headers.join("\t"))?;
    }

    let writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer_boxed);

    Ok(writer)
}

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

    let mut writer = build_tsv_writer(output)?;

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
            GenomicRangesParser::Bed4(iter) => {
                let gr = GRanges::from_iter(iter, &genome)?;
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
    Ok(CommandOutput::new((), Some(report)))
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

            Ok(CommandOutput::new((), None))
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

            Ok(CommandOutput::new((), None))
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

            Ok(CommandOutput::new((), None))
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

            Ok(CommandOutput::new((), None))
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
            GenomicRangesParser::Bed4(iter) => {
                let gr = if skip_missing {
                    GRanges::from_iter(iter.retain_seqnames(&seqnames), &genome)?
                } else {
                    GRanges::from_iter(iter, &genome)?
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
            let mut writer = build_tsv_writer(output)?;

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
                GenomicRangesParser::Bed4(_iter) => {
                    unimplemented!()
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
    Ok(CommandOutput::new((), None))
}

/// # Developer Notes
/// This function is a great way to see GRange's methods in action.
pub fn granges_map(
    seqlens: impl Into<PathBuf>,
    left_path: &PathBuf,
    right_path: &PathBuf,
    operations: Vec<FloatOperation>,
    output: Option<&PathBuf>,
    skip_missing: bool,
) -> Result<CommandOutput<()>, GRangesError> {
    let genome = read_seqlens(seqlens)?;
    let seqnames: Vec<String> = genome.keys().cloned().collect();

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

    Ok(CommandOutput::new((), None))
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
    Ok(CommandOutput::new((), None))
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

    Ok(CommandOutput::new((), None))
}

/// Merges all the genomic ranges if they overlap by `distance`.
#[derive(Parser)]
pub struct Merge {
    /// The input BED-like TSV file to merge.
    #[arg(short, long, required = true)]
    bedfile: PathBuf,

    /// The minimum distance at which to merge ranges. Like `bedtools merge`,
    /// `--distance 0` will merge "book-ended" ranges. Negative numbers
    /// will only merge ranges that overlap by that degree of overlap.
    #[clap(short, long, default_value_t = 0)]
    distance: PositionOffset,

    ///// Whether to "group by" feature name, i.e. overlapping ranges
    ///// with different feature names will not be merged.
    //#[clap(short, long)]
    //group_features: usize,
    /// Operation to do to summarize the score column.
    #[clap(short, long, value_parser = clap::value_parser!(FloatOperation))]
    func: Option<FloatOperation>,

    /// An optional output file (standard output will be used if not specified)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

impl Merge {
    // TODO optional genome file for validation?
    pub fn run(&self) -> Result<CommandOutput<()>, GRangesError> {
        let bedfile = &self.bedfile;
        let distance = &self.distance;
        let ranges_iter = GenomicRangesFile::parsing_iterator(bedfile)?;
        let func = &self.func;

        let mut writer = build_tsv_writer(self.output.as_ref())?;

        match ranges_iter {
            GenomicRangesParser::Bed3(iter) => {
                let merging_iter = MergingEmptyResultIterator::new(iter, *distance);
                for result in merging_iter {
                    let record = result?;
                    writer.serialize(record)?;
                }
                Ok(CommandOutput::new((), None))
            }
            GenomicRangesParser::Bed4(iter) => {
                let merging_iter = MergingResultIterator::new(iter, *distance, |data| {
                    data.into_iter()
                        .map(|x| x.name)
                        .collect::<Vec<_>>()
                        .join(",")
                });
                for result in merging_iter {
                    let record = result?;
                    writer.serialize(record)?;
                }
                Ok(CommandOutput::new((), None))
            }
            GenomicRangesParser::Bed5(iter) => {
                // merging iterator, where we extract scores and apply an operation to all merged genomic ranges' scores
                let merging_iter = MergingResultIterator::new(iter, *distance, |data| {
                    let mut scores: Vec<f64> = data
                        .into_iter()
                        .filter_map(|bed5_cols| bed5_cols.score)
                        .collect();
                    // this unwrap is safe -- if func is None, we use Bed3
                    func.as_ref().unwrap().run(&mut scores)
                });

                for result in merging_iter {
                    let record = result?;
                    writer.serialize(record)?;
                }
                Ok(CommandOutput::new((), None))
            }
            GenomicRangesParser::Bedlike(_iter) => {
                todo!()
            }
            GenomicRangesParser::Unsupported => {
                Err(GRangesError::UnsupportedGenomicRangesFileFormat)
            }
        }
    }
}

/// Filter out ranges not in the specified "genome" file.
#[derive(Parser)]
pub struct FilterChroms {
    /// A TSV genome file of chromosome names and their lengths
    #[arg(short, long, required = true)]
    genome: PathBuf,

    /// The input BED file.
    #[arg(short, long, required = true)]
    bedfile: PathBuf,

    /// An optional output file (standard output will be used if not specified)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

impl FilterChroms {
    pub fn run(&self) -> Result<CommandOutput<()>, GRangesError> {
        let bedfile = &self.bedfile;
        let genome = read_seqlens(&self.genome)?;
        let bedlike_iterator = BedlikeIterator::new(bedfile)?;

        let mut writer = build_tsv_writer(self.output.as_ref())?;

        // If we don't need to sort, use iterator-based streaming processing.
        for record in bedlike_iterator {
            let range = record?;
            let seqname = &range.seqname;
            let passes_filter = genome.contains_key(seqname);
            if passes_filter {
                writer.serialize(range)?;
            }
        }

        Ok(CommandOutput::new((), None))
    }
}

// tranpose two nested vecs
// thanks to this clever solution: https://stackoverflow.com/a/64499219/147427
fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

/// Calculate the density of features in a BED4 file, per window. There are two 
/// modes:
///
///   1. Without --exclusive (default): a basepair overlapping two features will
///   increment the counts of both. 
///
///   2. With --exclusive: a basepair overlapping two features will be assigned to
///      a new composite "feature set" of the two features, and increment the 
///      count of that. In contrast to the default case, this means every overlapping
///      feature in a window is assigned exclusively to one feature set.
///
/// # Tips
/// This is useful for a quick exploratory look at feature density
/// by window, e.g. with
///
///  $ granges hist-features --bedfile hg38_ncbiRefSeq.bed.gz --width 1000 --genome 
///    hg38_seqlens.tsv --headers  | xsv table -d'\t' | less
///
/// The xsv tool (https://github.com/BurntSushi/xsv) is useful for visualizing
/// the results more clearly.
///
#[derive(Parser)]
pub struct FeatureDensity {
    /// A TSV genome file of chromosome names and their lengths
    #[arg(short, long, required = true)]
    genome: PathBuf,

    /// The input BED file.
    #[arg(short, long, required = true)]
    bedfile: PathBuf,

    /// Width (in basepairs) of each window.
    #[arg(short, long)]
    width: Position,

    /// Step width (by default: window size).
    #[arg(short, long)]
    step: Option<Position>,

    /// If last window remainder is shorter than width, remove?
    #[arg(short, long)]
    chop: bool,

    /// Assign each basepair exclusively to a single "feature set"
    /// based on what features overlap it, and count the number of
    /// these per window. E.g. a region that overlaps "CDS" and "exon"
    /// will be called a "CDS,exon" feature. Columns are sorted
    /// in descending order by column sums. Headers will always be
    /// used in this case.
    ///
    /// ⚠️  Warning! This feature is *experimental* and currently does not
    /// have unit or integration tests.
    #[arg(short, long)]
    exclusive: bool,

    /// Whether to output the results as a TSV with headers.
    /// When --classify is set, headers will always be output.
    #[arg(short = 'l', long)]
    headers: bool,

    /// An optional output file (standard output will be used if not specified)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

type GRangesFeatureMatrix = GRanges<VecRangesIndexed, Vec<Vec<Position>>>;

impl FeatureDensity {
    /// Calculate feature density per window, with non-exclusive assignment 
    /// of basepairs to features. E.g. a basepair that overlaps
    /// "CDS" and "exon" features will be added to the tallies of both.
    pub fn feature_density(&self) -> Result<(GRangesFeatureMatrix, Vec<String>), GRangesError> {
        let bedfile = &self.bedfile;
        let genome = read_seqlens(&self.genome)?;
        let bed4_iter = Bed4Iterator::new(bedfile)?;

        // Split the elements in the iterator by feature into multiple GRanges objects.
        let mut records_by_features: HashMap<String, Vec<GenomicRangeRecordEmpty>> =
            HashMap::new();
        for result in bed4_iter {
            let range = result?;
            let feature = &range.data.name;
            let vec = records_by_features.entry(feature.to_string()).or_default();
            vec.push(range.into_empty()) // drop the feature, since we're hashing on it
        }

        // Load all the *merged* split records into memory as interval trees.
        let mut gr_by_features: HashMap<String, GRangesEmpty<COITreesEmpty>> = HashMap::new();
        for (feature, ranges) in records_by_features.into_iter() {
            // doing a streaming merge of each feature; bookend merge
            let merging_iter = MergingEmptyIterator::new(ranges, 0);

            // load into memory and convert to interval trees
            let gr = GRangesEmpty::from_iter_ok(merging_iter, &genome)?.into_coitrees()?;

            assert!(!gr_by_features.contains_key(&feature));
            gr_by_features.insert(feature, gr);
        }

        // Now, create windows.
        let windows = GRangesEmpty::from_windows(&genome, self.width, self.step, self.chop)?;

        let features: Vec<_> = gr_by_features.keys().cloned().collect();
        let mut feature_matrix = Vec::new();
        for (_feature, gr) in gr_by_features.into_iter() {
            // find the features that overlap each window
            // TODO/OPTIMIZE: how costly is this clone?
            let window_counts = windows
                .clone()
                .left_overlaps(&gr)?
                .map_joins(|joins| {
                    // these are merged, so this is the number of *unique* basepairs
                    let total_overlaps: Position = joins.join.overlaps().iter().cloned().sum();
                    total_overlaps
                })?
            .take_data()?;
            feature_matrix.push(window_counts);
        }

        let feature_matrix = transpose(feature_matrix);

        // Unite the windows with the feature matrix.
        let windows = GRangesEmpty::from_windows(&genome, self.width, self.step, self.chop)?;
        let windows = windows.into_granges_data(feature_matrix)?;
        Ok((windows, features))
    }

    /// Calculate feature density, but assign each basepair exclusively to a single 
    /// "feature set", which is the unique set of features that a particular 
    /// basepair overlaps.
    pub fn feature_density_exclusive(&self) -> Result<(GRangesFeatureMatrix, Vec<String>), GRangesError> {
        let bedfile = &self.bedfile;
        let genome = read_seqlens(&self.genome)?;
        let bed4_iter = Bed4Iterator::new(bedfile)?;

        // Create GRanges of the feature indices.
        // We do this manually to avoid creating a needless data container.
        let mut gr = GRanges::new_vec_keyed(&genome);
        for result in bed4_iter {
            let range = result?;

            // Note the trick here: we are *re-using* indices.
            gr.push_range_with_key(&range.seqname, range.start, range.end, &range.data.name)?;
        }
        let gr = gr.into_coitrees()?;

        // clone the feature map
        let feature_map = gr.data().ok_or(GRangesError::NoDataContainer)?.clone();

        // Create windows.
        let windows = GRangesEmpty::from_windows(&genome, self.width, self.step, self.chop)?;

        // "Reduce" ranges and
        let windows_overlaps = windows.left_overlaps(&gr)?.map_joins(|mut join| {
            // "reduce" the ranges to a minimum spanning set that contains a vec of indices
            let ranges = join.join.reduce_ranges();
            let mut feature_overlaps: HashMap<Vec<usize>, _> = HashMap::new();
            for range in ranges.iter() {
                let mut key: Vec<usize> = range.indices().iter().map(|x| x.unwrap()).collect();
                key.sort(); // sorted key is compound key
                *feature_overlaps.entry(key).or_insert(0) += range.width();
            }

            feature_overlaps
        })?;

        // Now, we go over the data and get the observed sets
        let mut observed_sets = UniqueIdentifier::new();
        let data = windows_overlaps
            .data()
            .ok_or(GRangesError::NoDataContainer)?;
        for feature_counts in data.iter() {
            for feature_set in feature_counts.keys() {
                observed_sets.get_or_insert(feature_set);
            }
        }
        let nsets = observed_sets.len(); // total feature sets

        // Go over data again, making the sparse hashmap into a full
        // matrix of overlappng basepairs.
        let window_counts = windows_overlaps.map_data(|feature_counts| {
            // fill out matrix
            let mut counts = vec![0; nsets];
            for (i, feature_set) in observed_sets.keys().enumerate() {
                counts[i] = *feature_counts.get(feature_set).unwrap_or(&0);
            }
            counts
        })?;

        // To make this prettier, let's sort by column totals (this
        // is a little pricey; we can make option later if needed).
        let matrix = window_counts.data().ok_or(GRangesError::NoDataContainer)?;
        let col_totals = column_totals(matrix);
        let sorted_indices = sorted_indices_by_values(&col_totals);

        let window_counts = window_counts.map_data(|row| {
            let row: Vec<_> = sorted_indices.iter().map(|i| row[*i]).collect();
            row
        })?;

        // Create the labels.
        let feature_sets: Vec<String> = observed_sets
            .keys()
            .map(|indices_key| {
                let labels: Vec<_> = indices_key
                    .iter()
                    .map(|idx| feature_map.get_key(*idx).unwrap().clone())
                    .collect();
                labels.join(",")
            })
        .collect();

        // Re-order the headers too by the sorted indices.
        let feature_sets: Vec<_> = sorted_indices
            .iter()
            .map(|i| feature_sets[*i].clone())
            .collect();

        Ok((window_counts, feature_sets))
    }

    /// Run this command given the command line interface.
    pub fn run(&self) -> Result<CommandOutput<()>, GRangesError> {
        if !self.exclusive {
            let (window_counts, features) = self.feature_density()?;

            // Write everything.
            if !self.headers {
                // we have to *manually* add headers (serde needs flat structs otherwise)
                window_counts.write_to_tsv(self.output.as_ref(), &BED_TSV)?;
            } else {
                let mut headers = vec!["chrom".to_string(), "start".to_string(), "end".to_string()];
                headers.extend(features);

                let config = TsvConfig {
                    no_value_string: "NA".to_string(),
                    headers: Some(headers),
                };
                window_counts.write_to_tsv(self.output.as_ref(), &config)?;
            }
        } else {
            let (window_counts, feature_sets) = self.feature_density_exclusive()?;
            let mut headers = vec!["chrom".to_string(), "start".to_string(), "end".to_string()];
            headers.extend(feature_sets);

            let config = TsvConfig {
                no_value_string: "NA".to_string(),
                headers: Some(headers),
            };
            window_counts.write_to_tsv(self.output.as_ref(), &config)?;
        }
        Ok(CommandOutput::new((), None))
    }
}

// get column totals
fn column_totals(matrix: &Vec<Vec<Position>>) -> Vec<Position> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return Vec::new();
    }

    let num_columns = matrix[0].len();
    let mut totals = vec![0; num_columns];

    for row in matrix {
        for (i, &value) in row.iter().enumerate() {
            if i < totals.len() {
                totals[i] += value;
            }
        }
    }

    totals
}

// get descending indices order
fn sorted_indices_by_values(values: &[Position]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..values.len()).collect();
    indices.sort_by_key(|&i| std::cmp::Reverse(values[i]));
    indices
}
