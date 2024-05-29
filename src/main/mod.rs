use std::path::PathBuf;

use clap::{Parser, Subcommand};
use granges::{
    commands::{
        granges_adjust, granges_filter, granges_flank, granges_map, granges_windows,
        FeatureDensity, FilterChroms, Merge, ProcessingMode,
    },
    data::operations::FloatOperation,
    prelude::GRangesError,
    Position, PositionOffset,
};

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[cfg(feature = "dev-commands")]
use granges::commands::granges_random_bed;

const INFO: &str = r#"
granges: genomic range operations built off of the GRanges library
usage: granges [--help] <subcommand>

Subcommands:
  
  adjust:             Adjust each genomic range, e.g. to add a kilobase to each end.

  filter:             Filter the left ranges based on whether they have at least one
                      overlap with a right range. This is equivalent to a filtering
                      "semi-join" in SQL terminology. 

  feature-density     Calculate the density of features per window, e.g. how many 
                      basepairs are "exon", "CDS", etc. With --exclusive, this will assign
                      each basepair exclusively to a unique "feature set", such that
                      if a basepair overlaps "CDS" and "exon" ranges, the overlapping
                      number of basepairs will be added to a new composite "CDS,exon" 
                      feature set.

  map:                Compute the left grouped overlaps between the left genomic ranges
                      and right genomic ranges, and apply one or more operations to the 
                      score column of the right BED5 file.

  merge:              Merge ranges that are within a minimum distance of each other.
          
  windows:            Create a set of genomic windows of the specified width (in 
                      basepairs), stepping the specified step size (the width, by 
                      default).
          

NOTE: granges is under active development. It is not currently meant to be
a full replacement for other genomic ranges software, such as bedtools. The
command line functionality currently used for testing and benchmarking.

Please request features and report any issues on GitHub: 
https://github.com/vsbuffalo/granges/issues

"#;

#[derive(Parser)]
#[clap(name = "granges")]
#[clap(about = INFO)]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Adjust the start, end, or both coordinates of each range by some
    /// specified amount.
    Adjust {
        /// A TSV genome file of chromosome names and their lengths
        #[arg(short, long, required = true)]
        genome: PathBuf,

        /// An input BED-like TSV file
        #[arg(required = true)]
        bedfile: PathBuf,

        /// Number of basepairs to expand the range start and end positions by
        #[arg(short, long)]
        both: PositionOffset,

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Sort the ranges after adjusting their start and end positions
        #[arg(short, long)]
        sort: bool,
        // TODO add skip_missing here
    },
    FilterChroms(FilterChroms),
    /// Filter out the left ranges that do not have overlaps with any
    /// right ranges. This is a "semi-join" in SQL terminology.
    Filter {
        /// A TSV genome file of chromosome names and their lengths
        #[arg(short, long, required = true)]
        genome: PathBuf,

        /// The "left" BED-like TSV file
        #[arg(short, long, required = true)]
        left: PathBuf,

        /// The "right" BED-like TSV file
        #[arg(short, long, required = true)]
        right: PathBuf,

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Skip ranges from sequences (e.g. chromosomes) not present in the genome file.
        /// By default, ranges with sequence names not in the genome file will raise an error.
        #[arg(short, long)]
        skip_missing: bool,
    },
    /// Compute the flanking regions for each range.
    Flank {
        /// A TSV genome file of chromosome names and their lengths
        #[arg(short, long, required = true)]
        genome: PathBuf,

        /// An input BED-like TSV file
        #[arg(required = true)]
        bedfile: PathBuf,

        /// Width (in basepairs) of flank regions to create on both sides of each range
        #[arg(short, long)]
        both: Option<Position>,

        /// Width (in basepairs) of flank regions to create on the left side of each range
        #[arg(short, long)]
        left: Option<Position>,

        /// Width (in basepairs) of flank regions to create on the right side of each range
        #[arg(short, long)]
        right: Option<Position>,

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Skip ranges from sequences (e.g. chromosomes) not present in the genome file.
        /// By default, ranges with sequence names not in the genome file will raise an error.
        #[arg(long)]
        skip_missing: bool,

        /// Processing mode
        #[arg(long)]
        in_mem: bool,
    },
    FeatureDensity(FeatureDensity),
    /// Do a "left grouped join", on the specified left and right genomic ranges,
    /// and apply one or more functions to the BED5 scores for all right genomic
    /// ranges.
    ///
    /// This is analogous to 'bedtools map'.
    Map {
        /// A TSV genome file of chromosome names and their lengths
        #[arg(short, long, required = true)]
        genome: PathBuf,

        /// The "left" BED-like TSV file
        #[arg(short, long, required = true)]
        left: PathBuf,

        /// The "right" BED-like TSV file
        #[arg(short, long, required = true)]
        right: PathBuf,

        /// Operation
        #[clap(short, long, value_parser = clap::value_parser!(FloatOperation), use_value_delimiter = true, value_delimiter = ',')]
        func: Vec<FloatOperation>,

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Skip ranges from sequences (e.g. chromosomes) not present in the genome file.
        /// By default, ranges with sequence names not in the genome file will raise an error.
        #[arg(short, long)]
        skip_missing: bool,
    },
    Merge(Merge),
    /// Create a set of genomic windows ranges using the specified width
    /// and step size, and output to BED3.
    ///
    /// If --chop is set, the "remainder" windows at the end of a chromosome
    /// that would have width less than that specified by --width are chopped
    /// off.
    ///
    /// This is analogous to 'bedtools makewindows'.
    Windows {
        /// A TSV genome file of chromosome names and their lengths
        #[arg(short, long, required = true)]
        genome: PathBuf,

        /// Width (in basepairs) of each window.
        #[arg(short, long)]
        width: Position,

        /// Step width (by default: window size).
        #[arg(short, long)]
        step: Option<Position>,

        /// If last window remainder is shorter than width, remove?
        #[arg(short, long)]
        chop: bool,

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    #[cfg(feature = "dev-commands")]
    RandomBed {
        /// a TSV genome file of chromosome names and their lengths
        #[arg(required = true)]
        genome: PathBuf,

        /// number of random ranges to generate
        #[arg(short, long, required = true)]
        num: usize,

        /// an optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// sort the ranges
        #[arg(short, long)]
        sort: bool,

        /// add an additional score columns
        #[arg(short = 'c', long)]
        scores: bool,
    },
}

fn run() -> Result<(), GRangesError> {
    let cli = Cli::parse();
    let result = match &cli.command {
        Some(Commands::Adjust {
            bedfile,
            genome,
            both,
            output,
            sort,
        }) => granges_adjust(bedfile, genome, *both, output.as_ref(), *sort),
        Some(Commands::Filter {
            genome,
            left,
            right,
            output,
            skip_missing,
        }) => granges_filter(genome, left, right, output.as_ref(), *skip_missing),
        Some(Commands::FilterChroms(filter_chroms)) => filter_chroms.run(),
        Some(Commands::Flank {
            genome,
            bedfile,
            both,
            left,
            right,
            output,
            skip_missing,
            in_mem,
        }) => {
            if both.is_some() && (left.is_some() || right.is_some()) {
                let error = clap::Error::raw(
                    clap::error::ErrorKind::ArgumentConflict,
                    "set either --both, or --left and/or --right",
                );
                return Err(error.into());
            }
            let left = left.or(*both);
            let right = right.or(*both);
            if left.is_none() && right.is_none() {
                let error = clap::Error::raw(
                    clap::error::ErrorKind::ArgumentConflict,
                    "set either --both, or --left and/or --right",
                );
                return Err(error.into());
            }
            let mode = if *in_mem {
                ProcessingMode::InMemory
            } else {
                ProcessingMode::Streaming
            };
            granges_flank(
                genome,
                bedfile,
                left,
                right,
                output.as_ref(),
                *skip_missing,
                mode,
            )
        }
        Some(Commands::Map {
            genome,
            left,
            right,
            func,
            output,
            skip_missing,
        }) => {
            if func.is_empty() {
                return Err(GRangesError::NoOperationSpecified);
            }
            granges_map(
                genome,
                left,
                right,
                func.to_vec(),
                output.as_ref(),
                *skip_missing,
            )
        }
        // NOTE: this is the new API, so clean!
        Some(Commands::FeatureDensity(density)) => density.run(),
        Some(Commands::Merge(merge)) => merge.run(),
        Some(Commands::Windows {
            genome,
            width,
            step,
            chop,
            output,
        }) => granges_windows(genome, *width, *step, *chop, output.as_ref()),
        #[cfg(feature = "dev-commands")]
        Some(Commands::RandomBed {
            genome,
            num,
            output,
            sort,
            scores,
        }) => granges_random_bed(genome, *num, output.as_ref(), *sort, *scores),
        None => {
            println!("{}\n", INFO);
            std::process::exit(1);
        }
    };
    let _output = result?;
    // TODO do something with reporting, etc.
    Ok(())
}

fn main() {
    match run() {
        Ok(_) => {}
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}
