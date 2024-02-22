use std::path::PathBuf;

use clap::{Parser, Subcommand};
use granges::{
    commands::{granges_adjust, granges_filter, granges_flank, ProcessingMode},
    prelude::GRangesError,
    Position, PositionOffset,
};

#[cfg(feature = "dev-commands")]
use granges::commands::granges_random_bed;

const INFO: &str = r#"
granges: genomic range operations built off of the GRanges library
usage: granges [--help] <subcommand>

Subcommands:
  
  adjust: Adjust each genomic range, e.g. to add a kilobase to each end.

  filter: Filter the left ranges based on whether they have at least one
          overlap with a right range. This is equivalent to a filtering
          "semi-join" in SQL terminology. 
          
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
        #[arg(short = 'f', long)]
        skip_missing: bool,
    },
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
        #[arg(short = 'f', long)]
        skip_missing: bool,

        /// Processing mode
        #[arg(long)]
        in_mem: bool,
    },
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

        /// An optional output file (standard output will be used if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Skip ranges from sequences (e.g. chromosomes) not present in the genome file.
        /// By default, ranges with sequence names not in the genome file will raise an error.
        #[arg(short = 'f', long)]
        skip_missing: bool,

        /// Processing mode
        #[arg(long)]
        in_mem: bool,
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
            output,
            skip_missing,
            in_mem,
        }) => {

        }
 
        #[cfg(feature = "dev-commands")]
        Some(Commands::RandomBed {
            genome,
            num,
            output,
            sort,
        }) => granges_random_bed(genome, *num, output.as_ref(), *sort),
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
