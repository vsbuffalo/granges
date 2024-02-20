use std::path::PathBuf;

use clap::{Parser, Subcommand};
use granges::{
    commands::{granges_adjust, granges_filter},
    prelude::GRangesError,
    PositionOffset,
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

        /// Sort the ranges after adjusting their start and end positions
        #[arg(short, long)]
        sort: bool,

        /// Skip ranges from sequences (e.g. chromosomes) not present in the genome file.
        /// By default, ranges with sequence names not in the genome file will raise an error.
        #[arg(short = 'f', long)]
        skip_missing: bool,
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
            sort,
            skip_missing,
        }) => granges_filter(genome, left, right, output.as_ref(), *skip_missing, *sort),
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
