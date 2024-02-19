use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use granges::{
    commands::{granges_adjust, granges_filter},
    prelude::{GRangesError, GRanges, BedlikeIterator},
    PositionOffset, traits::{GenomicRangesOperationsExtended, RangeContainer}, Position,
};

#[cfg(feature = "dev-commands")]
use granges::commands::granges_random_bed;
use indexmap::IndexMap;

const INFO: &str = "\
granges: genomic range operations built off of the GRanges library
usage: granges [--help] <subcommand>

Subcommands:
  
  adjust: adjust each genomic range, e.g. to add a kilobase to each end.
 
";

#[derive(Parser)]
#[clap(name = "granges")]
#[clap(about = INFO)]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Args, Clone)]
#[group(required = true, multiple = false)]
enum FileType {
    /// use a BED3 file (i.e. only containing the 3 genomic range columns) 
    #[arg(long)]
    Bed3(PathBuf),

    /// use a BED-like file (a TSV that extends a BED3 file with more columns)
    #[arg(long)]
    Bedl(PathBuf),
}



impl FileType {
    pub fn to_granges<R: RangeContainer, T>(self, seqlens: &IndexMap<String, Position>) -> Result<GRanges<R, Vec<T>>, GRangesError> 
    where GRanges<R, T>: GenomicRangesOperationsExtended<R> {
        match self {
            FileType::Bed3(path) => {
                let iter = BedlikeIterator::new(path)?.drop_data();
                let gr = GRanges::from_iter(iter, seqlens);
                gr
            },
            FileType::Bedl(path) => {
                let iter = BedlikeIterator::new(path)?;
                let gr = GRanges::from_iter(iter, seqlens);
                gr
            }
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    Adjust {
        /// a TSV genome file of chromosome names and their lengths
        #[arg(long, required = true)]
        seqlens: PathBuf,
        /// an input BED-like TSV file
        #[command(flatten)]
        bedfile: PathBuf,
        /// number of basepairs to expand the range start and end positions by
        #[arg(long)]
        both: PositionOffset,
        /// an optional output file (standard output will be used if not specified)
        #[arg(long)]
        output: Option<PathBuf>,
        /// sort the ranges after adjusting their start and end positions
        #[arg(long)]
        sort: bool,
    },
    Filter {
        /// a TSV genome file of chromosome names and their lengths
        #[arg(long, required = true)]
        seqlens: PathBuf,

        /// an input BED-like TSV file
        #[command(flatten)]
        left: FileType,

        /// the "right" BED-like TSV file
        #[arg(long, required = true)]
        right: FileType,

        /// an optional output file (standard output will be used if not specified)
        #[arg(long)]
        output: Option<PathBuf>,
        /// sort the ranges after adjusting their start and end positions
        #[arg(long)]
        sort: bool,
    },

    #[cfg(feature = "dev-commands")]
    RandomBed {
        /// a TSV genome file of chromosome names and their lengths
        #[arg(required = true)]
        seqlens: PathBuf,
        /// number of random ranges to generate
        #[arg(long, required = true)]
        num: u32,
        /// an optional output file (standard output will be used if not specified)
        #[arg(long)]
        output: Option<PathBuf>,
        /// sort the ranges
        #[arg(long)]
        sort: bool,
    },
}

fn run() -> Result<(), GRangesError> {
    let cli = Cli::parse();
    let result = match &cli.command {
        Some(Commands::Adjust {
            bedfile,
            seqlens,
            both,
            output,
            sort,
        }) => granges_adjust(bedfile, seqlens, *both, output.as_ref(), *sort),
        Some(Commands::Filter {
            seqlens,
            left,
            right,
            output,
            sort,
        }) => granges_filter(seqlens, left, right, output.as_ref(), *sort),
        #[cfg(feature = "dev-commands")]
        Some(Commands::RandomBed {
            seqlens,
            num,
            output,
            sort,
        }) => granges_random_bed(seqlens, *num, output.as_ref(), *sort),
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
