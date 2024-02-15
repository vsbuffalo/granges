use std::path::PathBuf;

use clap::{Parser, Subcommand};
use granges::prelude::GRangesError;

pub mod commands;
use crate::commands::granges_adjust;

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

#[derive(Subcommand)]
enum Commands {
    Adjust {
        /// a TSV genome file of chromosome names and their lengths
        #[arg(long, required = true)]
        seqlens: PathBuf,
        /// an input BED-like TSV file
        #[arg(required = true)]
        bedfile: PathBuf,
        /// number of basepairs to expand the range start and end positions by
        #[arg(long)]
        both: i32,
        /// an optional output file (standard output will be used if not specified)
        output: Option<PathBuf>,
    },
}

fn run() -> Result<(), GRangesError> {
    let cli = Cli::parse();
    match &cli.command {
        Some(Commands::Adjust {
            bedfile,
            seqlens,
            both,
            output,
        }) => granges_adjust(bedfile, seqlens, both, output.as_ref()),
        None => {
            println!("{}\n", INFO);
            std::process::exit(1);
        }
    }
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
