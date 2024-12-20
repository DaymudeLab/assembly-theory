use std::path::PathBuf;

use assembly::{depth, index, search_space};
use clap::{Parser, ValueEnum};

// Molecule definition, joining operation
mod molecule;

// Data IO
mod loader;

// The hard bit: compute assembly index
mod assembly;

// Utility functions
mod utils;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum Measure {
    NaiveDepth,
    NaiveIndex,
    SearchSpace,
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    path: PathBuf,

    #[arg(short, long)]
    measure: Option<Measure>,
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let molecule = loader::parse(&cli.path)?;
    let ix = if let Some(m) = cli.measure {
        match m {
            Measure::NaiveIndex => index(&molecule),
            Measure::NaiveDepth => depth(&molecule),
            Measure::SearchSpace => search_space(&molecule),
        }
    } else {
        index(&molecule)
    };
    println!("{ix}");
    Ok(())
}
