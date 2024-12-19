use std::path::PathBuf;

use clap::Parser;

// Molecule definition, joining operation
mod molecule;

// Data IO
mod loader;

// Canonizer
mod canonize;

// The hard bit: compute assembly index
mod assembly;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(short = 'p', long)]
    path: Option<PathBuf>,
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();

    if let Some(path) = cli.path {
        let molecule = loader::parse(&path)?;
        println!("Final Molecule Signature: {}", canonize::canonize(&molecule));
        let index = assembly::index(&molecule);
        println!("{}", index);
    }
    Ok(())
}
