use std::path::PathBuf;

use clap::Parser;

// Molecule definition, joining operation
mod molecule;

// Data IO
mod loader;

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
        let index = assembly::index(&molecule);
        println!("{}", index);
    }
    Ok(())
}
