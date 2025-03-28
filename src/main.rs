use std::fs;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use assembly_theory::assembly::{index_search, serial_index_search, Bound};
use assembly_theory::{loader, molecule::Molecule};
use clap::{Args, Parser, ValueEnum};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum Bounds {
    Log,
    IntChain,
    VecChain,
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    path: PathBuf,

    #[arg(short, long)]
    /// Print out search space, duplicate subgraphs, and assembly index
    verbose: bool,

    #[command(flatten)]
    boundgroup: Option<BoundGroup>,

    #[arg(long)]
    /// Dump out molecule graph
    molecule_info: bool,

    #[arg(long)]
    /// Disable all parallelism
    serial: bool,
}

#[derive(Args, Debug)]
#[group(required = false, multiple = false)]
struct BoundGroup {
    #[arg(long)]
    /// Run branch-and-bound index search with no bounds
    no_bounds: bool,

    #[arg(long, num_args = 1..)]
    /// Run branch-and-bound index search with only specified bounds
    bounds: Vec<Bounds>,
}

fn make_boundlist(u: &[Bounds]) -> Vec<Bound> {
    let mut boundlist = u
        .iter()
        .flat_map(|b| match b {
            Bounds::Log => vec![Bound::Log],
            Bounds::IntChain => vec![Bound::IntChain],
            Bounds::VecChain => vec![Bound::VecChainSimple, Bound::VecChainSmallFrags],
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

fn index_message(mol: &Molecule, bounds: &[Bound], verbose: bool, serial: bool) -> String {
    let (index, duplicates, space) = if serial {
        serial_index_search(mol, bounds)
    } else {
        index_search(mol, bounds)
    };
    if verbose {
        let mut message = String::new();
        message.push_str(&format!("Assembly Index: {index}\n"));
        message.push_str(&format!("Duplicate subgraph pairs: {duplicates}\n"));
        message.push_str(&format!("Search Space: {space}"));
        message
    } else {
        index.to_string()
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let molfile = fs::read_to_string(&cli.path).context("Cannot read input file.")?;
    let molecule = loader::parse_molfile_str(&molfile).context("Cannot parse molfile.")?;
    if molecule.is_malformed() {
        bail!("Bad input! Molecule has self-loops or doubled edges")
    }

    if cli.molecule_info {
        println!("{}", molecule.info());
        return Ok(());
    }

    let output = match cli.boundgroup {
        None => index_message(
            &molecule,
            &[
                Bound::IntChain,
                Bound::VecChainSimple,
                Bound::VecChainSmallFrags,
            ],
            cli.verbose,
            cli.serial,
        ),
        Some(BoundGroup {
            no_bounds: true, ..
        }) => index_message(&molecule, &[], cli.verbose, cli.serial),
        Some(BoundGroup {
            no_bounds: false,
            bounds,
        }) => index_message(&molecule, &make_boundlist(&bounds), cli.verbose, cli.serial),
    };

    println!("{output}");

    Ok(())
}
