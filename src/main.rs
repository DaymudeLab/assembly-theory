use std::fs;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, ValueEnum};
use orca::assembly::{index, index_search, Bound};
use orca::{loader, molecule::Molecule};

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

fn verbose_index(mol: &Molecule, bounds: &[Bound]) -> String {
    let (index, duplicates, space) = index_search(mol, bounds);
    let mut message = String::new();
    message.push_str(&format!("Assembly Index: {index}\n"));
    message.push_str(&format!("Duplicate subgraph pairs: {duplicates}\n"));
    message.push_str(&format!("Search Space: {space}"));
    message
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

    let output = match (cli.verbose, cli.boundgroup) {
        (false, None) => index(&molecule).to_string(),
        (
            false,
            Some(BoundGroup {
                no_bounds: true, ..
            }),
        ) => index_search(&molecule, &[]).0.to_string(),
        (
            false,
            Some(BoundGroup {
                no_bounds: false,
                bounds,
            }),
        ) => index_search(&molecule, &make_boundlist(&bounds))
            .0
            .to_string(),
        (true, None) => verbose_index(
            &molecule,
            &make_boundlist(&[Bounds::IntChain, Bounds::VecChain]),
        ),
        (
            true,
            Some(BoundGroup {
                no_bounds: true, ..
            }),
        ) => verbose_index(&molecule, &[]),
        (
            true,
            Some(BoundGroup {
                no_bounds: false,
                bounds,
            }),
        ) => verbose_index(&molecule, &make_boundlist(&bounds)),
    };

    println!("{output}");

    Ok(())
}
