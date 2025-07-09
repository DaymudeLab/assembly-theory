use std::fs;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use assembly_theory::assembly::{clique_index_search, index_search, serial_index_search, Bound, Kernel};
use assembly_theory::{loader, molecule::Molecule};
use clap::{Args, Parser, ValueEnum};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum Bounds {
    Log,
    IntChain,
    VecChain,
    Weight,
    Color,
    Cover,
    Fragment,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, ValueEnum, Ord, Debug)]
pub enum KernelOption {
    Never,
    Once,
    Depth1,
    All,
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

    #[arg(long)]
    kernel_method: Option<KernelOption>,

    #[arg(long)]
    no_clique: bool,
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
            Bounds::Weight => vec![Bound::Weight],
            Bounds::Color => vec![Bound::Color],
            Bounds::Cover => vec![Bound::CoverNoSort, Bound::CoverSort],
            Bounds::Fragment => vec![Bound::Fragment],
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

fn index_message(mol: &Molecule, bounds: &[Bound], verbose: bool, serial: bool, no_clique: bool, kernel: Kernel) -> String {
    let (index, duplicates, space) = if serial {
        serial_index_search(mol, bounds)
    } else if no_clique {
        index_search(mol, bounds)
    } else {
        clique_index_search(mol, bounds, kernel)
    };
    if verbose {
        let mut message = String::new();
        message.push_str(&format!("Assembly index: {index}\n"));
        message.push_str(&format!("Duplicate subgraph pairs: {duplicates}\n"));
        message.push_str(&format!("Search space: {space}"));
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

    let kernel = match cli.kernel_method {
        Some(x) => match x {
            KernelOption::Never => Kernel::Never,
            KernelOption::Once => Kernel::Once,
            KernelOption::Depth1 => Kernel::Depth1,
            KernelOption::All => Kernel::All,
        }
        None => Kernel::Once,
    };

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
            cli.no_clique,
            kernel
        ),
        Some(BoundGroup {
            no_bounds: true, ..
        }) => index_message(&molecule, &[], cli.verbose, cli.serial, cli.no_clique, kernel),
        Some(BoundGroup {
            no_bounds: false,
            bounds,
        }) => index_message(&molecule, &make_boundlist(&bounds), cli.verbose, cli.serial, cli.no_clique, kernel),
    };

    println!("{output}");

    Ok(())
}
