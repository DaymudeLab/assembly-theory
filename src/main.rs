use std::{fs, path::PathBuf};

use anyhow::{bail, Context, Result};
use assembly_theory::assembly::{
    ParallelMode,
    KernelMode,
    assembly_depth,
    index_search,
    serial_index_search,
};
use assembly_theory::bounds::Bound;
use assembly_theory::loader::parse_molfile_str;
use assembly_theory::molecule::{EnumerateMode, CanonizeMode};
use clap::{Args, Parser};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to .mol file to compute the assembly index for.
    molpath: PathBuf,

    /// Print molecule graph information, skipping assembly index calculation.
    #[arg(long)]
    molinfo: bool,

    /// Calculate and print the molecule's assembly depth.
    #[arg(long)]
    depth: bool,

    /// Print the assembly index, assembly depth, number of non-overlapping
    /// isomorphic subgraph pairs, and size of the search space. Note that the
    /// search space size is nondeterministic owing to some `HashMap` details.
    #[arg(long)]
    verbose: bool,

    /// Strategy for enumerating connected, non-induced subgraphs.
    #[arg(long, value_enum, default_value_t = EnumerateMode::GrowErode)]
    enumerate: EnumerateMode,

    /// Algorithm for graph canonization.
    #[arg(long, value_enum, default_value_t = CanonizeMode::Nauty)]
    canonize: CanonizeMode,

    /// Parallelization strategy for the search phase.
    #[arg(long, value_enum, default_value_t = ParallelMode::Always)]
    parallel: ParallelMode,

    /// Use dynamic programming memoization in the search phase.
    #[arg(long)]
    memoize: bool,

    /// Bounding strategies to apply in the search phase.
    #[command(flatten)]
    boundsgroup: Option<BoundsGroup>,

    /// Strategy for performing graph kernelization during the search phase.
    #[arg(long, value_enum, default_value_t = KernelMode::None)]
    kernel: KernelMode,
}

#[derive(Args, Debug)]
#[group(required = false, multiple = false)]
struct BoundsGroup {
    /// Do not use any bounding strategy during the search phase.
    #[arg(long)]
    no_bounds: bool,

    /// Apply the specified bounding strategies during the search phase.
    #[arg(long, num_args = 1..)]
    bounds: Vec<Bound>,
}

fn main() -> Result<()> {
    // Parse command line arguments.
    let cli = Cli::parse();

    // Load the .mol file as a molecule::Molecule.
    let molfile = fs::read_to_string(&cli.molpath)
        .context("Cannot read input file.")?;
    let mol = parse_molfile_str(&molfile)
        .context("Cannot parse molfile.")?;
    if mol.is_malformed() {
        bail!("Bad input! Molecule has self-loops or multi-edges.")
    }

    // If --molinfo is set, print molecule graph and exit.
    if cli.molinfo {
        println!("{}", mol.info());
        return Ok(());
    }

    // If --depth is set, calculate and print assembly depth and exit.
    if cli.depth {
        println!("{}", assembly_depth(&mol));
        return Ok(());
    }

    // TODO: Do something about EnumerateModes.
    match cli.enumerate {
        EnumerateMode::Bfs => {
            println!("WARNING: Ignoring EnumerateMode::Bfs; not implemented yet");
        }
        EnumerateMode::BfsPrune => {
            println!("WARNING: Ignoring EnumerateMode::BfsPrune; not implemented yet");
        }
        EnumerateMode::GrowErode => {
            println!("Using recursive grow-erode for subgraph enumeration");
        }
        EnumerateMode::GrowErodeIterative => {
            println!("WARNING: Ignoring EnumerateMode::GrowErodeIterative; not implemented yet");
        }
    }

    // TODO: Do something about CanonizeModes.
    match cli.canonize {
        CanonizeMode::Nauty => {
            println!("Using nauty for graph canonization");
        }
        CanonizeMode::Faulon => {
            println!("WARNING: Ignoring CanonizeMode::Faulon; not implemented yet");
        }
        CanonizeMode::TreeNauty => {
            println!("WARNING: Ignoring CanonizeMode::TreeNauty; not implemented yet");
        }
        CanonizeMode::TreeFaulon => {
            println!("WARNING: Ignoring CanonizeMode::TreeFaulon; not implemented yet");
        }
    }

    // TODO: Do something about ParallelModes.
    let parallel = match cli.parallel {
        ParallelMode::None => false,
        ParallelMode::DepthOne => {
            println!("WARNING: Ignoring ParallelMode::DepthOne; not implemented yet");
            false
        },
        ParallelMode::Always => true,
    };

    // Handle bounding strategy CLI arguments.
    let boundlist: &[Bound] = match cli.boundsgroup {
        // By default, use a combination of the integer and vector bounds.
        None => &[
            Bound::Int,
            Bound::VecSimple,
            Bound::VecSmallFrags,
        ],
        // If --no-bounds is set, do not use any bounds.
        Some(BoundsGroup {
            no_bounds: true, ..
        }) => &[],
        // Otherwise, use the bounds that were specified.
        Some(BoundsGroup {
            no_bounds: false, bounds,
        }) => &bounds.clone(),
    };

    // Call index calculation with all the various options.
    // TODO: Rework with the full list of options.
    let (index, dup_pairs, search_size) = if parallel {
        index_search(&mol, boundlist)
    } else {
        serial_index_search(&mol, boundlist)
    };

    // Print final output, depending on --verbose.
    if cli.verbose {
        println!("Assembly Index: {index}");
        println!("Non-Overlapping Isomorphic Subgraph Pairs: {dup_pairs}");
        println!("Search Space Size: {search_size}");
    } else {
        println!("{index}");
    }

    Ok(())
}
