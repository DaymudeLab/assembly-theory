use std::{fs, path::PathBuf};

use anyhow::{bail, Context, Result};
use assembly_theory::{
    assembly::{depth, index_search, ParallelMode},
    bounds::Bound,
    canonize::CanonizeMode,
    kernels::KernelMode,
    loader::parse_molfile_str,
    memoize::MemoizeMode,
};
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

    /// Print the molecule's assembly index, number of matches (edge-disjoint
    /// isomorphic subgraph pairs), number of states searched, molecule graph
    /// structure, and assembly pathways (if --pathways is given).
    #[arg(long)]
    verbose: bool,

    /// Timeout duration in milliseconds after which search is stopped and the
    /// best assembly index found so far is returned, or `None` (default) if
    /// search is run until the true assembly index is found.
    #[arg(long)]
    timeout: Option<u64>,

    /// Algorithm for graph canonization.
    #[arg(long, value_enum, default_value_t = CanonizeMode::TreeNauty)]
    canonize: CanonizeMode,

    /// Parallelism strategy for the search phase. When parallelism is enabled, the number of
    /// states searched (visible with --verbose) and the minimum assembly pathways reconstructed
    /// (obtained with --pathways) are nondeterministic.
    #[arg(long, value_enum, default_value_t = ParallelMode::DepthOne)]
    parallel: ParallelMode,

    /// Strategy for memoizing assembly states in the search phase.
    #[arg(long, value_enum, default_value_t = MemoizeMode::CanonIndex)]
    memoize: MemoizeMode,

    /// Bounding strategies to apply in the search phase.
    #[command(flatten)]
    boundsgroup: Option<BoundsGroup>,

    /// Strategy for performing graph kernelization during the search phase.
    #[arg(long, value_enum, default_value_t = KernelMode::None)]
    kernel: KernelMode,

    /// Maximum number of minimum assembly pathways to reconstruct, or `None`
    /// (default) to disable pathway reconstruction. Set this option to 0 to
    /// reconstruct all minimum assembly pathways found during search. If not
    /// `None`, output will automatically use --verbose formatting.
    #[arg(long)]
    pathways: Option<usize>,
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
    let molfile = fs::read_to_string(&cli.molpath).context("Cannot read input file.")?;
    let mol = parse_molfile_str(&molfile).context("Cannot parse molfile.")?;
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
        println!("{}", depth(&mol));
        return Ok(());
    }

    // Handle bounding strategy CLI arguments.
    let boundlist: &[Bound] = match cli.boundsgroup {
        // By default, use a combination of the integer and vector bounds.
        None => &[Bound::Int, Bound::MatchableEdges],
        // If --no-bounds is set, do not use any bounds.
        Some(BoundsGroup {
            no_bounds: true, ..
        }) => &[],
        // Otherwise, use the bounds that were specified.
        Some(BoundsGroup {
            no_bounds: false,
            bounds,
        }) => &bounds.clone(),
    };

    // Call index calculation with all the various options.
    let (index, num_matches, states_searched, pathways) = index_search(
        &mol,
        cli.timeout,
        cli.canonize,
        cli.parallel,
        cli.memoize,
        cli.kernel,
        boundlist,
        cli.pathways,
    );

    // Print final output depending on --verbose, pathways, and timeout status.
    let timed_out = states_searched.is_none();
    if timed_out {
        if cli.verbose {
            // Timed out and verbose.
            println!("Assembly Index:  <= {index} (timed out)");
            println!("Matches:         {num_matches}");
            println!("States Searched: not computed on timeout");
            println!("Pathways Found:  not computed on timeout");
            println!("\nMolecule Graph: {}", mol.info());
        } else {
            // Timed out but not verbose.
            println!("<= {index} (timed out)");
        }
    } else {
        if cli.verbose || !cli.pathways.is_none() {
            // Did not time out AND [verbose OR at least one pathway].
            println!("Assembly Index:  {index}");
            println!("Matches:         {num_matches}");
            println!("States Searched: {}", states_searched.unwrap());
            println!("Pathways Found:  {}", pathways.len());
            println!("\nMolecule Graph: {}", mol.info());
            for (i, pathway) in pathways.iter().enumerate() {
                println!("Pathway {i}: {pathway}");
            }
        } else {
            // Did not time out and only assembly index is desired.
            println!("{index}");
        }
    }

    Ok(())
}
