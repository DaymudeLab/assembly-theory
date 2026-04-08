use std::{fs, path::PathBuf};

use anyhow::{bail, Context, Result};
use assembly_theory::{
    assembly::{depth, index_search, ParallelMode},
    bounds::Bound,
    canonize::CanonizeMode,
    kernels::KernelMode,
    loader::parse_molfile_str,
    matches::Matches,
    memoize::MemoizeMode,
    pathway::edges_as_vertex_pairs,
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

    /// Print the assembly index, assembly depth, number of edge-disjoint
    /// isomorphic subgraph pairs, and size of the search space. Note that the
    /// search space size is nondeterministic owing to some `HashMap` details.
    #[arg(long)]
    verbose: bool,

    /// Timeout duration in milliseconds after which search is stopped and the
    /// best assembly index found so far is returned, or `None` if search is
    /// run until the true assembly index is found.
    #[arg(long)]
    timeout: Option<u64>,

    /// Algorithm for graph canonization.
    #[arg(long, value_enum, default_value_t = CanonizeMode::TreeNauty)]
    canonize: CanonizeMode,

    /// Parallelization strategy for the search phase.
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

    /// Extract and print the assembly pathway (bottom-up reconstruction).
    #[arg(long)]
    extract_pathway: bool,

    /// Collect and print alternative assembly pathways that use different
    /// duplicates/remnants. Implies --extract-pathway.
    #[arg(long)]
    alternative_pathways: bool,

    /// Maximum number of alternative pathways to collect (default: 10).
    /// Only meaningful with --alternative-pathways.
    #[arg(long, default_value_t = 10)]
    max_pathways: usize,
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
    let extract = cli.extract_pathway || cli.alternative_pathways;
    let max_pathways = if cli.alternative_pathways {
        Some(cli.max_pathways)
    } else {
        None
    };
    let (
        index,
        num_matches,
        states_searched,
        pathway,
        decomposition,
        alternative_pathways,
        _alternative_decompositions,
    ) = index_search(
        &mol,
        cli.timeout,
        cli.canonize,
        cli.parallel,
        cli.memoize,
        cli.kernel,
        boundlist,
        extract,
        max_pathways,
    );

    // Print final output, depending on --verbose.
    match (cli.verbose, states_searched) {
        // Found the exact assembly index.
        (true, Some(states_searched)) => {
            println!("Assembly Index:  {index}");
            println!("Matches:         {num_matches}");
            println!("States Searched: {states_searched}");
        }
        (false, Some(_)) => {
            println!("{index}");
        }
        // Search timed out and returned an upper bound.
        (true, None) => {
            println!("Assembly Index:  <= {index} (timed out)");
            println!("Matches:         {num_matches}");
            println!("States Searched: not computed on timeout");
        }
        (false, None) => {
            println!("<= {index} (timed out)");
        }
    }

    // Print the assembly pathway if extracted.
    if let Some(ref steps) = pathway {
        println!("\nAssembly Pathway ({} steps):", steps.len());
        for (i, step) in steps.iter().enumerate() {
            println!("  Step {}: {}", i + 1, step);
        }
    }

    // Print the raw decomposition (C++ comparable) if available.
    if let Some((ref match_seq, _)) = decomposition {
        let matches = Matches::new(&mol, cli.canonize);
        println!("\nDecomposition (duplicates largest-first):");
        for &mix in match_seq.iter() {
            let (h1, h2) = matches.match_fragments(mix);
            let left = edges_as_vertex_pairs(&mol, h1);
            let right = edges_as_vertex_pairs(&mol, h2);
            println!("  Left:  {:?}", left);
            println!("  Right: {:?}", right);
        }
        let mut present = bit_set::BitSet::with_capacity(mol.graph().edge_count());
        present.extend(mol.graph().edge_indices().map(|ix| ix.index()));
        for &mix in match_seq.iter() {
            let (_, h2) = matches.match_fragments(mix);
            present.difference_with(h2);
        }
        println!("  Remnant: {:?}", edges_as_vertex_pairs(&mol, &present));

        let dup_sizes: Vec<usize> = match_seq
            .iter()
            .map(|&mix| matches.match_fragments(mix).0.len())
            .collect();
        println!("\n  Summary:");
        println!(
            "    Duplicates: {} (sizes: {:?})",
            dup_sizes.len(),
            dup_sizes
        );
        println!("    Remnant edges: {}", present.len());
    }

    // Print alternative pathways if collected. These are distinct optimal
    // decompositions (same assembly index, different set of duplicates used).
    if let Some(ref alternatives) = alternative_pathways {
        let total = alternatives.len();
        for (idx, steps) in alternatives.iter().enumerate() {
            println!(
                "\nAlternative Pathway {}/{} ({} steps):",
                idx + 1,
                total,
                steps.len()
            );
            for (i, step) in steps.iter().enumerate() {
                println!("  Step {}: {}", i + 1, step);
            }
        }
    }

    Ok(())
}
