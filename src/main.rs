use std::{fs, path::PathBuf};

use anyhow::{bail, Context, Result};
use assembly_theory::{
    assembly::{depth, index_search, ParallelMode},
    bounds::{Bound},
    canonize::CanonizeMode,
    enumerate::EnumerateMode,
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

    /// Print the assembly index, assembly depth, number of non-overlapping
    /// isomorphic subgraph pairs, and size of the search space. Note that the
    /// search space size is nondeterministic owing to some `HashMap` details.
    #[arg(long)]
    verbose: bool,

    /// Strategy for enumerating connected, non-induced subgraphs.
    #[arg(long, value_enum, default_value_t = EnumerateMode::GrowErode)]
    enumerate: EnumerateMode,

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

    #[arg(long)]
    tree: bool,
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
        None => &[Bound::Int, Bound::VecSimple, Bound::VecSmallFrags],
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


    if cli.tree {
        use std::ffi::{OsStr, OsString};
        use std::path::Path;
        use assembly_theory::molecule::Molecule;
        use rayon::iter::{ParallelIterator, IntoParallelRefIterator, IndexedParallelIterator};

        let dataset = "gdb17_200";
        fs::create_dir_all(Path::new("eval_out").join(dataset)).expect("Create dir error");
        let paths = fs::read_dir(Path::new("data").join(dataset)).unwrap();
        let mut mol_list: Vec<Molecule> = Vec::new();
        let mut names: Vec<OsString> = Vec::new();

        for path in paths {
            let name = path.unwrap().path();
            if name.extension().and_then(OsStr::to_str) != Some("mol") {
                continue;
            }
            names.push(name.file_name().unwrap().to_os_string());
            mol_list.push(
                parse_molfile_str(
                    &fs::read_to_string(name.clone())
                        .expect(&format!("Could not read file {name:?}")),
                )
                .expect(&format!("Failed to parse {name:?}")),
            );
        }

        names.par_iter().zip(mol_list.par_iter()).for_each(|(n, mol)| {
            if n == "46.mol" {
                println!("Skipping 46");
            }
            else {
                let path = Path::new("eval_out").join(dataset).join(n);
                fs::create_dir_all(&path).expect("Create dir error");
                let index = index_search(
                    &mol,
                    cli.enumerate,
                    cli.canonize,
                    cli.parallel,
                    cli.memoize,
                    cli.kernel,            
                    boundlist,
                    true,
                    Some(&path),
                );
                println!("{:?}: MA: {} Space: {}", n, index.0, index.1);
            }
        });
        std::process::exit(1);
    }
    

    // Call index calculation with all the various options.
    let (index, num_matches, states_searched) = index_search(
        &mol,
        cli.enumerate,
        cli.canonize,
        cli.parallel,
        cli.memoize,
        cli.kernel,
        boundlist,
        cli.tree,
        None,
    );

    // Print final output, depending on --verbose.
    if cli.verbose {
        println!("Assembly Index: {index}");
        println!("Non-Overlapping Isomorphic Subgraph Pairs: {num_matches}");
        println!("Assembly States Searched: {states_searched}");
    } else {
        println!("{index}");
    }

    Ok(())
}
