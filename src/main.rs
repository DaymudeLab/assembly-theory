use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use orca::assembly::{
    addition_bound, depth, index, index_and_states, log_bound, naive_index, search_space,
    vec_bound_simple, vec_bound_small_frags, Bound,
};
use orca::loader;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum Measure {
    Index,
    NaiveDepth,
    NaiveIndex,
    SearchSpace,
    MoleculeInfo,
    Naive,
    NoBounds,
    LogBound,
    AdditionBound,
    VectorBound,
    AllBounds,
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
    if molecule.is_malformed() {
        panic!("Bad input! Molecule has self-loops or doubled edges")
    }

    let output = if let Some(m) = cli.measure {
        let log_bound = Bound::Log(log_bound);
        let add_bound = Bound::Addition(addition_bound);
        let vec_simple = Bound::Vector(vec_bound_simple);
        let vec_small = Bound::Vector(vec_bound_small_frags);
        match m {
            Measure::Index => format!("{:?}", index(&molecule)),
            Measure::NaiveDepth => depth(&molecule).to_string(),
            Measure::NaiveIndex => naive_index(&molecule).to_string(),
            Measure::SearchSpace => search_space(&molecule).to_string(),
            Measure::MoleculeInfo => molecule.info(),
            Measure::Naive => naive_index(&molecule).to_string(),
            Measure::NoBounds => format!("{:?}", index_and_states(&molecule, &[])),
            Measure::LogBound => format!("{:?}", index_and_states(&molecule, &[log_bound])),
            Measure::AdditionBound => format!("{:?}", index_and_states(&molecule, &[add_bound])),
            Measure::VectorBound => format!(
                "{:?}",
                index_and_states(&molecule, &[vec_simple, vec_small])
            ),
            Measure::AllBounds => format!(
                "{:?}",
                index_and_states(&molecule, &[add_bound, vec_simple, vec_small])
            ),
        }
    } else {
        format!("{:?}", index(&molecule))
    };
    println!("{output}");
    Ok(())
}
