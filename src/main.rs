use std::path::PathBuf;

use assembly::{depth, index, search_space};
use clap::{Parser, ValueEnum};
use petgraph::visit::{IntoEdges, IntoEdgesDirected};

// Molecule definition, joining operation
mod molecule;

// Data IO
mod loader;

// The hard bit: compute assembly index
mod assembly;

// Utility functions
mod utils;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum Measure {
    NaiveDepth,
    NaiveIndex,
    SearchSpace,
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

    /*let G = molecule.graph(); 
    let mut nodes: Vec<molecule::Element> = Vec::new();
    for (i,w) in  G.node_weights().enumerate() {
        println!("{i}: {:?}", w.element);
        nodes.push(w.element);
    }
    println!();
    for (i,w) in  G.edge_weights().enumerate() {
        println!("{i}: {:?}", w);
    }
    println!();
    for idx in G.edge_indices().zip(G.edge_weights()) {
        let (e1, e2) = molecule.graph().edge_endpoints(idx.0).expect("bad");
        println!("{}: {:?}, ({}, {}), ({:?}, {:?})", idx.0.index(), idx.1, e1.index(), e2.index(), nodes[e1.index()], nodes[e2.index()]);
    }
    println!();*/

    let ix = if let Some(m) = cli.measure {
        match m {
            Measure::NaiveIndex => index(&molecule),
            Measure::NaiveDepth => depth(&molecule),
            Measure::SearchSpace => search_space(&molecule),
        }
    } else {
        index(&molecule)
    };
    println!("{ix}");
    Ok(())
}
