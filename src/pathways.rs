use petgraph::{Directed, Graph};

use crate::molecule::{Index, Molecule};

pub fn assembly_pathway(mol: &Molecule) -> Graph<Molecule, (), Directed, Index> {
    Graph::new()
}