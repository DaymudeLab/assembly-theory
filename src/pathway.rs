//! Reconstruct minimum assembly pathways from recursive search information.

use crate::{matches::Matches, molecule::Molecule};

pub struct Pathway {}

impl Pathway {
    pub fn new(mol: &Molecule, matches: &Matches, removal_order: &Vec<usize>) -> Self {
        Self {}
    }
}
