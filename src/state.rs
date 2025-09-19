use bit_set::BitSet;

use crate::{matches::Matches, molecule::Molecule};

pub struct State {
    fragments: Vec<BitSet>,
    clique_subgraph: Option<BitSet>,
    index: usize,
    removal_order: Vec<usize>,
    last_removed: usize,
}

impl State {
    pub fn new(mol: &Molecule) -> Self {
        Self {
            fragments: {
                let mut init = BitSet::new();
                init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
                vec![init]
            },
            clique_subgraph: None,
            index: mol.graph().edge_count() - 1,
            removal_order: Vec::new(),
            last_removed: 0,
        }
    }

    pub fn update(&self, matches: &Matches, fragments: Vec<BitSet>, mat: &(usize, usize), order_idx: usize, subgraph: &Option<BitSet>) -> Self {        
        let match_len = matches.get_frag(mat.0).len();
        
        Self {
            fragments,
            index: self.index - match_len + 1,
            removal_order: {
                let mut clone = self.removal_order.clone();
                clone.push(order_idx);
                clone
            },
            last_removed: matches.match_id(mat).unwrap(),
            clique_subgraph: {
                if *subgraph != None && match_len <= matches.clique_max_len() {
                    Some(subgraph.as_ref().unwrap().clone())
                }
                else {
                    None
                }
            }
        }
    }

    pub fn fragments(&self) -> &[BitSet] {
        &self.fragments
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn removal_order(&self) -> &Vec<usize> {
        &self.removal_order
    }

    pub fn last_removed(&self) -> usize {
        self.last_removed
    }

    pub fn use_subgraph(&self) -> bool {
        self.clique_subgraph != None
    }
}