use bit_set::BitSet;

pub struct State {
    fragments: Vec<BitSet>,
    index: usize,
    removal_order: Vec<usize>,
    largest_removed: usize,
}

impl State {
    pub fn fragments(&self) -> &Vec<BitSet> {
        &self.fragments
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn largest_removed(&self) -> usize {
        self.largest_removed
    }
}