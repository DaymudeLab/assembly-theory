use petgraph::graph::Graph;

pub enum Element {
    Hydrogen,
    Carbon,
    Nitrogen,
    Oxygen,
}

pub struct Atom {
    element: Element,
    capacity: u32,
}

pub enum Bond {
    Single,
    Double,
}

pub struct Molecule {
    graph: Graph<Atom, Bond>,
}

impl Molecule {
    pub fn join(&self, other: &Molecule) -> Option<Molecule> {
        None
    }

    pub fn from_graph(g: Graph<Atom, Bond>) -> Self {
        Self { graph: g }
    }
}
