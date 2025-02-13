use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fmt::Display,
    str::FromStr,
};

use fixedbitset::FixedBitSet;
use petgraph::{
    algo::{is_isomorphic, is_isomorphic_subgraph, subgraph_isomorphisms_iter},
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::utils::{edge_induced_subgraph, edges_contained_within, is_subset_connected};

pub(crate) type Index = u32;
pub(crate) type MGraph = Graph<Atom, Bond, Undirected, Index>;
type MSubgraph = Graph<Atom, Option<Bond>, Undirected, Index>;
type EdgeSet = BTreeSet<EdgeIndex<Index>>;
type NodeSet = BTreeSet<NodeIndex<Index>>;

macro_rules! periodic_table {
    ( $(($element:ident, $name:literal),)* ) => {
        #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
        pub enum Element {
            $( $element, )*
        }

        impl Display for Element {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match &self {
                    $( Element::$element => write!(f, "{}", $name), )*
                }
            }
        }

        impl FromStr for Element {
            type Err = ParseElementError;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s {
                    $( $name => Ok(Element::$element), )*
                    _ => Err(ParseElementError),
                }
            }
        }
    };
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ParseElementError;

periodic_table!(
    (Hydrogen, "H"),
    (Helium, "He"),
    (Lithium, "Li"),
    (Beryllium, "Be"),
    (Boron, "B"),
    (Carbon, "C"),
    (Nitrogen, "N"),
    (Oxygen, "O"),
    (Fluorine, "F"),
    (Neon, "Ne"),
    (Sodium, "Na"),
    (Magnesium, "Mg"),
    (Aluminum, "Al"),
    (Silicon, "Si"),
    (Phosphorus, "P"),
    (Sulfur, "S"),
    (Chlorine, "Cl"),
    (Argon, "Ar"),
    (Potassium, "K"),
    (Calcium, "Ca"),
    (Scandium, "Sc"),
    (Titanium, "Ti"),
    (Vanadium, "V"),
    (Chromium, "Cr"),
    (Manganese, "Mn"),
    (Iron, "Fe"),
    (Cobalt, "Co"),
    (Nickel, "Ni"),
    (Copper, "Cu"),
    (Zinc, "Zn"),
    (Gallium, "Ga"),
    (Germanium, "Ge"),
    (Arsenic, "As"),
    (Selenium, "Se"),
    (Bromine, "Br"),
    (Krypton, "Kr"),
    (Rubidium, "Rb"),
    (Strontium, "Sr"),
    (Yttrium, "Y"),
    (Zirconium, "Zr"),
    (Niobium, "Nb"),
    (Molybdenum, "Mo"),
    (Technetium, "Tc"),
    (Ruthenium, "Ru"),
    (Rhodium, "Rh"),
    (Palladium, "Pd"),
    (Silver, "Ag"),
    (Cadmium, "Cd"),
    (Indium, "In"),
    (Tin, "Sn"),
    (Antimony, "Sb"),
    (Tellurium, "Te"),
    (Iodine, "I"),
    (Xenon, "Xe"),
    (Cesium, "Cs"),
    (Barium, "Ba"),
    (Lanthanum, "La"),
    (Cerium, "Ce"),
    (Praseodymium, "Pr"),
    (Neodymium, "Nd"),
    (Promethium, "Pm"),
    (Samarium, "Sm"),
    (Europium, "Eu"),
    (Gadolinium, "Gd"),
    (Terbium, "Tb"),
    (Dysprosium, "Dy"),
    (Holmium, "Ho"),
    (Erbium, "Er"),
    (Thulium, "Tm"),
    (Ytterbium, "Yb"),
    (Lutetium, "Lu"),
    (Hafnium, "Hf"),
    (Tantalum, "Ta"),
    (Wolfram, "W"),
    (Rhenium, "Re"),
    (Osmium, "Os"),
    (Iridium, "Ir"),
    (Platinum, "Pt"),
    (Gold, "Au"),
    (Mercury, "Hg"),
    (Thallium, "Tl"),
    (Lead, "Pb"),
    (Bismuth, "Bi"),
    (Polonium, "Po"),
    (Astatine, "At"),
    (Radon, "Rn"),
    (Francium, "Fr"),
    (Radium, "Ra"),
    (Actinium, "Ac"),
    (Thorium, "Th"),
    (Protactinium, "Pa"),
    (Uranium, "U"),
    (Neptunium, "Np"),
    (Plutonium, "Pu"),
    (Americium, "Am"),
    (Curium, "Cm"),
    (Berkelium, "Bk"),
    (Californium, "Cf"),
    (Einsteinium, "Es"),
    (Fermium, "Fm"),
    (Mendelevium, "Md"),
    (Nobelium, "No"),
    (Lawrencium, "Lr"),
    (Rutherfordium, "Rf"),
    (Dubnium, "Db"),
    (Seaborgium, "Sg"),
    (Bohrium, "Bh"),
    (Hassium, "Hs"),
    (Meitnerium, "Mt"),
    (Darmstadtium, "Ds"),
    (Roentgenium, "Rg"),
    (Copernicium, "Cn"),
    (Nihonium, "Nh"),
    (Flerovium, "Fl"),
    (Moscovium, "Mc"),
    (Livermorium, "Lv"),
    (Tennessine, "Ts"),
    (Oganesson, "Og"),
);

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Atom {
    element: Element,
    capacity: u32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Bond {
    Single,
    Double,
    Triple,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ParseBondError;

impl TryFrom<usize> for Bond {
    type Error = ParseBondError;
    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Bond::Single),
            2 => Ok(Bond::Double),
            3 => Ok(Bond::Triple),
            _ => Err(ParseBondError),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Molecule {
    graph: MGraph,
}

impl Atom {
    pub fn new(element: Element, capacity: u32) -> Self {
        Self { element, capacity }
    }

    pub fn element(&self) -> Element {
        self.element
    }
}

impl Molecule {
    /// Join self with other on `on`
    pub fn join(
        &self,
        other: &Molecule,
        on: impl IntoIterator<Item = (NodeIndex<Index>, NodeIndex<Index>)>,
    ) -> Option<Molecule> {
        let mut output_graph = self.clone();

        let mut v_set = NodeSet::new();
        let mut io_map = HashMap::<NodeIndex<Index>, NodeIndex<Index>>::new();

        for (u, v) in on.into_iter() {
            v_set.insert(v);
            io_map.insert(v, u);
        }

        for ix in other.graph.node_indices() {
            if !v_set.contains(&ix) {
                let w = *other.graph.node_weight(ix)?;
                let out = output_graph.graph.add_node(w);
                io_map.insert(ix, out);
            }
        }

        for ix in other.graph.edge_indices() {
            let (u, v) = other.graph.edge_endpoints(ix)?;
            let um = io_map.get(&u)?;
            let vm = io_map.get(&v)?;
            let w = *other.graph.edge_weight(ix)?;

            output_graph.graph.add_edge(*um, *vm, w);
        }

        Some(output_graph)
    }

    pub fn is_isomorphic_to(&self, other: &Molecule) -> bool {
        is_isomorphic(&self.graph, &other.graph)
    }

    pub fn is_subgraph_of(&self, other: &Molecule) -> bool {
        is_isomorphic_subgraph(&self.graph, &other.graph)
    }

    pub fn enumerate_subgraphs(&self) -> impl Iterator<Item = NodeSet> {
        let mut solutions = HashSet::new();
        let remainder = BTreeSet::from_iter(self.graph.node_indices());
        self.generate_connected_subgraphs(
            remainder,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        solutions.into_iter().filter(|s| !s.is_empty())
    }

    // Check if a molecule has self-loops or doubled edges
    pub fn is_malformed(&self) -> bool {
        let mut uniq = HashSet::new();
        !self.graph.edge_indices().all(|ix| {
            uniq.insert(ix)
                && self
                    .graph
                    .edge_endpoints(ix)
                    .is_some_and(|(src, dst)| src != dst)
        })
    }

    // From
    // https://stackoverflow.com/a/15722579
    // https://stackoverflow.com/a/15658245
    fn generate_connected_subgraphs(
        &self,
        mut remainder: NodeSet,
        mut subset: NodeSet,
        mut neighbors: NodeSet,
        solutions: &mut HashSet<NodeSet>,
    ) {
        let candidates = if subset.is_empty() {
            remainder.clone()
        } else {
            remainder.intersection(&neighbors).cloned().collect()
        };

        if let Some(v) = candidates.first() {
            remainder.remove(v);

            self.generate_connected_subgraphs(
                remainder.clone(),
                subset.clone(),
                neighbors.clone(),
                solutions,
            );

            subset.insert(*v);
            neighbors.extend(self.graph.neighbors(*v));
            self.generate_connected_subgraphs(remainder, subset, neighbors, solutions);
        } else if subset.len() > 2 {
            solutions.insert(subset);
        }
    }

    pub fn matches(&self) -> impl Iterator<Item = (FixedBitSet, FixedBitSet)> {
        let mut matches = BTreeSet::new();
        for subgraph in self.enumerate_subgraphs() {
            let mut h = self.graph().clone();
            h.retain_nodes(|_, n| subgraph.contains(&n));

            let h_prime = self.graph().map(
                |_, n| *n,
                |i, e| {
                    let (src, dst) = self.graph.edge_endpoints(i).unwrap();
                    (!subgraph.contains(&src) || !subgraph.contains(&dst)).then_some(*e)
                },
            );

            for cert in isomorphic_subgraphs_of(&h, &h_prime) {
                let cert = BTreeSet::from_iter(cert);

                let mut c = FixedBitSet::new();
                c.extend(edges_contained_within(&self.graph, &cert).map(|e| e.index()));

                let mut h = FixedBitSet::new();
                h.extend(edges_contained_within(&self.graph, &subgraph).map(|e| e.index()));

                matches.insert(if c < h { (c, h) } else { (h, c) });
            }
        }
        matches.into_iter()
    }

    pub fn partitions(&self) -> Option<impl Iterator<Item = (Molecule, Molecule)> + '_> {
        let mut solutions = HashSet::new();
        let remaining_edges = self.graph.edge_indices().collect();
        self.backtrack(
            remaining_edges,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        Some(solutions.into_iter().map(|(left, right)| {
            (
                Molecule {
                    graph: edge_induced_subgraph(self.graph.clone(), &left),
                },
                Molecule {
                    graph: edge_induced_subgraph(self.graph.clone(), &right),
                },
            )
        }))
    }

    fn backtrack(
        &self,
        mut remaining_edges: Vec<EdgeIndex<Index>>,
        left: EdgeSet,
        right: EdgeSet,
        solutions: &mut HashSet<(EdgeSet, EdgeSet)>,
    ) {
        if let Some(suffix) = remaining_edges.pop() {
            let mut lc = left.clone();
            lc.insert(suffix);

            let mut rc = right.clone();
            rc.insert(suffix);

            self.backtrack(remaining_edges.clone(), lc, right, solutions);
            self.backtrack(remaining_edges, left, rc, solutions);
        } else if self.is_valid_partition(&left, &right) {
            solutions.insert((left, right));
        }
    }

    fn is_valid_partition(&self, left: &EdgeSet, right: &EdgeSet) -> bool {
        !left.is_empty()
            && !right.is_empty()
            && is_subset_connected(&self.graph, left)
            && is_subset_connected(&self.graph, right)
    }

    pub fn from_graph(g: MGraph) -> Self {
        Self { graph: g }
    }

    pub fn single_bond() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        g.add_edge(u, v, Bond::Single);
        Self { graph: g }
    }

    pub fn double_bond() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        g.add_edge(u, v, Bond::Double);
        Self { graph: g }
    }

    pub fn carbonyl() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Double);
        Self { graph: g }
    }

    pub fn hydroxyl() -> Self {
        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Single);
        Self { graph: g }
    }

    pub fn is_basic_unit(&self) -> bool {
        self.is_isomorphic_to(&Molecule::hydroxyl())
            || self.is_isomorphic_to(&Molecule::carbonyl())
            || self.is_isomorphic_to(&Molecule::single_bond())
            || self.is_isomorphic_to(&Molecule::double_bond())
    }

    pub fn graph(&self) -> &MGraph {
        &self.graph
    }

    pub fn info(&self) -> String {
        let mut info = String::new();
        let g = self.graph();
        let mut nodes: Vec<Element> = Vec::new();
        for (i, w) in g.node_weights().enumerate() {
            info.push_str(&format!("{i}: {:?}\n", w.element));
            nodes.push(w.element);
        }
        info.push('\n');
        for idx in g.edge_indices().zip(g.edge_weights()) {
            let (e1, e2) = self.graph().edge_endpoints(idx.0).expect("bad");
            info.push_str(&format!(
                "{}: {:?}, ({}, {}), ({:?}, {:?})\n",
                idx.0.index(),
                idx.1,
                e1.index(),
                e2.index(),
                nodes[e1.index()],
                nodes[e2.index()]
            ));
        }

        info
    }
}

// This is bad, BIG TODO: replace with handwritten vf3
pub fn isomorphic_subgraphs_of(pattern: &MGraph, target: &MSubgraph) -> Vec<Vec<NodeIndex<Index>>> {
    if let Some(iter) =
        subgraph_isomorphisms_iter(&pattern, &target, &mut |n0, n1| n1 == n0, &mut |e0, e1| {
            e1.is_some_and(|e| *e0 == e)
        })
    {
        iter.map(|v| v.into_iter().map(NodeIndex::new).collect())
            .collect()
    } else {
        Vec::new()
    }
}

mod tests {
    #![allow(unused_imports)]
    use super::*;

    #[test]
    fn join_carbonyl_hydroxyl() {
        let c = Molecule::carbonyl();
        let h = Molecule::hydroxyl();
        let j = c
            .join(&h, [(NodeIndex::new(0), NodeIndex::new(0))])
            .unwrap();

        let mut g = Graph::default();
        let u = g.add_node(Atom {
            capacity: 4,
            element: Element::Carbon,
        });
        let v = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        let w = g.add_node(Atom {
            capacity: 2,
            element: Element::Oxygen,
        });
        g.add_edge(u, v, Bond::Single);
        g.add_edge(u, w, Bond::Double);

        assert!(j.is_isomorphic_to(&Molecule { graph: g }));
    }

    #[test]
    fn element_to_string() {
        assert!(Element::Hydrogen.to_string() == "H")
    }

    #[test]
    fn element_from_string() {
        assert!(str::parse("H") == Ok(Element::Hydrogen));
        assert!(str::parse::<Element>("Foo").is_err());
    }
}
