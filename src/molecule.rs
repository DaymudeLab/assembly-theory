//! Graph-theoretic description of a molecule.
//!
//! This module provides functions for fragmenting and joining partial molecules. It is not
//! possible to manually construct a molecule. You can only construct a molecule by parsing a
//! `.mol` file
use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fmt::Display,
    str::FromStr,
};

use bit_set::BitSet;
use petgraph::{
    algo::{is_isomorphic, is_isomorphic_subgraph},
    dot::Dot,
    graph::{EdgeIndex, Graph, NodeIndex},
    visit::{EdgeCount, NodeCount},
    Undirected,
};

use crate::{
    utils::{edge_induced_subgraph, is_subset_connected},
    vf3::noninduced_subgraph_isomorphism_iter,
};

pub(crate) type Index = u32;
pub(crate) type MGraph = Graph<Atom, Bond, Undirected, Index>;
type EdgeSet = BTreeSet<EdgeIndex<Index>>;
type NodeSet = BTreeSet<NodeIndex<Index>>;

macro_rules! periodic_table {
    ( $(($element:ident, $name:literal),)* ) => {
        #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
        /// Represents a chemical element.
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

/// Thrown by [`Element::from_str`] if the string does not represent a chemical element.
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

/// Atoms are the vertices of a [`Molecule`] graph.
///
/// Atoms contain an element and have a (currently unused) `capacity` field.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Atom {
    element: Element,
    capacity: u32,
}

/// Bonds are the edges of a [`Molecule`] graph.
///
/// The `.mol` file spec describes seven types of bonds, but assembly theory literature
/// only considers single, double, and triple bonds. Notably, aromatic rings are represented
/// by alternating single and double bonds, instead of the aromatic bond type.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Bond {
    Single,
    Double,
    Triple,
}

/// Thrown when `from::<usize>()` does not recieve a 1, 2, or 3.
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

/// A molecule is a simple, loopless graph where nodes are elements and edges are bonds.
///
/// Assembly theory literature ignores hydrogen atoms by default. Molecules can hydrogen atoms
/// inserted into them, but by default are constructed without hydrogen atoms or bonds to hydrogen
/// atoms.
#[derive(Debug, Clone)]
pub struct Molecule {
    graph: MGraph,
}

impl Atom {
    /// Constructor function for an atom
    pub fn new(element: Element, capacity: u32) -> Self {
        Self { element, capacity }
    }

    /// Returns the element of a particular atom
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

    /// Return `true` if self is isomorphic to other
    pub fn is_isomorphic_to(&self, other: &Molecule) -> bool {
        is_isomorphic(&self.graph, &other.graph)
    }

    /// Return `true` if self is a subgraph of other
    pub fn is_subgraph_of(&self, other: &Molecule) -> bool {
        is_isomorphic_subgraph(&self.graph, &other.graph)
    }

    /// Return set of all subgraphs of self as an iterable data structure
    pub fn enumerate_induced_subgraphs(&self) -> impl Iterator<Item = NodeSet> {
        let mut solutions = HashSet::new();
        let remainder = BTreeSet::from_iter(self.graph.node_indices());
        self.generate_connected_induced_subgraphs(
            remainder,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        solutions.into_iter().filter(|s| !s.is_empty())
    }

    pub fn enumerate_noninduced_subgraphs(&self) -> impl Iterator<Item = BitSet> {
        let mut solutions = HashSet::new();
        let remainder = BTreeSet::from_iter(self.graph.edge_indices());
        self.generate_connected_noninduced_subgraphs(
            remainder,
            BTreeSet::new(),
            BTreeSet::new(),
            &mut solutions,
        );
        solutions.into_iter().filter_map(|s| {
            (!s.is_empty()).then_some(BitSet::from_iter(s.iter().map(|i| i.index())))
        })
    }

    /// Return `true` if self is not formed in a valid way
    ///
    /// In particular, a molecule is considered to be malformed if it contains
    /// multiple edges between the same source and destinations, or if there
    /// are edges from a source to itself
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
    fn generate_connected_induced_subgraphs(
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

            self.generate_connected_induced_subgraphs(
                remainder.clone(),
                subset.clone(),
                neighbors.clone(),
                solutions,
            );

            subset.insert(*v);
            neighbors.extend(self.graph.neighbors(*v));
            self.generate_connected_induced_subgraphs(remainder, subset, neighbors, solutions);
        } else if subset.len() > 2 {
            solutions.insert(subset);
        }
    }

    fn generate_connected_noninduced_subgraphs(
        &self,
        mut remainder: EdgeSet,
        mut subset: EdgeSet,
        mut neighbors: EdgeSet,
        solutions: &mut HashSet<EdgeSet>,
    ) {
        let candidates = if subset.is_empty() {
            remainder.clone()
        } else {
            remainder.intersection(&neighbors).cloned().collect()
        };

        if let Some(e) = candidates.first() {
            remainder.remove(e);

            self.generate_connected_noninduced_subgraphs(
                remainder.clone(),
                subset.clone(),
                neighbors.clone(),
                solutions,
            );

            subset.insert(*e);
            let (src, dst) = self.graph.edge_endpoints(*e).expect("malformed input");
            neighbors.extend(
                self.graph
                    .neighbors(src)
                    .filter_map(|n| self.graph.find_edge(src, n)),
            );
            neighbors.extend(
                self.graph
                    .neighbors(dst)
                    .filter_map(|n| self.graph.find_edge(dst, n)),
            );

            self.generate_connected_noninduced_subgraphs(remainder, subset, neighbors, solutions);
        } else if subset.len() > 1 && subset.len() < self.graph.edge_count() / 2 + 1 {
            solutions.insert(subset);
        }
    }

    #[allow(dead_code)]
    fn print_edgelist(&self, list: &[EdgeIndex], name: &str) {
        println!(
            "{name}: {:?}",
            list.iter()
                .map(|e| (
                    e.index(),
                    self.graph
                        .edge_endpoints(*e)
                        .map(|(i, j)| (
                            i.index(),
                            self.graph.node_weight(i).unwrap().element(),
                            j.index(),
                            self.graph.node_weight(j).unwrap().element(),
                        ))
                        .unwrap()
                ))
                .collect::<Vec<_>>()
        );
    }

    /// Return an iterator of bitsets from self containing all duplicate and
    /// non-overlapping pairs of isomorphic subgraphs
    pub fn matches(&self) -> impl Iterator<Item = (BitSet, BitSet)> {
        let mut matches = BTreeSet::new();
        for subgraph in self.enumerate_noninduced_subgraphs() {
            let mut h = self.graph().clone();
            h.retain_edges(|_, e| subgraph.contains(e.index()));

            let mut h_prime = self.graph.clone();
            h_prime.retain_edges(|_, e| !subgraph.contains(e.index()));
            h_prime.retain_nodes(|g, n| g.neighbors(n).count() > 0);

            for cert in noninduced_subgraph_isomorphism_iter(&h, &h_prime) {
                matches.insert(if subgraph < cert {
                    (subgraph.clone(), cert)
                } else {
                    (cert, subgraph.clone())
                });
            }
        }
        matches.into_iter()
    }

    /// Returns all ways to partition self as an iterable data structure of
    /// molecule pairs
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

    /// Constructor for Molecule using an MGraph
    pub(crate) fn from_graph(g: MGraph) -> Self {
        Self { graph: g }
    }

    pub fn is_basic_unit(&self) -> bool {
        self.graph.edge_count() == 1 && self.graph.node_count() == 2
    }

    /// Return the graph of self as an MGraph
    pub(crate) fn graph(&self) -> &MGraph {
        &self.graph
    }

    /// Pretty-printable string representation of self
    pub fn info(&self) -> String {
        let dot = Dot::new(&self.graph);
        format!("{:?}", dot)
    }
}

mod tests {
    #![allow(unused_imports)]
    use super::*;

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
