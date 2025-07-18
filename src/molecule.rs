//! Graph-theoretic representation of a molecule.
//!
//! Load molecules from `.mol` files, fragment and join partial molecules, and
//! perform various graph-theoretic tasks on molecules (e.g., enumerate all
//! subgraphs, test for isomorphisms, etc.).

use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fmt::Display,
    str::FromStr,
};

use petgraph::{
    dot::Dot,
    graph::{EdgeIndex, Graph, NodeIndex},
    Undirected,
};

use crate::utils::{edge_induced_subgraph, is_subset_connected};

pub(crate) type Index = u32;
pub(crate) type MGraph = Graph<Atom, Bond, Undirected, Index>;
type EdgeSet = BTreeSet<EdgeIndex<Index>>;
type NodeSet = BTreeSet<NodeIndex<Index>>;

/// Thrown by [`Element::from_str`] if the string does not represent a valid
/// chemical element.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ParseElementError;

macro_rules! periodic_table {
    ( $(($element:ident, $name:literal),)* ) => {
        #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

/// The nodes of a [`Molecule`] graph.
///
/// Atoms contain an element and have a (currently unused) `capacity` field.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Atom {
    element: Element,
    capacity: u32,
}

impl Atom {
    /// Construct an [`Atom`] of type `element` with capacity `capacity`.
    pub fn new(element: Element, capacity: u32) -> Self {
        Self { element, capacity }
    }

    /// Return this [`Atom`]'s element.
    pub fn element(&self) -> Element {
        self.element
    }
}

/// The edges of a [`Molecule`] graph.
///
/// The `.mol` file spec describes seven types of bonds, but the assembly
/// theory literature only considers single, double, and triple bonds. Notably,
/// aromatic rings are represented by alternating single and double bonds
/// instead of the aromatic bond type.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Bond {
    Single,
    Double,
    Triple,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomOrBond {
    Atom(Atom),
    Bond(Bond),
}

/// Thrown by [`Bond::try_from`] when given anything other than a 1, 2, or 3.
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

/// A simple, loopless graph with [`Element`]s as nodes and [`Bond`]s as edges.
///
/// Assembly theory literature ignores hydrogen atoms by default. Molecules can
/// have hydrogen atoms inserted into them, but by default are constructed
/// without hydrogen atoms or bonds to hydrogen atoms.
#[derive(Debug, Clone)]
pub struct Molecule {
    graph: MGraph,
}

impl Molecule {
    /// Construct a [`Molecule`] from an existing `MGraph`.
    pub(crate) fn from_graph(g: MGraph) -> Self {
        Self { graph: g }
    }

    /// Return a representation of this molecule as an `MGraph`.
    pub(crate) fn graph(&self) -> &MGraph {
        &self.graph
    }

    /// Return a pretty-printable representation of this molecule.
    pub fn info(&self) -> String {
        let dot = Dot::new(&self.graph);
        format!("{dot:?}")
    }

    /// Return `true` iff this molecule contains self-loops or multiple edges
    /// between any pair of nodes.
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

    /// Return `true` iff this molecule comprises only one bond (of any type).
    pub fn is_basic_unit(&self) -> bool {
        self.graph.edge_count() == 1 && self.graph.node_count() == 2
    }

    /// Join this molecule with `other` on edge `on`.
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

    /// Return an iterator over all ways of partitioning this molecule into two
    /// submolecules.
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
