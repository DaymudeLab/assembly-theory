//! Graph-theoretic representation of a molecule.

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

/// Thrown by [`Element::from_str`] if the string is not a valid chemical
/// element.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ParseElementError;

macro_rules! periodic_table {
    ( $(($element:ident, $name:literal, $atomicweight:literal),)* ) => {
        #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        /// A chemical element on the periodic table.
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

        impl Element {
            pub fn repr(self) -> u8 {
                match &self {
                    $( Element::$element => $atomicweight, )*
                }
            }
        }
    };
}

periodic_table!(
    (Hydrogen, "H", 1),
    (Helium, "He", 2),
    (Lithium, "Li", 3),
    (Beryllium, "Be", 4),
    (Boron, "B", 5),
    (Carbon, "C", 6),
    (Nitrogen, "N", 7),
    (Oxygen, "O", 8),
    (Fluorine, "F", 9),
    (Neon, "Ne", 10),
    (Sodium, "Na", 11),
    (Magnesium, "Mg", 12),
    (Aluminum, "Al", 13),
    (Silicon, "Si", 14),
    (Phosphorus, "P", 15),
    (Sulfur, "S", 16),
    (Chlorine, "Cl", 17),
    (Argon, "Ar", 18),
    (Potassium, "K", 19),
    (Calcium, "Ca", 20),
    (Scandium, "Sc", 21),
    (Titanium, "Ti", 22),
    (Vanadium, "V", 23),
    (Chromium, "Cr", 24),
    (Manganese, "Mn", 25),
    (Iron, "Fe", 26),
    (Cobalt, "Co", 27),
    (Nickel, "Ni", 28),
    (Copper, "Cu", 29),
    (Zinc, "Zn", 30),
    (Gallium, "Ga", 31),
    (Germanium, "Ge", 32),
    (Arsenic, "As", 33),
    (Selenium, "Se", 34),
    (Bromine, "Br", 35),
    (Krypton, "Kr", 36),
    (Rubidium, "Rb", 37),
    (Strontium, "Sr", 38),
    (Yttrium, "Y", 39),
    (Zirconium, "Zr", 40),
    (Niobium, "Nb", 41),
    (Molybdenum, "Mo", 42),
    (Technetium, "Tc", 43),
    (Ruthenium, "Ru", 44),
    (Rhodium, "Rh", 45),
    (Palladium, "Pd", 46),
    (Silver, "Ag", 47),
    (Cadmium, "Cd", 48),
    (Indium, "In", 49),
    (Tin, "Sn", 50),
    (Antimony, "Sb", 51),
    (Tellurium, "Te", 52),
    (Iodine, "I", 53),
    (Xenon, "Xe", 54),
    (Cesium, "Cs", 55),
    (Barium, "Ba", 56),
    (Lanthanum, "La", 57),
    (Cerium, "Ce", 58),
    (Praseodymium, "Pr", 59),
    (Neodymium, "Nd", 60),
    (Promethium, "Pm", 61),
    (Samarium, "Sm", 62),
    (Europium, "Eu", 63),
    (Gadolinium, "Gd", 64),
    (Terbium, "Tb", 65),
    (Dysprosium, "Dy", 66),
    (Holmium, "Ho", 67),
    (Erbium, "Er", 68),
    (Thulium, "Tm", 69),
    (Ytterbium, "Yb", 70),
    (Lutetium, "Lu", 71),
    (Hafnium, "Hf", 72),
    (Tantalum, "Ta", 73),
    (Wolfram, "W", 74),
    (Rhenium, "Re", 75),
    (Osmium, "Os", 76),
    (Iridium, "Ir", 77),
    (Platinum, "Pt", 78),
    (Gold, "Au", 79),
    (Mercury, "Hg", 80),
    (Thallium, "Tl", 81),
    (Lead, "Pb", 82),
    (Bismuth, "Bi", 83),
    (Polonium, "Po", 84),
    (Astatine, "At", 85),
    (Radon, "Rn", 86),
    (Francium, "Fr", 87),
    (Radium, "Ra", 88),
    (Actinium, "Ac", 89),
    (Thorium, "Th", 90),
    (Protactinium, "Pa", 91),
    (Uranium, "U", 92),
    (Neptunium, "Np", 93),
    (Plutonium, "Pu", 94),
    (Americium, "Am", 95),
    (Curium, "Cm", 96),
    (Berkelium, "Bk", 97),
    (Californium, "Cf", 98),
    (Einsteinium, "Es", 99),
    (Fermium, "Fm", 100),
    (Mendelevium, "Md", 101),
    (Nobelium, "No", 102),
    (Lawrencium, "Lr", 103),
    (Rutherfordium, "Rf", 104),
    (Dubnium, "Db", 105),
    (Seaborgium, "Sg", 106),
    (Bohrium, "Bh", 107),
    (Hassium, "Hs", 108),
    (Meitnerium, "Mt", 109),
    (Darmstadtium, "Ds", 110),
    (Roentgenium, "Rg", 111),
    (Copernicium, "Cn", 112),
    (Nihonium, "Nh", 113),
    (Flerovium, "Fl", 114),
    (Moscovium, "Mc", 115),
    (Livermorium, "Lv", 116),
    (Tennessine, "Ts", 117),
    (Oganesson, "Og", 118),
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

/// Either an [`Atom`] or a [`Bond`].
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

impl Bond {
    pub fn repr(self) -> u8 {
        match self {
            Bond::Single => 1,
            Bond::Double => 2,
            Bond::Triple => 3,
        }
    }
}

/// A simple, loopless graph with [`Atom`]s as nodes and [`Bond`]s as edges.
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
    ///
    /// Only public for benchmarking purposes.
    pub fn graph(&self) -> &MGraph {
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
