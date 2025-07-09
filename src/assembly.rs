//! Compute assembly indices of molecules.
//! # Example
//! ```
//! # use std::fs;
//! # use std::path::PathBuf;
//! # use assembly_theory::*;
//! # fn main() -> Result<(), std::io::Error> {
//! # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
//! // Read a molecule data file
//! let molfile = fs::read_to_string(path)?;
//! let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
//!
//! // Compute assembly index of benzene
//! assert_eq!(assembly::index(&benzene), 3);
//! # Ok(())
//! # }
//! ```
use std::{
    collections::BTreeSet, sync::{
        atomic::{AtomicUsize, Ordering::Relaxed},
        Arc,
    }
};

use bit_set::BitSet;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    molecule::Bond, molecule::Element, molecule::Molecule, utils::connected_components_under_edges
};

static PARALLEL_MATCH_SIZE_THRESHOLD: usize = 100;
static ADD_CHAIN: &[usize] = &[0, 1, 2, 2, 3, 3, 4, 3, 4, 4];

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeType {
    bond: Bond,
    ends: (Element, Element),
}

/// Enum to represent the different bounds available during the computation of molecular assembly
/// indices.
/// Bounds are used by `index_search()` to speed up assembly index computations.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Bound {
    /// `Log` bounds by the logarithm base 2 of remaining edges
    Log,
    /// `IntChain` bounds by the length of the smallest addition chain to create the remaining
    /// fragments
    IntChain,
    /// 'VecChainSimple' bounds using addition chain length with the information of the edge types
    /// in a molecule
    VecChainSimple,
    /// 'VecChainSmallFrags' bounds using information on the number of fragments of size 2 in the
    /// molecule
    VecChainSmallFrags,
    Weight,
    Color,
    CoverNoSort,
    CoverSort,
    Fragment,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kernel {
    Never,
    Once,
    Depth1,
    All
}

#[derive(Debug)]
struct CompatGraph {
    graph: Vec<BitSet>,
    weights: Vec<usize>,
    matches: Vec<(BitSet, BitSet)>,
}

impl CompatGraph {
    pub fn new(mut init_matches: Vec<(BitSet, BitSet)>) -> Self{
        let size = init_matches.len();
        init_matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

        // Initialize weights and empty graph
        let mut init_graph: Vec<BitSet> = Vec::with_capacity(size);
        let mut init_weights: Vec<usize> = Vec::with_capacity(size);
        for m in init_matches.iter() {
            init_graph.push(BitSet::with_capacity(size));
            init_weights.push(m.0.len() - 1);
        }

        // Populate graph
        for (idx1, (h1, h2)) in init_matches.iter().enumerate() {
            for (idx2, (h1p, h2p)) in init_matches[idx1 + 1..].iter().enumerate() {
                let idx2 = idx2 + idx1 + 1;

                let forward_compatible = {
                    h2.is_disjoint(h1p) && 
                    h2.is_disjoint(h2p) &&
                    (h1.is_disjoint(h1p) || h1.is_superset(h1p)) &&
                    (h1.is_disjoint(h2p) || h1.is_superset(h2p))
                };

                if forward_compatible {
                    init_graph[idx1].insert(idx2);
                    init_graph[idx2].insert(idx1);
                }
            }
        }

        Self {
            graph: init_graph,
            weights: init_weights,
            matches: init_matches,
        }
    }

    pub fn savings_ground_truth(&self, subgraph: &BitSet) -> usize {
        self.savings_ground_truth_recurse(0, 0, subgraph)
    }

    fn savings_ground_truth_recurse(&self, ix: usize, mut best: usize, subgraph: &BitSet) -> usize {
        if subgraph.len() == 0 {
            return ix;
        }
        let mut cx = ix;

        /*if ix + subgraph.iter().count() <= best && ix + self.remaining_weight_bound(&subgraph) <= best {
            return ix;
        }*/
        /*if ix + self.color_bound(&subgraph) <= best {
            return ix;
        }*/
        /*if ix + self.cover_bound(&subgraph, true) <= best{
            return ix;
        }*/

        // Search for duplicatable fragment
        for v in subgraph.iter() {
            let subgraph_clone = self.forward_neighbors(v, &subgraph);

            cx = cx.max(self.savings_ground_truth_recurse(
                ix + self.weights[v],
                best,
                &subgraph_clone,
            ));
            best = best.max(cx);
        }

        cx
    }

    pub fn len(&self) -> usize {
        self.matches.len()
    }

    pub fn degree(&self, v: usize, subgraph: &BitSet) -> usize {
        self.graph[v].intersection(subgraph).count()
    }

    pub fn neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        
        neighbors
    }

    pub fn forward_neighbors(&self, v: usize, subgraph: &BitSet) -> BitSet {
        let mut neighbors = self.graph[v].clone();
        neighbors.intersect_with(subgraph);
        let mut to_remove = vec![];
        for u in neighbors.iter() {
            if u <= v {
                to_remove.push(u);
            }
            if u > v {
                break;
            }
        }
        for u in to_remove {
            neighbors.remove(u);
        }

        neighbors
    }

    pub fn are_adjacent(&self, v: usize, u: usize) -> bool {
        self.graph[v].contains(u)
    }

    pub fn remaining_weight_bound(&self, subgraph: &BitSet) -> usize {
        let deg_sum = subgraph.iter().map(|v| self.degree(v, subgraph)).sum::<usize>() as f32;
        let max_clique = ((1_f32 + (4_f32 * deg_sum + 1_f32).sqrt()) / 2_f32).floor() as usize;
        let mut sum = 0;
        let mut iter = subgraph.iter();
        for _ in 0..max_clique {
            sum += self.weights[iter.next().unwrap()];
        };

        sum
    }

    pub fn get_match(&self, v: usize) -> &(BitSet, BitSet) {
        &self.matches[v]
    }  

    pub fn color_bound(&self, subgraph: &BitSet) -> usize{
        // Greedy coloring
        let mut colors: Vec<i32> = vec![-1; self.len()];
        let mut num_colors = 0;
        let mut largest: Vec<usize> = Vec::new();
        

        for v in (0..self.matches.len()).rev() {
            if !subgraph.contains(v) {
                continue;
            }

            let mut used: Vec<usize> = vec![0; num_colors];

            for u in subgraph.intersection(&self.graph[v]) {
                if colors[u] != -1 {
                    used[colors[u] as usize] = 1;
                }
            }

            let mut max = 0;
            let mut max_idx = num_colors;
            for i in 0..num_colors {
                if used[i] == 0 && largest[i] > max {
                    max = largest[i];
                    max_idx = i;
                }
            }

            if max_idx == num_colors {
                num_colors += 1;
                largest.push(0);
            }
            if self.weights[v] > largest[max_idx] {
                largest[max_idx] = self.weights[v]
            }

            colors[v] = max_idx as i32;
        }

        largest.iter().sum::<usize>()
    }

    pub fn cover_bound(&self, subgraph: &BitSet, sort: bool) -> usize {
        // Sort vertices
        if sort {
            let mut vertices: Vec<(usize, usize)> = Vec::with_capacity(subgraph.len());
            for v in subgraph {
                vertices.push((v, self.degree(v, subgraph)));
            }   
            vertices.sort_by(|a, b| b.1.cmp(&a.1));
            self.cover_bound_helper(subgraph, vertices.iter().map(|(v, _)| *v))
        }
        else {
            let vertices = (0..self.matches.len()).rev().filter(|v| subgraph.contains(*v));
            self.cover_bound_helper(subgraph, vertices)
        }
    }

    fn cover_bound_helper(&self, subgraph: &BitSet, iter: impl Iterator<Item = usize>) -> usize {
        let mut colors: Vec<Option<Vec<usize>>> = vec![None; self.len()];
        let mut col_weights = vec![];
        let mut num_col = 0;

        for v in iter {
            let mut v_col = Vec::new();
            let mut used = vec![0; num_col];

            // Find colors used in neighborhood of v
            for u in subgraph.intersection(&self.graph[v]) {
                let Some(u_col) = &colors[u] else {
                    continue;
                };

                for c in u_col {
                    used[*c] = 1;
                }
            }

            let mut total_weight = 0;
            let v_val = self.weights[v];
            // Find colors to give to v
            for c in 0..num_col {
                if used[c] == 1 {
                    continue;
                }

                v_col.push(c);
                total_weight += col_weights[c];

                if total_weight >= v_val {
                    break;
                }
            }

            if total_weight == 0 {
                v_col.push(num_col);
                col_weights.push(v_val);
                num_col += 1
            }
            else if total_weight < v_val {
                let mut k = num_col - 1;
                while used[k] == 1 {
                    k -= 1
                }
                col_weights[k] += v_val - total_weight
            }

            colors[v] = Some(v_col);
        };

        col_weights.iter().sum()
    }

    pub fn frag_bound(&self, subgraph: &BitSet, fragments: &[BitSet]) -> usize {
        let total_bonds = fragments.iter().map(|x| x.len()).sum::<usize>();
        let mut bound = 0;
        let sizes = {
            let mut vec = vec![];
            let mut prev = 0;
            for s in subgraph.iter().map(|v| self.weights[v] + 1) {
                if s != prev {
                    vec.push(s);
                    prev = s;
                }
            }

            vec
        };
        
        for i in sizes {
            let mut bound_temp = 0;
            let mut has_bonds = fragments.len();
            let mut num_bonds: Vec<usize> = fragments.iter().map(|x| x.len()).collect();
            let mut smallest_remove = i;

            for v in subgraph.iter() {
                if has_bonds == 0 {
                    break;
                }
                if self.weights[v] + 1 > i {
                    continue;
                }


                let dup = &self.matches[v];
                let bond = dup.1.iter().next().unwrap();
                let mut j = 0;
                while !fragments[j].contains(bond) {
                    j += 1;
                }

                if num_bonds[j] > 0 {
                    let remove = std::cmp::min(dup.0.len(), num_bonds[j]);
                    bound_temp += 1;
                    num_bonds[j] -= remove;
                    smallest_remove = std::cmp::min(smallest_remove, remove);

                    if num_bonds[j] == 0 {
                        has_bonds -= 1;
                    }
                }
            }

            let leftover = num_bonds.iter().sum::<usize>();
            let log = {
                if leftover > 0 {
                    0
                }
                else {
                    (smallest_remove as f32).log2().ceil() as usize
                }
            };
            bound = std::cmp::max(bound, total_bonds - bound_temp - leftover - log);
        }

        bound
    }
}

pub fn naive_assembly_depth(mol: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in mol.partitions().unwrap() {
        let l = if left.is_basic_unit() {
            0
        } else {
            naive_assembly_depth(&left)
        };

        let r = if right.is_basic_unit() {
            0
        } else {
            naive_assembly_depth(&right)
        };

        ix = ix.min(l.max(r) + 1)
    }
    ix
}

fn recurse_naive_index_search(
    mol: &Molecule,
    matches: &BTreeSet<(BitSet, BitSet)>,
    fragments: &[BitSet],
    ix: usize,
) -> usize {
    let mut cx = ix;
    for (h1, h2) in matches {
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            continue;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
            fractures.push(h1.clone());
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));

            fractures.push(h1.clone());
        }
        cx = cx.min(recurse_naive_index_search(
            mol,
            matches,
            &fractures,
            ix - h1.len() + 1,
        ));
    }
    cx
}

/// Calculates the assembly index of a molecule without using any bounding strategy or
/// parallelization. This function is very inefficient and should only be used as a performance
/// benchmark against other strategies.
pub fn naive_index_search(mol: &Molecule) -> u32 {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    recurse_naive_index_search(
        mol,
        &mol.matches().collect(),
        &[init],
        mol.graph().edge_count() - 1,
    ) as u32
}

#[allow(clippy::too_many_arguments)]
fn recurse_index_search(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    ix: usize,
    largest_remove: usize,
    mut best: usize,
    bounds: &[Bound],
    states_searched: &mut usize,
    last_removed: i32
) -> usize {
    let mut cx = ix;

    *states_searched += 1;

    // Branch and Bound
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::IntChain => ix - addition_bound(fragments, largest_remove) >= best,
            Bound::VecChainSimple => ix - vec_bound_simple(fragments, largest_remove, mol) >= best,
            Bound::VecChainSmallFrags => {
                ix - vec_bound_small_frags(fragments, largest_remove, mol) >= best
            }
            _ => false
        };
        if exceeds {
            return ix;
        }
    }

    // Search for duplicatable fragment
    for (i, (h1, h2)) in matches.iter().enumerate() {
        let i = i as i32;
        if i <= last_removed {continue;}
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let largest_remove = h1.len();

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            continue;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));
        }

        fractures.retain(|i| i.len() > 1);
        fractures.push(h1.clone());

        cx = cx.min(recurse_index_search(
            mol,
            &matches,
            &fractures,
            ix - h1.len() + 1,
            largest_remove,
            best,
            bounds,
            states_searched,
            i
        ));
        best = best.min(cx);
    }

    cx
}


fn recurse_clique_index_search(mol: &Molecule,
    fragments: &[BitSet],
    ix: usize,
    mut best: usize,
    bounds: &[Bound],
    states_searched: &mut usize,
    subgraph: BitSet,
    matches_graph: &CompatGraph,
    depth: usize,
    must_include: &Vec<usize>,
    kernel_method: &Kernel,
) -> usize {
    if subgraph.len() == 0 {
        return ix;
    }
    let mut cx = ix;
    let largest_remove = matches_graph.weights[subgraph.iter().next().unwrap()] + 1;
    *states_searched += 1;

    // Bounds
    for bound_type in bounds {
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::IntChain => ix - addition_bound(fragments, largest_remove) >= best,
            Bound::VecChainSimple => ix - vec_bound_simple(fragments, largest_remove, mol) >= best,
            Bound::VecChainSmallFrags => {
                ix - vec_bound_small_frags(fragments, largest_remove, mol) >= best
            },
            Bound::Weight => {
                ix >= best + subgraph.iter().count() && 
                ix >= best + matches_graph.remaining_weight_bound(&subgraph)
            },
            Bound::Color => ix >= best + matches_graph.color_bound(&subgraph),
            Bound::CoverNoSort => ix >= best + matches_graph.cover_bound(&subgraph, false),
            Bound::CoverSort => ix >= best + matches_graph.cover_bound(&subgraph, true),
            Bound::Fragment => ix >= best + matches_graph.frag_bound(&subgraph, fragments),
        };
        if exceeds {
            return ix;
        }
    }

    // Search for duplicatable fragment
    for v in subgraph.iter() {
        let (h1, h2) = matches_graph.get_match(v);

        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            continue;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));
        }

        fractures.retain(|i| i.len() > 1);
        fractures.push(h1.clone());

        let mut subgraph_clone = matches_graph.forward_neighbors(v, &subgraph);
        let mut must_include_clone = must_include.clone();

        // Kernelize
        if *kernel_method == Kernel::All || (*kernel_method == Kernel:: Depth1 && depth == 1) {
            subgraph_clone = deletion_kernel(matches_graph, subgraph_clone);
            must_include_clone.append(&mut inclusion_kernel(matches_graph, &subgraph_clone));
        }

        cx = cx.min(recurse_clique_index_search(
            mol,
            &fractures,
            ix - matches_graph.weights[v],
            best,
            bounds,
            states_searched,
            subgraph_clone,
            matches_graph,
            depth + 1,
            &must_include_clone,
            kernel_method,
        ));
        best = best.min(cx);

        if must_include.contains(&v) {
            return cx;
        }
    }

    cx
}

pub fn clique_index_search(mol: &Molecule, bounds: &[Bound], kernel_method: Kernel) -> (u32, u32, usize) {
    // Graph Initialization
    let matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    let num_matches = matches.len();
    let matches_graph = CompatGraph::new(matches);

    let mut subgraph = BitSet::with_capacity(num_matches);
    for i in 0..num_matches {
        subgraph.insert(i);
    }

    // Kernelization
    if kernel_method != Kernel::Never {
        subgraph = deletion_kernel(&matches_graph, subgraph);
    }

    // Search
    let mut total_search = 0;
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
    let edge_count = mol.graph().edge_count();

    let index = recurse_clique_index_search(
        mol, 
        &[init], 
        edge_count - 1, 
        edge_count - 1,
        bounds,
        &mut total_search,
        subgraph.clone(),
        &matches_graph,
        1,
    &vec![],
        &kernel_method,);

    (index as u32, num_matches as u32, total_search)
}

pub fn clique_index_search_bench(mol: &Molecule, matches: Vec<(BitSet, BitSet)>, kernel_method: Kernel) -> (u32, u32, usize) {
    // Graph Initialization
    let num_matches = matches.len();
    let matches_graph = CompatGraph::new(matches);

    let mut subgraph = BitSet::with_capacity(num_matches);
    for i in 0..num_matches {
        subgraph.insert(i);
    }

    // Kernelization
    if kernel_method != Kernel::Never {
        subgraph = deletion_kernel(&matches_graph, subgraph);
    }

    // Search
    let mut total_search = 0;
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));
    let edge_count = mol.graph().edge_count();
    let bounds = vec![
        Bound::IntChain,
        Bound::VecChainSimple,
        Bound::VecChainSmallFrags,
    ];

    let index = recurse_clique_index_search(
        mol, 
        &[init], 
        edge_count - 1, 
        edge_count - 1,
        &bounds,
        &mut total_search,
        subgraph.clone(),
        &matches_graph,
        1,
        &vec![],
        &kernel_method,);

    (index as u32, num_matches as u32, total_search)
}

fn deletion_kernel(g: &CompatGraph, mut subgraph: BitSet) -> BitSet {
    let subgraph_copy = subgraph.clone();

    for v in subgraph_copy.iter() {
        let v_val = g.weights[v];
        let neighbors_v = g.neighbors(v, &subgraph);

        let Some(w1) = neighbors_v.iter().next() else {
            continue;
        };
        let Some(w2) = neighbors_v.iter().last() else {
            continue;
        };

        let mut s = subgraph.clone();
        s.intersect_with(&g.graph[w1]);
        for u in s.intersection(&&g.graph[w2]) {
            if g.are_adjacent(v, u) || v == u {
                continue;
            }

            let u_val = g.weights[u];
            if v_val > u_val {
                continue;
            }

            let neighbors_u = g.neighbors(u, &subgraph);

            if neighbors_v.is_subset(&neighbors_u) {
                subgraph.remove(v);
                break;
            }
        }
    }

    subgraph
}

fn inclusion_kernel(g: &CompatGraph, subgraph: &BitSet) -> Vec<usize> {
    let mut kernel = Vec::new();
    let tot = subgraph.iter().map(|v| g.weights[v]).sum::<usize>();

    'outer: for v in subgraph {
        let vw = g.weights[v];
        let nw = g.neighbors(v, subgraph).iter().map(|u| g.weights[u]).sum::<usize>();
        if vw >= tot - nw - vw { 
            kernel.push(v);
            continue;
        }

        let mut neighbors: Vec<usize> = vec![];

        for u in subgraph.difference(&g.graph[v]) {
            if u == v {
                continue;
            }
            if g.weights[u] > vw {
                continue 'outer;
            }

            for w in neighbors.iter() {
                if g.are_adjacent(u, *w) {
                    continue 'outer;
                }   
            }

            neighbors.push(u);
        }

        kernel.push(v);
    }

    kernel
}


#[allow(clippy::too_many_arguments)]
fn parallel_recurse_index_search(
    mol: &Molecule,
    matches: &[(BitSet, BitSet)],
    fragments: &[BitSet],
    ix: usize,
    largest_remove: usize,
    best: AtomicUsize,
    bounds: &[Bound],
    states_searched: Arc<AtomicUsize>,
) -> usize {
    let cx = AtomicUsize::from(ix);

    states_searched.fetch_add(1, Relaxed);

    // Branch and Bound
    for bound_type in bounds {
        let best = best.load(Relaxed);
        let exceeds = match bound_type {
            Bound::Log => ix - log_bound(fragments) >= best,
            Bound::IntChain => ix - addition_bound(fragments, largest_remove) >= best,
            Bound::VecChainSimple => ix - vec_bound_simple(fragments, largest_remove, mol) >= best,
            Bound::VecChainSmallFrags => {
                ix - vec_bound_small_frags(fragments, largest_remove, mol) >= best
            },
            _ => false
        };
        if exceeds {
            return ix;
        }
    }

    // Search for duplicatable fragment
    matches.par_iter().enumerate().for_each(|(i, (h1, h2))| {
        let mut fractures = fragments.to_owned();
        let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
        let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

        let largest_remove = h1.len();

        let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) else {
            return;
        };

        // All of these clones are on bitsets and cheap enough
        if i1 == i2 {
            let mut union = h1.clone();
            union.union_with(h2);
            let mut difference = f1.clone();
            difference.difference_with(&union);
            let c = connected_components_under_edges(mol.graph(), &difference);
            fractures.extend(c);
            fractures.swap_remove(i1);
        } else {
            let mut f1r = f1.clone();
            f1r.difference_with(h1);
            let mut f2r = f2.clone();
            f2r.difference_with(h2);

            let c1 = connected_components_under_edges(mol.graph(), &f1r);
            let c2 = connected_components_under_edges(mol.graph(), &f2r);

            fractures.extend(c1);
            fractures.extend(c2);

            fractures.swap_remove(i1.max(i2));
            fractures.swap_remove(i1.min(i2));
        }

        fractures.retain(|i| i.len() > 1);
        fractures.push(h1.clone());

        let output = parallel_recurse_index_search(
            mol,
            &matches[i + 1..],
            &fractures,
            ix - h1.len() + 1,
            largest_remove,
            best.load(Relaxed).into(),
            bounds,
            states_searched.clone(),
        );
        cx.fetch_min(output, Relaxed);

        best.fetch_min(cx.load(Relaxed), Relaxed);
    });

    cx.load(Relaxed)
}

/// Computes information related to the assembly index of a molecule using the provided bounds.
///
/// The first result in the returned tuple is the assembly index of the molecule. The second result
/// gives the number of duplicatable subgraphs (pairs of disjoint and isomorphic subgraphs) in the
/// molecule. The third result is the number of states searched where a new state is considered to
/// be searched each time a duplicatable subgraph is removed.
///
/// If the search space of the molecule is large (>100) parallelization will be used.
///
/// Bounds will be used in the order provided in the `bounds` slice. Execution along a search path
/// will halt immediately after finding a bound that exceeds the current best assembly pathway. It
/// is generally better to provide bounds that are quick to compute first.
///
/// # Example
/// ```
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// use assembly_theory::assembly::{Bound, index_search};
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene naively, with no bounds.
/// let (slow_index, _, _) = index_search(&benzene, &[]);
///
/// // Compute assembly index of benzene with the log and integer chain bounds
/// let (fast_index, _, _) = index_search(&benzene, &[Bound::Log, Bound::IntChain]);
///
/// assert_eq!(slow_index, 3);
/// assert_eq!(fast_index, 3);
/// # Ok(())
/// # }
/// ```
pub fn index_search(mol: &Molecule, bounds: &[Bound]) -> (u32, u32, usize) {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Create and sort matches array
    let mut matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

    let edge_count = mol.graph().edge_count();

    let (index, total_search) = if matches.len() > PARALLEL_MATCH_SIZE_THRESHOLD {
        let total_search = Arc::new(AtomicUsize::from(0));
        let index = parallel_recurse_index_search(
            mol,
            &matches,
            &[init],
            edge_count - 1,
            edge_count,
            (edge_count - 1).into(),
            bounds,
            total_search.clone(),
        );
        let total_search = total_search.load(Relaxed);
        (index as u32, total_search)
    } else {
        let mut total_search = 0;
        let index = recurse_index_search(
            mol,
            &matches,
            &[init],
            edge_count - 1,
            edge_count,
            edge_count - 1,
            bounds,
            &mut total_search,
            -1
        );
        (index as u32, total_search)
    };

    (index, matches.len() as u32, total_search)
}

/// Like [`index_search`], but no parallelism is used.
///
/// # Example
/// ```
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// use assembly_theory::assembly::{Bound, serial_index_search};
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene naively, with no bounds.
/// let (slow_index, _, _) = serial_index_search(&benzene, &[]);
///
/// // Compute assembly index of benzene with the log and integer chain bounds
/// let (fast_index, _, _) = serial_index_search(&benzene, &[Bound::Log, Bound::IntChain]);
///
/// assert_eq!(slow_index, 3);
/// assert_eq!(fast_index, 3);
/// # Ok(())
/// # }
/// ```
pub fn serial_index_search(mol: &Molecule, bounds: &[Bound]) -> (u32, u32, usize) {
    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Create and sort matches array
    let mut matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    matches.sort_by(|e1, e2| e2.0.len().cmp(&e1.0.len()));

    let edge_count = mol.graph().edge_count();
    let mut total_search = 0;

    let index = recurse_index_search(
        mol,
        &matches,
        &[init],
        edge_count - 1,
        edge_count,
        edge_count - 1,
        bounds,
        &mut total_search,
        -1
    );

    (index as u32, matches.len() as u32, total_search)
}

fn log_bound(fragments: &[BitSet]) -> usize {
    let mut size = 0;
    for f in fragments {
        size += f.len();
    }

    size - (size as f32).log2().ceil() as usize
}

fn addition_bound(fragments: &[BitSet], m: usize) -> usize {
    let mut max_s: usize = 0;
    let mut frag_sizes: Vec<usize> = Vec::new();

    for f in fragments {
        frag_sizes.push(f.len());
    }

    let size_sum: usize = frag_sizes.iter().sum();

    // Test for all sizes m of largest removed duplicate
    for max in 2..m + 1 {
        let log = {
            if max <= 10 {
                ADD_CHAIN[max - 1]
            }
            else {
                (max as f32).log2().ceil() as usize
            }
        };
        let mut aux_sum: usize = 0;

        for len in &frag_sizes {
            aux_sum += (len / max) + (len % max != 0) as usize
        }

        max_s = max_s.max(size_sum - log - aux_sum);
    }

    max_s
}

// Count number of unique edges in a fragment
// Helper function for vector bounds
fn unique_edges(fragment: &BitSet, mol: &Molecule) -> Vec<EdgeType> {
    let g = mol.graph();
    let mut nodes: Vec<Element> = Vec::new();
    for v in g.node_weights() {
        nodes.push(v.element());
    }
    let edges: Vec<petgraph::prelude::EdgeIndex> = g.edge_indices().collect();
    let weights: Vec<Bond> = g.edge_weights().copied().collect();

    // types will hold an element for every unique edge type in fragment
    let mut types: Vec<EdgeType> = Vec::new();
    for idx in fragment.iter() {
        let bond = weights[idx];
        let e = edges[idx];

        let (e1, e2) = g.edge_endpoints(e).expect("bad");
        let e1 = nodes[e1.index()];
        let e2 = nodes[e2.index()];
        let ends = if e1 < e2 { (e1, e2) } else { (e2, e1) };

        let edge_type = EdgeType { bond, ends };

        if types.contains(&edge_type) {
            continue;
        } else {
            types.push(edge_type);
        }
    }

    types
}

fn vec_bound_simple(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
    // Calculate s (total number of edges)
    // Calculate z (number of unique edges)
    let mut s = 0;
    for f in fragments {
        s += f.len();
    }

    let mut union_set = BitSet::new();
    for f in fragments {
        union_set.union_with(f);
    }
    let z = unique_edges(&union_set, mol).len();

    (s - z) - ((s - z) as f32 / m as f32).ceil() as usize
}

fn vec_bound_small_frags(fragments: &[BitSet], m: usize, mol: &Molecule) -> usize {
    let mut size_two_fragments: Vec<BitSet> = Vec::new();
    let mut large_fragments: Vec<BitSet> = fragments.to_owned();
    let mut indices_to_remove: Vec<usize> = Vec::new();

    // Find and remove fragments of size 2
    for (i, frag) in fragments.iter().enumerate() {
        if frag.len() == 2 {
            indices_to_remove.push(i);
        }
    }
    for &index in indices_to_remove.iter().rev() {
        let removed_bitset = large_fragments.remove(index);
        size_two_fragments.push(removed_bitset);
    }

    // Compute z = num unique edges of large_fragments NOT also in size_two_fragments
    let mut fragments_union = BitSet::new();
    let mut size_two_fragments_union = BitSet::new();
    for f in fragments {
        fragments_union.union_with(f);
    }
    for f in size_two_fragments.iter() {
        size_two_fragments_union.union_with(f);
    }
    let z = unique_edges(&fragments_union, mol).len()
        - unique_edges(&size_two_fragments_union, mol).len();

    // Compute s = total number of edges in fragments
    // Compute sl = total number of edges in large fragments
    let mut s = 0;
    let mut sl = 0;
    for f in fragments {
        s += f.len();
    }
    for f in large_fragments {
        sl += f.len();
    }

    // Find number of unique size two fragments
    let mut size_two_types: Vec<(EdgeType, EdgeType)> = Vec::new();
    for f in size_two_fragments.iter() {
        let mut types = unique_edges(f, mol);
        types.sort();
        if types.len() == 1 {
            size_two_types.push((types[0], types[0]));
        } else {
            size_two_types.push((types[0], types[1]));
        }
    }
    size_two_types.sort();
    size_two_types.dedup();

    s - (z + size_two_types.len() + size_two_fragments.len())
        - ((sl - z) as f32 / m as f32).ceil() as usize
}

/// Computes the assembly index of a molecule using an effecient bounding strategy
/// # Example
/// ```
/// # use std::fs;
/// # use std::path::PathBuf;
/// # use assembly_theory::*;
/// # fn main() -> Result<(), std::io::Error> {
/// # let path = PathBuf::from(format!("./data/checks/benzene.mol"));
/// // Read a molecule data file
/// let molfile = fs::read_to_string(path).expect("Cannot read input file.");
/// let benzene = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
///
/// // Compute assembly index of benzene
/// assert_eq!(assembly::index(&benzene), 3);
/// # Ok(())
/// # }
/// ```
pub fn index(m: &Molecule) -> u32 {
    index_search(
        m,
        &[
            Bound::IntChain,
            Bound::VecChainSimple,
            Bound::VecChainSmallFrags,
        ],
    )
    .0
}
