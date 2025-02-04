use std::{cmp::Ordering, collections::BTreeSet};

use bit_set::BitSet;

use crate::{
    molecule::Bond, molecule::Element, molecule::Molecule, utils::connected_components_under_edges,
};

fn top_down_search(mol: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in mol.partitions().unwrap() {
        let l = if left.is_basic_unit() {
            0
        } else {
            top_down_search(&left)
        };

        let r = if right.is_basic_unit() {
            0
        } else {
            top_down_search(&right)
        };

        ix = ix.min(l.max(r) + 1)
    }
    ix
}

fn naive_search(mol: &Molecule) -> u32 {
    fn recurse(
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

            // All of these clones are on bitsets and cheap enough
            if let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) {
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
                cx = cx.min(recurse(mol, matches, &fractures, ix - h1.len() + 1));
            }
        }
        cx
    }

    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    recurse(
        mol,
        &mol.matches().collect(),
        &[init],
        mol.graph().edge_count() - 1,
    ) as u32
}

fn remnant_search(mol: &Molecule) -> (u32, u32) {
    fn recurse(
        mol: &Molecule,
        matches: &[(BitSet, BitSet)],
        fragments: &[BitSet],
        ix: usize,
        largest_remove: usize,
        mut best: usize,
        search_space: &mut u32,
    ) -> usize {
        let mut cx = ix;

        *search_space += 1;
        // Branch and Bound
        // Seet function

        let exceeds_add_chain_bound = ix - addition_chain_bound(largest_remove, fragments) >= best;
        let exceeds_vec_chain_bound = ix - vec_chain_bound(largest_remove, fragments, mol) >= best;
        let exceeds_vec_chain_bound2 =
            ix - vec_chain_bound2(largest_remove, fragments, mol) >= best;

        if exceeds_add_chain_bound || exceeds_vec_chain_bound || exceeds_vec_chain_bound2 {
            return ix;
        }

        for (i, (h1, h2)) in matches.iter().enumerate() {
            let mut fractures = fragments.to_owned();
            let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
            let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

            let largest_remove = h1.len();

            // All of these clones are on bitsets and cheap enough
            if let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) {
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

                cx = cx.min(recurse(
                    mol,
                    &matches[i + 1..],
                    &fractures,
                    ix - h1.len() + 1,
                    largest_remove,
                    best,
                    search_space,
                ));
                best = best.min(cx);
            }
        }

        cx
    }

    let mut init = BitSet::new();
    init.extend(mol.graph().edge_indices().map(|ix| ix.index()));

    // Create and sort matches array
    let mut matches: Vec<(BitSet, BitSet)> = mol.matches().collect();
    matches.sort_by(my_sort);

    fn my_sort(e1: &(BitSet, BitSet), e2: &(BitSet, BitSet)) -> Ordering {
        e2.0.len().cmp(&e1.0.len())
    }

    let mut total_search = 0;

    let ans = recurse(
        mol,
        &matches,
        &[init],
        mol.graph().edge_count() - 1,
        mol.graph().edge_count(),
        mol.graph().edge_count() - 1,
        &mut total_search,
    ) as u32;

    (ans, total_search)
}

// Compute the assembly index of a molecule
pub fn index(m: &Molecule) -> u32 {
    remnant_search(m).0
}

pub fn naive_index(m: &Molecule) -> u32 {
    naive_search(m)
}

pub fn depth(m: &Molecule) -> u32 {
    top_down_search(m)
}

pub fn search_space(m: &Molecule) -> u32 {
    m.matches().count() as u32
}

fn addition_chain_bound(m: usize, fragments: &[BitSet]) -> usize {
    let mut max_s: usize = 0;
    let mut frag_sizes: Vec<usize> = Vec::new();

    for f in fragments {
        frag_sizes.push(f.len());
    }

    let size_sum: usize = frag_sizes.iter().sum();

    for max in 2..m + 1 {
        let log = (max as f32).log2().ceil();
        let mut aux_sum: usize = 0;

        for len in &frag_sizes {
            aux_sum += (len / max) + (len % max != 0) as usize
        }

        max_s = std::cmp::max(max_s, size_sum - log as usize - aux_sum);
    }

    max_s
}

// Count number of unique edges in a fragment
fn unique_edges(fragment: &BitSet, mol: &Molecule) -> usize {
    let g = mol.graph();
    let mut nodes: Vec<Element> = Vec::new();
    for v in g.node_weights() {
        nodes.push(v.element);
    }
    let edges: Vec<petgraph::prelude::EdgeIndex> = g.edge_indices().collect();
    let weights: Vec<&Bond> = g.edge_weights().collect();

    let mut z: usize = 0;
    let mut types: Vec<(&Bond, (Element, Element))> = Vec::new();
    for idx in fragment.iter() {
        let bond = weights[idx];
        let e = edges[idx];

        let (e1, e2) = g.edge_endpoints(e).expect("bad");
        let e1 = nodes[e1.index()];
        let e2 = nodes[e2.index()];
        let ends = if e1 < e2 { (e1, e2) } else { (e2, e1) };

        let edge_type = (bond, ends);

        if types.iter().any(|&t| t == edge_type) {
            continue;
        } else {
            z += 1;
            types.push(edge_type);
        }
    }

    z
}

fn vec_chain_bound(m: usize, fragments: &[BitSet], mol: &Molecule) -> usize {
    let mut max_s: usize = 0;
    let mut frag_sizes: Vec<usize> = Vec::new();

    for f in fragments {
        frag_sizes.push(f.len());
    }

    let size_sum: usize = frag_sizes.iter().sum();

    let mut union_set = BitSet::new();
    for f in fragments {
        union_set.union_with(f);
    }
    let z = unique_edges(&union_set, mol);

    for max in 2..m + 1 {
        let mut aux_sum = 0;
        for f in fragments {
            let zi = unique_edges(f, mol);
            aux_sum += ((f.len() - zi) as f32 / max as f32).ceil() as usize;
        }

        let log = ((max as f32).log2() - (z as f32).log2()).ceil() as usize;
        max_s = max_s.max(size_sum + 1 - aux_sum - z - log);
    }

    max_s
}

// TODO: rename this function
fn vec_chain_bound2(m: usize, fragments: &[BitSet], mol: &Molecule) -> usize {
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
    let z = unique_edges(&union_set, mol);

    (s - z) - ((s - z) as f32 / m as f32).ceil() as usize
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, path::PathBuf};

    use csv::ReaderBuilder;

    use crate::loader;

    use super::*;

    // Read Master CSV
    fn read_master() -> HashMap<String, u32> {
        let mut reader = ReaderBuilder::new()
            .from_path("./data/master.csv")
            .expect("data/master.csv does not exist.");
        let mut master_records = HashMap::new();
        for result in reader.records() {
            let record = result.expect("master.csv is malformed.");
            let record = record.iter().collect::<Vec<_>>();
            master_records.insert(
                record[0].to_string(),
                record[1]
                    .to_string()
                    .parse::<u32>()
                    .expect("Assembly index is not an integer."),
            );
        }
        master_records
    }

    // Read Test CSV
    fn test_suite<F>(function: F, filename: &str)
    where
        F: Fn(&Molecule) -> u32,
    {
        let mut reader = ReaderBuilder::new()
            .from_path(filename)
            .expect("Test file does not exist.");
        let mut molecule_names: Vec<String> = Vec::new();
        for result in reader.records() {
            let record = result.expect("Cannot read test file.");
            for field in &record {
                molecule_names.push(field.to_string());
            }
        }
        let master_dataset: HashMap<String, u32> = read_master();
        for name in molecule_names {
            let path = PathBuf::from(format!("./data/{}", name));
            let molecule = loader::parse(&path).unwrap_or_else(|_| {
                panic!("Cannot generate assembly index for molecule: {}.", name)
            });
            let index = function(&molecule);
            assert_eq!(index, *master_dataset.get(&name).unwrap());
        }
    }

    #[test]
    fn bound_test_small() {
        test_suite(index, "./tests/suite1.csv");
    }

    #[test]
    fn bound_test_medium() {
        test_suite(index, "./tests/suite2.csv");
    }

    #[test]
    fn naive_test_small() {
        test_suite(naive_index, "./tests/suite1.csv");
    }

    #[test]
    fn naive_test_medium() {
        test_suite(naive_index, "./tests/suite2.csv");
    }
}
