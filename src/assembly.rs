use std::collections::BTreeSet;

use bit_set::BitSet;

use crate::{molecule::Molecule, utils::connected_components_under_edges};

fn top_down_search(m: &Molecule) -> u32 {
    let mut ix = u32::MAX;
    for (left, right) in m.partitions().unwrap() {
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

fn remnant_search(m: &Molecule) -> u32 {
    fn recurse(
        m: &Molecule,
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
                if f1 == f2 {
                    let mut union = h1.clone();
                    union.union_with(h2);
                    let mut difference = f1.clone();
                    difference.difference_with(&union);
                    let c = connected_components_under_edges(m.graph(), &difference);
                    fractures.extend(c);
                    fractures.swap_remove(i1);
                    fractures.push(h1.clone());
                } else {
                    let mut f1r = f1.clone();
                    f1r.difference_with(h1);
                    let mut f2r = f2.clone();
                    f2r.difference_with(h2);

                    let c1 = connected_components_under_edges(m.graph(), &f1r);
                    let c2 = connected_components_under_edges(m.graph(), &f2r);

                    fractures.extend(c1);
                    fractures.extend(c2);

                    fractures.swap_remove(i1.max(i2));
                    fractures.swap_remove(i1.min(i2));

                    fractures.push(h1.clone());
                }
                cx = cx.min(recurse(
                    m,
                    matches,
                    &fractures,
                    ix - h1.len() + 1,
                ));
            }
        }
        cx
    }

    let mut init = BitSet::new();
    init.extend(m.graph().edge_indices().map(|ix| ix.index()));

    recurse(
        m,
        &m.matches().collect(),
        &[init],
        m.graph().edge_count() - 1,
    ) as u32
}

// Compute the assembly index of a molecule
pub fn index(m: &Molecule) -> u32 {
    remnant_search(m)
}

pub fn depth(m: &Molecule) -> u32 {
    top_down_search(m)
}

pub fn search_space(m: &Molecule) -> u32 {
    m.matches().count() as u32
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
    fn test_setup(filename: &str) {
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
            let molecule = loader::parse(&path).expect(&format!(
                "Cannot generate assembly index for molecule: {}.",
                name
            ));
            let index = index(&molecule);
            assert_eq!(index, *master_dataset.get(&name).unwrap());
        }
    }

    #[test]
    fn test_small() {
        test_setup("./tests/suite1.csv");
    }

    #[test]
    fn test_medium() {
        test_setup("./tests/suite2.csv");
    }

    #[test]
    fn test_large() {
        test_setup("./tests/suite3.csv");
    }
}
