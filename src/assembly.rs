use std::collections::BTreeSet;

use petgraph::graph::NodeIndex;

use crate::{
    molecule::{isomorphic_subgraphs_of, Molecule},
    utils::connected_components_under,
};

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
    let mut matches = Vec::new();
    for subgraph in m.enumerate_subgraphs() {
        let mut h = m.graph().clone();
        h.retain_nodes(|_, n| subgraph.contains(&n));
        let h_prime = m
            .graph()
            .map(|i, n| (!subgraph.contains(&i)).then_some(*n), |_, e| *e);

        for cert in isomorphic_subgraphs_of(&h, &h_prime) {
            let cert = BTreeSet::from_iter(cert);
            let comp = subgraph
                .clone()
                .into_iter()
                .collect::<BTreeSet<NodeIndex>>();
            matches.push((cert.clone(), comp.clone()));
            matches.push((comp, cert));
        }
    }
    fn recurse(
        m: &Molecule,
        matches: &Vec<(BTreeSet<NodeIndex>, BTreeSet<NodeIndex>)>,
        fragments: &Vec<BTreeSet<NodeIndex>>,
        mut ix: u32,
    ) -> u32 {
        for (h1, h2) in matches {
            let mut fractures = fragments.clone();
            let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
            let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

            if let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) {
                let f1r = f1.difference(h1);
                let f2r = f2.difference(h2);

                let c1 = connected_components_under(m.graph(), &f1r.cloned().collect());
                let c2 = connected_components_under(m.graph(), &f2r.cloned().collect());

                fractures.extend(c1);
                fractures.extend(c2);

                fractures.swap_remove(i1);
                fractures.swap_remove(i2);

                fractures.push(h1.clone())
            }

            ix = ix.min(recurse(m, matches, &fractures, ix - h1.len() as u32 + 1))
        }
        0
    }

    recurse(
        m,
        &matches,
        &vec![m.graph().node_indices().collect()],
        m.graph().edge_count() as u32,
    )
}

// Compute the assembly index of a molecule
pub fn index(m: &Molecule) -> u32 {
    remnant_search(m)
}

pub fn depth(m: &Molecule) -> u32 {
    top_down_search(m)
}


#[cfg(test)]
mod tests {
    use std::{collections::HashMap, path::PathBuf};

    use csv::ReaderBuilder;

    use crate::loader;

    use super::*;

    /*
        Read Master CSV
     */

    fn read_master() -> HashMap<String,u32> {
        let mut rdr = ReaderBuilder::new().from_path("./data/master.csv").expect("Master CSV should exist to run tests!");
        let mut master_records: HashMap<String,u32> = HashMap::new();
        for result in rdr.records() {
            let record = result.expect("Error while reading Master file content");
            let record_vec: Vec<&str> = record.iter().collect();
            master_records.insert(record_vec[0].to_string(), record_vec[1].to_string().parse::<u32>().expect("Error while reading Master file content: Assembly Index is not Integer"));
        }
        master_records
    }

    fn test_setup(filename: String) -> Vec<String> {
        let mut rdr = ReaderBuilder::new().from_path(filename).expect("Given Test CSV should exist to run tests!");
        let mut mol_names: Vec<String> = Vec::new();
        for result in rdr.records() {
            let record = result.expect("Error while reading Test file content");
            for field in &record {
                mol_names.push(field.to_string());
            }
        }
        mol_names
    }

    #[test]
    fn test_small() {
        let master_dataset: HashMap<String, u32> = read_master();
        let test_mol_names: Vec<String> = test_setup("./tests/suite1.csv".to_string());

        for mol in test_mol_names {
            let path = PathBuf::from(format!("./data/{}", mol));
            let molecule = loader::parse(&path).expect(&format!("Error while generating assembly index for molecule: {}", mol));
            let index = index(&molecule);
            assert_eq!(index, *master_dataset.get(&mol).unwrap());
        }
    }

    #[test]
    fn test_medium() {
        let master_dataset: HashMap<String, u32> = read_master();
        let test_mol_names: Vec<String> = test_setup("./tests/suite2.csv".to_string());

        for mol in test_mol_names {
            let path = PathBuf::from(format!("./data/{}", mol));
            let molecule = loader::parse(&path).expect(&format!("Error while generating assembly index for molecule: {}", mol));
            let index = index(&molecule);
            assert_eq!(index, *master_dataset.get(&mol).unwrap());
        }
    }

    #[test]
    fn test_large() {
        let master_dataset: HashMap<String, u32> = read_master();
        let test_mol_names: Vec<String> = test_setup("./tests/suite3.csv".to_string());

        for mol in test_mol_names {
            let path = PathBuf::from(format!("./data/{}", mol));
            let molecule = loader::parse(&path).expect(&format!("Error while generating assembly index for molecule: {}", mol));
            let index = index(&molecule);
            assert_eq!(index, *master_dataset.get(&mol).unwrap());
        }
    }
}