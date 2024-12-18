use petgraph::graph::NodeIndex;

use crate::molecule::{isomorphic_subgraphs_of, Molecule};

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
    let mut precompute = Vec::new();
    for subgraph in m.enumerate_subgraphs() {
        let mut h = m.graph().clone();
        h.retain_nodes(|_, n| subgraph.contains(&n));
        let h_prime = m
            .graph()
            .map(|i, n| (!subgraph.contains(&i)).then_some(*n), |_, e| *e);

        for cert in isomorphic_subgraphs_of(&h, &h_prime) {
            let comp = subgraph.clone().into_iter().collect::<Vec<NodeIndex>>();
            precompute.push((cert.clone(), comp.clone()));
            precompute.push((comp, cert));
        }
    }

    for e in precompute {
        println!("{:?}", e);
    }
    0
}

// Compute the assembly index of a molecule
pub fn index(m: &Molecule) -> u32 {
    remnant_search(m)
}

pub fn depth(m: &Molecule) -> u32 {
    top_down_search(m)
}
