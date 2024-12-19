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
