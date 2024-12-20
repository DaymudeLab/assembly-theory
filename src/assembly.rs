use std::{collections::BTreeSet, u32};

use petgraph::graph::EdgeIndex;

use crate::{
    molecule::{isomorphic_subgraphs_of, Molecule},
    utils::{connected_components_under_edges, edges_contained_within},
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
<<<<<<< HEAD
    let mut t1:u32 = 0;
    let mut t2:u32 = 0;

    let mut matches = BTreeSet::new();
    for subgraph in m.enumerate_subgraphs() {
        let mut h = m.graph().clone();
        h.retain_nodes(|_, n| subgraph.contains(&n));

        let h_prime = m.graph().map(
            |_, n| *n,
            |i, e| {
                let (src, dst) = m.graph().edge_endpoints(i).unwrap();
                (!subgraph.contains(&src) || !subgraph.contains(&dst)).then_some(*e)
            },
        );

        for cert in isomorphic_subgraphs_of(&h, &h_prime) {
            let cert = BTreeSet::from_iter(cert);
            let cert = BTreeSet::from_iter(edges_contained_within(m.graph(), &cert));
            let comp = BTreeSet::from_iter(edges_contained_within(m.graph(), &subgraph));
            matches.insert(if cert < comp {
                (cert, comp)
            } else {
                (comp, cert)
            });

            t1 = t1 + 1;
            if(t1 == 1000000) {
                t1 = 0;
                t2 = t2 + 1;
                println!("{}\n", t2);
            }
        }
    }

=======
>>>>>>> d507ecf2d38e3f2687aa3f7bcb3748dc40ba9457
    fn recurse(
        m: &Molecule,
        matches: &BTreeSet<(BTreeSet<EdgeIndex>, BTreeSet<EdgeIndex>)>,
        fragments: &Vec<BTreeSet<EdgeIndex>>,
        ix: usize,
        depth: usize,
    ) -> usize {
        let mut cx = ix;
        for (h1, h2) in matches {
            let mut fractures = fragments.clone();
            let f1 = fragments.iter().enumerate().find(|(_, c)| h1.is_subset(c));
            let f2 = fragments.iter().enumerate().find(|(_, c)| h2.is_subset(c));

            if let (Some((i1, f1)), Some((i2, f2))) = (f1, f2) {
                if f1 == f2 {
                    let remainder = f1
                        .difference(&h1.union(h2).cloned().collect::<BTreeSet<EdgeIndex>>())
                        .cloned()
                        .collect::<BTreeSet<EdgeIndex>>();
                    let c = connected_components_under_edges(m.graph(), &remainder);
                    fractures.extend(c);
                    fractures.swap_remove(i1);
                    fractures.push(h1.clone());
                } else {
                    let f1r = f1.difference(h1).cloned().collect::<BTreeSet<EdgeIndex>>();
                    let f2r = f2.difference(h2).cloned().collect::<BTreeSet<EdgeIndex>>();

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
                    depth + 1,
                ));
            }
        }
        cx
    }


    println!("Done!\n");
    recurse(
        m,
        &m.matches().collect(),
        &vec![m.graph().edge_indices().collect()],
        m.graph().edge_count() - 1,
        0,
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

pub fn addition_chain_bound(m: usize, fragments: &Vec<BTreeSet<EdgeIndex>>) -> usize{
    let mut max_s: usize = 0;
    let mut frag_sizes: Vec<usize> = Vec::new();

    for f in fragments {
        frag_sizes.push(f.len());
    }

    let size_sum: usize = frag_sizes.iter().sum();

    for min_frag in 2..m {
        let log= usize::BITS - min_frag.leading_zeros() - 1;
        let mut aux_sum: usize = 0;

        for len in &frag_sizes {
            aux_sum += (len / min_frag) + (len % min_frag != 0) as usize
        }

        max_s = std::cmp::max(max_s, size_sum - log as usize - aux_sum);
        
    }

    return max_s;
}
