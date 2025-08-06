use std::os::raw::{c_int, c_void};

use nauty_Traces_sys::{
    densenauty, empty_graph, groupautomproc, grouplevelproc, groupptr, grouprec, makecosetreps,
    optionblk, statsblk, FALSE,
};
use petgraph::{EdgeType, Graph};

use crate::nauty::{canon::bit_adj_to_graph, dense::DenseGraph, dense::Nodes};

#[repr(C)]
#[derive(Debug)]
pub struct AutoGroups<N> {
    /// The automorphism group of the graph.
    pub data: Vec<Vec<i32>>,

    /// The number of automorphisms of the graph.
    pub count: usize,

    /// The number of nodes in the graph.
    pub n: usize,

    /// The canonical labelling of the graph.
    pub canon: Vec<u64>,

    /// The `Nodes` of the graph.
    pub nodes: Nodes<N>,
}
impl<N> AutoGroups<N>
where
    N: Ord + Eq + Clone,
{
    pub fn new(n: usize, canon: Vec<u64>, nodes: Nodes<N>) -> Self {
        Self {
            data: Vec::with_capacity(n * 100),
            count: 0,
            n,
            canon,
            nodes,
        }
    }

    pub fn from_petgraph<Ty: EdgeType>(graph: &Graph<N, (), Ty>) -> Self {
        let dense = DenseGraph::<N>::from_petgraph(graph);
        Self::from_dense(dense, graph.is_directed())
    }

    pub fn from_dense(mut dense: DenseGraph<N>, is_directed: bool) -> Self {
        let mut stats = statsblk::default();
        let mut options = autom_opts(is_directed);
        let mut canon = empty_graph(dense.m, dense.n);

        unsafe {
            densenauty(
                dense.g.as_mut_ptr(),
                dense.nodes.lab.as_mut_ptr(),
                dense.nodes.ptn.as_mut_ptr(),
                dense.nodes.orbits.as_mut_ptr(),
                &mut options,
                &mut stats,
                dense.m as c_int,
                dense.n as c_int,
                canon.as_mut_ptr(),
            );
            let mut autogroups = AutoGroups::new(dense.n, canon, dense.nodes);
            let group = groupptr(FALSE);
            makecosetreps(group);

            allgroup3(
                group,
                Some(writeautom3::<N>),
                &mut autogroups as *mut AutoGroups<N> as *mut c_void,
            );
            autogroups
        }
    }

    pub fn orbits(&self) -> &[i32] {
        &self.nodes.orbits
    }

    pub fn automorphisms(&self) -> &Vec<Vec<i32>> {
        &self.data
    }

    pub fn size(&self) -> usize {
        self.count
    }

    pub fn n_nodes(&self) -> usize {
        self.n
    }

    pub fn canonical(&self) -> &Vec<u64> {
        &self.canon
    }
}

fn autom_opts(is_directed: bool) -> optionblk {
    optionblk {
        writeautoms: 0,
        getcanon: 1,
        digraph: is_directed.into(),
        userautomproc: Some(groupautomproc),
        userlevelproc: Some(grouplevelproc),
        ..optionblk::default()
    }
}

impl<N, Ty> From<AutoGroups<N>> for Graph<(), (), Ty>
where
    Ty: EdgeType,
{
    fn from(ag: AutoGroups<N>) -> Self {
        bit_adj_to_graph(&ag.canon, ag.n * ag.n, ag.n)
    }
}

impl<N, Ty> From<&AutoGroups<N>> for Graph<(), (), Ty>
where
    Ty: EdgeType,
{
    fn from(ag: &AutoGroups<N>) -> Self {
        bit_adj_to_graph(&ag.canon, ag.n * ag.n, ag.n)
    }
}

// bindgen gets this wrong somehow? linker won't find it.
extern "C" {
    pub fn allgroup3(
        arg1: *mut grouprec,
        arg2: ::std::option::Option<
            unsafe extern "C" fn(
                arg1: *mut ::std::os::raw::c_int,
                arg2: ::std::os::raw::c_int,
                arg3: *mut ::std::os::raw::c_int,
                arg4: *mut ::std::os::raw::c_void,
            ),
        >,
        arg3: *mut ::std::os::raw::c_void,
    );
}

extern "C" fn writeautom3<N>(p: *mut i32, n: i32, _abort: *mut i32, userptr: *mut c_void) {
    unsafe {
        // convert void pointer to an AutoGroups pointer
        let autogroups_ptr = userptr as *mut AutoGroups<N>;

        // convert the raw pointer to a mutable reference
        let autogroups = &mut *autogroups_ptr;

        // push the node id to the vector
        autogroups.data.push(Vec::with_capacity(n as usize));
        for i in 0..n {
            autogroups.data[autogroups.count].push(*p.offset(i as isize));
        }
        autogroups.count += 1;
    }
}

#[cfg(test)]
mod testing {

    use super::*;
    use petgraph::{Directed, Undirected};

    #[test]
    fn autogroups_directed() {
        let e1 = vec![(0, 1), (1, 2), (2, 0)];
        let e2 = vec![(0, 1), (1, 2), (2, 0), (1, 0), (2, 1), (0, 2)];
        let g1 = Graph::<(), (), Directed>::from_edges(&e1);
        let g2 = Graph::<(), (), Directed>::from_edges(&e2);
        let a1 = AutoGroups::from_petgraph(&g1);
        let a2 = AutoGroups::from_petgraph(&g2);
        let expected_1 = vec![vec![0, 1, 2], vec![1, 2, 0], vec![2, 0, 1]];
        let expected_2 = vec![
            vec![0, 1, 2],
            vec![1, 0, 2],
            vec![2, 0, 1],
            vec![0, 2, 1],
            vec![1, 2, 0],
            vec![2, 1, 0],
        ];

        assert_eq!(a1.size(), 3);
        assert_eq!(a2.size(), 6);
        assert_eq!(a1.n_nodes(), 3);
        assert_eq!(a2.n_nodes(), 3);

        for e in expected_1 {
            assert!(a1.automorphisms().contains(&e));
        }
        for e in expected_2 {
            assert!(a2.automorphisms().contains(&e));
        }
    }

    #[test]
    fn autogroups_undirected() {
        let e1 = vec![(0, 1), (1, 2), (2, 0)];
        let e2 = vec![(0, 1), (1, 2), (2, 0), (1, 0), (2, 1), (0, 2)];
        let g1 = Graph::<(), (), Undirected>::from_edges(&e1);
        let g2 = Graph::<(), (), Undirected>::from_edges(&e2);
        let a1 = AutoGroups::from_petgraph(&g1);
        let a2 = AutoGroups::from_petgraph(&g2);
        let expected = vec![
            vec![0, 1, 2],
            vec![1, 0, 2],
            vec![2, 0, 1],
            vec![0, 2, 1],
            vec![1, 2, 0],
            vec![2, 1, 0],
        ];

        assert_eq!(a1.size(), 6);
        assert_eq!(a2.size(), 6);
        assert_eq!(a1.n_nodes(), 3);
        assert_eq!(a2.n_nodes(), 3);

        for e in expected {
            assert!(a1.automorphisms().contains(&e));
            assert!(a2.automorphisms().contains(&e));
        }
    }
}
