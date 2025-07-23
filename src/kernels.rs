//! Kernelize match-compatibility graphs to improve top-down search efficiency.
//!
//! The problem of computing the minimum assembly index of a molecule can be
//! reduced to finding the maximum weight clique in the compatibility graph
//! definined in ```reductions.rs```. Strucutral properties of this graph can
//! be used to determine vertices (each corresponding to a duplicatable subgraph
//! match) that will either never be used in an optimal solution or will be used
//! in some optimal solution. We call the process of identifying these vertices 
//! kernelization. Uses the strategies of neighborhood removal, isolated vertex 
//! removal, and domination as described in section 5.2 of
//! [Lamm et al. (2019)](https://doi.org/10.1137/1.9781611975499.12). (Note that
//! they solve the equivalent problem of weighted independent set).

use clap::ValueEnum;

/// Graph kernelization strategy when searching using the clique reduction.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum KernelMode {
    /// Do not apply any kernelizations.
    None,
    /// Apply kernels only after the initial construction of the compatibility
    /// graph.
    Once,
    /// Apply kernels after the initial construction of the compability 
    /// graph, and again after any fragmentations of the full molecule.
    DepthOne,
    /// Apply kernels after every fragmentation step.
    Always,
}
