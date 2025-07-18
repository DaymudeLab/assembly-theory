//! Kernelize molecular graphs to improve top-down search efficiency.
//!
//! TODO: Longer explanation of what that means.

use clap::ValueEnum;

/// Graph kernelization strategy when searching using the clique reduction.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum KernelMode {
    /// No kernelization.
    None,
    /// Only kernelize the original molecule.
    Once,
    /// Kernelize the original molecule and the recursion's first level only.
    DepthOne,
    /// Perform kernelization at every recursive step.
    Always,
}
