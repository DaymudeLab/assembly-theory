//! Expose public `assembly_theory` functionality to a Python package using
//! [`pyo3`](https://docs.rs/pyo3/latest/pyo3/).
//!
//! The package is available [on PyPI](https://pypi.org/project/assembly-theory/);
//! see that README for installation and usage instructions. To build the
//! Python package directly from this crate's source code, see the instructions
//! in the [GitHub README](https://github.com/DaymudeLab/assembly-theory).
//!
//! Note that all Rust functions in this module have the form `_fn_name`, which
//! correspond to the actual Rust function `fn_name` elsewhere in the crate and
//! are exposed to the Python package as `fn_name`.
//!
//! # Python Example
//!
//! ```custom,{class=language-python}
//! import assembly_theory as at
//!
//! # Load a mol block from file...
//! with open('data/checks/anthracene.mol') as f:
//!     mol_block = f.read()
//!
//! # ...or define the mol block directly.
//! mol_block = """
//!
//!
//!  14 16  0  0  0  0  0  0  0  0999 V2000
//!    25.2202  -16.2366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    25.2202  -17.6385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    26.4373  -18.3394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    26.4373  -15.5356    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    27.6471  -16.2366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    27.6412  -17.6385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    28.8523  -18.3446    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    28.8644  -15.5409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    30.0755  -16.2469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    30.1327  -17.6453    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    31.2674  -18.3552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    32.4846  -17.6672    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    32.4973  -16.2688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!    31.2927  -15.5589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!   9  8  2  0     0  0
//!   8  5  1  0     0  0
//!   5  4  2  0     0  0
//!   4  1  1  0     0  0
//!   1  2  2  0     0  0
//!   2  3  1  0     0  0
//!   3  6  2  0     0  0
//!   9 10  1  0     0  0
//!  10 11  1  0     0  0
//!  11 12  2  0     0  0
//!  12 13  1  0     0  0
//!  13 14  2  0     0  0
//!  14  9  1  0     0  0
//!   5  6  1  0     0  0
//!   6  7  1  0     0  0
//!   7 10  2  0     0  0
//! M  END"""
//!
//! # Calculate the molecule's assembly index.
//! at.index(mol_block)  # 6
//! ```
//!
//! To customize assembly index calculation beyond the default strategy, set a
//! timeout for long-running calculations, and/or reconstruct minimum assembly
//! pathways, see [`_index_search`].

use std::str::FromStr;

use pyo3::{
    exceptions::{PyOSError, PyValueError},
    prelude::*,
    PyErr,
};

use crate::{
    assembly::{depth, index, index_search, ParallelMode},
    bounds::Bound as OurBound,
    canonize::CanonizeMode,
    kernels::KernelMode,
    loader::{parse_molfile_str, ParserError},
    memoize::MemoizeMode,
};

/// Implement a Python version of [`crate::loader::ParserError`].
impl From<ParserError> for PyErr {
    fn from(err: ParserError) -> PyErr {
        PyOSError::new_err(err.to_string())
    }
}

// TODO: Is there a clean way of avoiding the duplication of all our various
// algorithm variant enums?

/// Mirrors the [`CanonizeMode`] enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyCanonizeMode {
    Nauty,
    Faulon,
    TreeNauty,
    TreeFaulon,
}

/// Mirrors the [`ParallelMode`] enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyParallelMode {
    None,
    DepthOne,
    Always,
}

/// Mirrors the [`MemoizeMode`] enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyMemoizeMode {
    None,
    CanonIndex,
}

/// Mirrors the `kernels::KernelMode` enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyKernelMode {
    None,
    Once,
    DepthOne,
    Always,
}

/// Mirrors the `bounds::Bound` enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum PyBound {
    Log,
    Int,
    VecSimple,
    VecSmallFrags,
    MatchableEdges,
}

/// Converts bound options in `&str` format to `PyCanonizeMode`.
impl FromStr for PyCanonizeMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "nauty" => Ok(PyCanonizeMode::Nauty),
            "faulon" => Ok(PyCanonizeMode::Faulon),
            "tree-nauty" => Ok(PyCanonizeMode::TreeNauty),
            "tree-faulon" => Ok(PyCanonizeMode::TreeFaulon),
            _ => Err(PyValueError::new_err(format!(
                "Invalid canonization mode \"{s}\", options are: \
                [\"nauty\", \"faulon\", \"tree-nauty\", \"tree-faulon\"]"
            ))),
        }
    }
}

/// Converts bound options in `&str` format to `PyParallelMode`.
impl FromStr for PyParallelMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" => Ok(PyParallelMode::None),
            "depth-one" => Ok(PyParallelMode::DepthOne),
            "always" => Ok(PyParallelMode::Always),
            _ => Err(PyValueError::new_err(format!(
                "Invalid parallelization mode \"{s}\", options are: \
                [\"none\", \"depth-one\", \"always\"]"
            ))),
        }
    }
}

/// Converts bound options in `&str` format to `PyMemoizeMode`.
impl FromStr for PyMemoizeMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" => Ok(PyMemoizeMode::None),
            "canon-index" => Ok(PyMemoizeMode::CanonIndex),
            _ => Err(PyValueError::new_err(format!(
                "Invalid memoization mode \"{s}\", options are: \
                [\"none\", \"frags-index\", \"canon-index\"]"
            ))),
        }
    }
}

/// Converts bound options in `&str` format to `PyKernelMode`.
impl FromStr for PyKernelMode {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" => Ok(PyKernelMode::None),
            "once" => Ok(PyKernelMode::Once),
            "depth-one" => Ok(PyKernelMode::DepthOne),
            "always" => Ok(PyKernelMode::Always),
            _ => Err(PyValueError::new_err(format!(
                "Invalid kernelization mode \"{s}\", options are: \
                [\"none\", \"once\", \"depth-one\", \"always\"]"
            ))),
        }
    }
}

/// Converts bound options in `&str` format to `PyBound`.
impl FromStr for PyBound {
    type Err = PyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "log" => Ok(PyBound::Log),
            "int" => Ok(PyBound::Int),
            "vec-simple" => Ok(PyBound::VecSimple),
            "vec-small-frags" => Ok(PyBound::VecSmallFrags),
            "matchable-edges" => Ok(PyBound::MatchableEdges),
            _ => Err(PyValueError::new_err(format!(
                "Invalid bound \"{s}\", options are: \
                [\"log\", \"int\", \"vec-simple\", \"vec-small-frags\", \
                \"matchable-edges\"]"
            ))),
        }
    }
}

/// Converts a `Vec<String>` of bound Python strings into a `Vec<PyBound>`,
/// raising an error if any bound string is invalid.
fn process_bound_strs(bound_strs: Vec<String>) -> PyResult<Vec<PyBound>> {
    bound_strs
        .iter()
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
}

/// Converts a slice of `PyBound`s into a vector of `bounds::Bound`s.
fn make_boundlist(pybounds: &[PyBound]) -> Vec<OurBound> {
    let mut boundlist = pybounds
        .iter()
        .flat_map(|b| match b {
            PyBound::Log => vec![OurBound::Log],
            PyBound::Int => vec![OurBound::Int],
            PyBound::VecSimple => vec![OurBound::VecSimple],
            PyBound::VecSmallFrags => vec![OurBound::VecSmallFrags],
            PyBound::MatchableEdges => vec![OurBound::MatchableEdges],
        })
        .collect::<Vec<_>>();
    boundlist.dedup();
    boundlist
}

/// Return a [DOT-formatted](https://graphviz.org/doc/info/lang.html) string
/// representation of this molecule's graph structure.
///
/// Python version of [`crate::molecule::Molecule::info`].
///
/// # Python Parameters
/// - `mol_block`: The contents of a `.mol` file as a `str`.
///
/// # Python Returns
/// - A DOT-formatted `str` representing the molecule's atoms and bonds.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
///
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
///
/// # Print the molecule's graph structure.
/// print(at.mol_info(mol_block))
///
/// # graph {
/// #     0 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     1 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     2 [ label = "Atom { element: Carbon, capacity: 0 }" ]
/// #     ...
/// #     0 -- 1 [ label = "Double" ]
/// #     1 -- 2 [ label = "Single" ]
/// #     2 -- 5 [ label = "Double" ]
/// #     ...
/// # }
/// ```
#[pyfunction(name = "mol_info")]
pub fn _mol_info(mol_block: &str) -> PyResult<String> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol = parse_molfile_str(mol_block)?;

    // Return molecule info.
    Ok(mol.info())
}

/// Compute assembly depth; see
/// [Pagel et al. (2024)](https://arxiv.org/abs/2409.05993).
///
/// Python version of [`depth`].
///
/// # Python Parameters
/// - `mol_block`: The contents of a `.mol` file as a `str`.
///
/// # Python Returns
/// - The molecule's `int` assembly depth.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
///
/// # Load a mol block from file.
/// with open('data/checks/benzene.mol') as f:
///     mol_block = f.read()
///
/// # Calculate the molecule's assembly index.
/// at.depth(mol_block)  # 3
/// ```
#[pyfunction(name = "depth")]
pub fn _depth(mol_block: &str) -> PyResult<u32> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol = parse_molfile_str(mol_block)?;

    // Calculate assembly depth.
    Ok(depth(&mol))
}

/// Compute a molecule's assembly index using an efficient default strategy.
///
/// Python version of [`index`].
///
/// To customize assembly index calculation beyond the default strategy, set a
/// timeout for long-running calculations, and/or reconstruct minimum assembly
/// pathways, see [`_index_search`].
///
/// # Python Parameters
/// - `mol_block`: The contents of a `.mol` file as a `str`.
///
/// # Python Returns
/// - The molecule's `int` assembly index.
///
/// # Python Example
///
/// ```custom,{class=language-python}
/// import assembly_theory as at
///
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
///
/// # Calculate the molecule's assembly index.
/// at.index(mol_block)  # 6
/// ```
#[pyfunction(name = "index")]
pub fn _index(mol_block: &str) -> PyResult<u32> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol = parse_molfile_str(mol_block)?;

    // Calculate the assembly index.
    Ok(index(&mol))
}

/// Compute a molecule's assembly index and related information using a
/// top-down recursive algorithm, parameterized by the specified options.
///
/// Python version of [`index_search`].
///
/// # Python Parameters
///
/// - `mol_block`: The contents of a `.mol` file as a `str`.
/// - `timeout`: An `int` duration in milliseconds after which search is
///   stopped and the best upper bound on the assembly index found so far is
///   returned, or `None` (default) if search is run to completion.
/// - `canonize_str`: A canonization mode from [`"nauty"`, `"faulon"`,
///   `"tree-nauty"` (default), `"tree-faulon"`]. See [`CanonizeMode`] for
///   details.
/// - `parallel_str`: A parallelization mode from [`"none"`, `"depth-one"`
///   (default), `"always"`]. See [`ParallelMode`] for details.
/// - `memoize_str`: A memoization mode from [`none`, `canon-index` (default)].
///   See [`MemoizeMode`] for details.
/// - `kernel_str`: A kernelization mode from [`"none"` (default), `"once"`,
///   `"depth-one"`, `"always"`]. See [`KernelMode`] for details.
/// - `bound_strs`: A list of bounds containing zero or more of [`"log"`,
///   `"int"`, `"vec-simple"`, `"vec-small-frags"`, `"matchable-edges"`].
///   The default bounds are [`"int"`, `"matchable-edges"`]. See
///   [`crate::bounds::Bound`] for details.
/// - `max_pathways`: An `int` maximum number of minimum assembly pathways to
///   reconstruct, `0` for all such pathways, or `None` (default) to disable
///   pathway reconstruction.
///
/// # Python Returns
///
/// A 4-tuple containing:
/// - The molecule's `int` assembly index (or an upper bound if timed out).
/// - The molecule's `int` number of edge-disjoint isomorphic subgraph pairs.
/// - The `int` number of assembly states searched if search completes, or
///   `None` otherwise (i.e., if search timed out).
/// - A list of [DOT-formatted](https://graphviz.org/doc/info/lang.html) `str`
///   representations of the molecule's minimum assembly pathways. This list is
///   nonempty if and only if search did not time out and `max_pathways` is not
///   `None`.
///
/// # Python Examples
///
/// ## Customizing Assembly Index Search Options
///
/// The two examples below show different customizations of the assembly index
/// search options, enabling/disabling parallelism, memoization, and bounding.
/// ```custom,{class=language-python}
/// import assembly_theory as at
///
/// # Load a mol block from file.
/// with open('data/checks/anthracene.mol') as f:
///     mol_block = f.read()
///
/// # Calculate the molecule's assembly index without any search optimizations
/// # like parallelism, memoization, kernelization, or bounding.
/// (slow_index, num_matches, states_searched, _) = at.index_search(
///     mol_block,
///     timeout=None,
///     canonize_str="tree-nauty",
///     parallel_str="none",        # Disable parallelism.
///     memoize_str="none",         # Disable memoization.
///     kernel_str="none",          # Disable kernelization.
///     bound_strs=[],              # Disable bounding.
///     max_pathways=None)
/// print(f"Assembly Index:  {slow_index}")       # 6
/// print(f"Matches:         {num_matches}")      # 466
/// print(f"States Searched: {states_searched}")  # 133814
///
/// # Repeat the calculation with parallelism, memoization, and some bounds.
/// (fast_index, num_matches, states_searched, _) = at.index_search(
///     mol_block,
///     timeout=None,
///     canonize_str="tree-nauty",
///     parallel_str="depth-one",   # Enable parallelism.
///     memoize_str="canon-index",  # Enable memoization.
///     kernel_str="none",
///     bound_strs=["log", "int"],  # Enable some bounds.
///     max_pathways=None)
/// print(f"Assembly Index:  {fast_index}")       # 6, unchanged
/// print(f"Matches:         {num_matches}")      # 466, unchanged
/// print(f"States Searched: {states_searched}")  # 1607, e.g., but parallelism makes this nondeterministic
/// ```
///
/// ## Setting a Timeout for Long Calculations
///
/// Assembly index calculation can be slow on large molecules, even with search
/// optimizations enabled. If `timeout` is set, two outcomes are possible:
/// 1. Search completes before the timeout and everything behaves as usual.
/// 2. The timeout is reached before search completes, search is interrupted,
///    and the best upper bound on the assembly index found so far is returned.
///    The number of states searched and assembly pathways are not returned.
/// ```custom,{class=language-python}
/// # Limit search to 10 s, which is plenty for search to complete as usual,
/// # even with some search optimizations disabled.
/// (true_index, num_matches, states_searched, pathways) = at.index_search(
///     mol_block,
///     timeout=10000,              # Set a 10 s timeout.
///     canonize_str="tree-nauty",
///     parallel_str="none",        # Disable parallelism.
///     memoize_str="none",         # Disable memoization.
///     kernel_str="none",          # Disable kernelization.
///     bound_strs=["log", "int"],  # Enable some bounds.
///     max_pathways=1)             # Reconstruct one pathway.
/// print(f"Assembly Index:         {true_index}")       # 6
/// print(f"Matches:                {num_matches}")      # 466
/// print(f"States Searched:        {states_searched}")  # 2454
/// print(f"Pathways Reconstructed: {len(pathways)}")    # 1
///
/// # Limit search to 1 ms, which will time out without more optimizations.
/// (index_bound, num_matches, states_searched, pathways) = at.index_search(
///     mol_block,
///     timeout=1,                  # Set a 1 ms timeout.
///     canonize_str="tree-nauty",
///     parallel_str="none",        # Disable parallelism.
///     memoize_str="none",         # Disable memoization.
///     kernel_str="none",          # Disable kernelization.
///     bound_strs=["log", "int"],  # Enable some bounds.
///     max_pathways=1)             # Reconstruct one pathway.
/// print(f"Assembly Index:         {index_bound}")      # >= 6
/// print(f"Matches:                {num_matches}")      # 466, unaffected
/// print(f"States Searched:        {states_searched}")  # None, not computed on timeout
/// print(f"Pathways Reconstructed: {len(pathways)}")    # 0, not computed on timeout
/// ```
///
/// ## Reconstructing Minimum Assembly Pathways
///
/// An assembly pathway records the fragment join operations producing the
/// target molecule. Set `max_pathways` to `0` to reconstruct all minimum
/// assembly pathways discovered by the search process, `n` for at most the
/// first `n` such pathways, or `None` to disable pathway reconstruction.
/// See [`crate::pathway::Pathway::dag`] for details on the pathway data
/// structure.
///
/// Note that the definition of `max_pathway = 0`'s behavior is very carefully
/// worded: it finds all minimum assembly pathways *discovered by the search
/// process*. This is not to be taken as *all possible minimum assembly
/// pathways*. When search optimizations like parallelism, memoization, and
/// bounding are used, search prunes pathways that cannot yield a better
/// assembly index. This may include pathways that yield an *equal* (i.e., also
/// minimum) assembly index. Once pruned, these pathways will not survive to
/// pathway reconstruction.
/// ```custom,{class=language-python}
/// # Reconstruct a minimum assembly pathway discovered during search.
/// (_, _, _, pathways) = at.index_search(
///     mol_block,
///     timeout=None,
///     canonize_str="tree-nauty",
///     parallel_str="none",        # Disable parallelism for determinism.
///     memoize_str="canon-index",
///     kernel_str="none",
///     bound_strs=["log", "int"],
///     max_pathways=1)             # Reconstruct one pathway.
/// print(pathways[0])
///
/// # digraph {
/// #     0 [ label = "{14}" ]
/// #     1 [ label = "{15}" ]
/// #     2 [ label = "{14, 15}" ]
/// #     3 [ label = "{7, 14, 15}" ]
/// #     4 [ label = "{6, 7, 14, 15}" ]
/// #     5 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
/// #     6 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
/// #     7 [ label = "{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}" ]
/// #     0 -> 2 [ label = "{14}" ]
/// #     1 -> 2 [ label = "{15}" ]
/// #     0 -> 3 [ label = "{7}" ]
/// #     2 -> 3 [ label = "{14, 15}" ]
/// #     1 -> 4 [ label = "{6}" ]
/// #     3 -> 4 [ label = "{7, 14, 15}" ]
/// #     4 -> 5 [ label = "{6, 7, 14, 15}" ]
/// #     3 -> 5 [ label = "{3, 4, 5}" ]
/// #     2 -> 6 [ label = "{2, 13}" ]
/// #     5 -> 6 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
/// #     6 -> 7 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
/// #     5 -> 7 [ label = "{0, 1, 8, 9, 10, 11, 12}" ]
/// # }
/// #
/// ```
#[pyfunction(name = "index_search")]
#[pyo3(signature = (mol_block, timeout=None, canonize_str="tree-nauty", parallel_str="depth-one", memoize_str="canon-index", kernel_str="none", bound_strs=vec!["int".to_string(), "matchable-edges".to_string()], max_pathways=None), text_signature = "(mol_block, timeout=None, canonize_str=\"tree-nauty\", parallel_str=\"depth-one\", memoize_str=\"canon-index\", kernel_str=\"none\", bound_strs=[\"int\", \"matchable-edges\"], max_pathways=None)")]
pub fn _index_search(
    mol_block: &str,
    timeout: Option<u64>,
    canonize_str: &str,
    parallel_str: &str,
    memoize_str: &str,
    kernel_str: &str,
    bound_strs: Vec<String>,
    max_pathways: Option<usize>,
) -> PyResult<(u32, u32, Option<usize>, Vec<String>)> {
    // Parse the .mol file contents as a molecule::Molecule.
    let mol = parse_molfile_str(mol_block)?;

    // Parse the various modes and bound options.
    let canonize_mode = match PyCanonizeMode::from_str(canonize_str) {
        Ok(PyCanonizeMode::Nauty) => CanonizeMode::Nauty,
        Ok(PyCanonizeMode::Faulon) => CanonizeMode::Faulon,
        Ok(PyCanonizeMode::TreeNauty) => CanonizeMode::TreeNauty,
        Ok(PyCanonizeMode::TreeFaulon) => CanonizeMode::TreeFaulon,
        Err(e) => return Err(e),
    };
    let parallel_mode = match PyParallelMode::from_str(parallel_str) {
        Ok(PyParallelMode::None) => ParallelMode::None,
        Ok(PyParallelMode::DepthOne) => ParallelMode::DepthOne,
        Ok(PyParallelMode::Always) => ParallelMode::Always,
        Err(e) => return Err(e),
    };
    let memoize_mode = match PyMemoizeMode::from_str(memoize_str) {
        Ok(PyMemoizeMode::None) => MemoizeMode::None,
        Ok(PyMemoizeMode::CanonIndex) => MemoizeMode::CanonIndex,
        Err(e) => return Err(e),
    };
    let kernel_mode = match PyKernelMode::from_str(kernel_str) {
        Ok(PyKernelMode::None) => KernelMode::None,
        Ok(PyKernelMode::Once) => KernelMode::Once,
        Ok(PyKernelMode::DepthOne) => KernelMode::DepthOne,
        Ok(PyKernelMode::Always) => KernelMode::Always,
        Err(e) => return Err(e),
    };
    let pybounds = process_bound_strs(bound_strs)?;
    let boundlist = make_boundlist(&pybounds);

    // Compute assembly index.
    let (index, num_matches, states_searched, pathways) = index_search(
        &mol,
        timeout,
        canonize_mode,
        parallel_mode,
        memoize_mode,
        kernel_mode,
        &boundlist,
        max_pathways,
    );

    // Format assembly pathways as DOT strings.
    let mut pathways_dot = Vec::<String>::new();
    for pathway in pathways {
        pathways_dot.push(format!("{pathway}"));
    }

    Ok((index, num_matches, states_searched, pathways_dot))
}

/// A Python wrapper for the assembly_theory Rust crate.
// Registers the listed functions as a Python module named 'assembly_theory';
// the above line is used as a docstring.
#[pymodule(name = "assembly_theory")]
fn _assembly_theory(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_mol_info, m)?)?;
    m.add_function(wrap_pyfunction!(_depth, m)?)?;
    m.add_function(wrap_pyfunction!(_index, m)?)?;
    m.add_function(wrap_pyfunction!(_index_search, m)?)?;
    Ok(())
}
