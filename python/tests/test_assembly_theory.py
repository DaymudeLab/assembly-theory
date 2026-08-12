"""
test_python: Test all public functions in the assembly_theory Python package.
"""

import os.path as osp
import pytest

import assembly_theory as at

# These tests use the following molecules:
# anthracene: https://www.kegg.jp/entry/C14315
# benzene:    https://www.kegg.jp/entry/C01407


def test_mol_info():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    info = at.mol_info(mol_block)
    num_atoms = info.count('label = "Atom')
    num_single = info.count('label = "Single"')
    num_double = info.count('label = "Double"')
    num_triple = info.count('label = "Triple"')

    assert (num_atoms, num_single, num_double, num_triple) == (14, 9, 7, 0)


def test_mol_info_bad_molblock():
    with pytest.raises(OSError) as e:
        at.mol_info("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_depth():
    with open(osp.join("data", "checks", "benzene.mol")) as f:
        mol_block = f.read()

    assert at.depth(mol_block) == 3


def test_depth_bad_molblock():
    with pytest.raises(OSError) as e:
        at.depth("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    assert at.index(mol_block) == 6


def test_index_bad_molblock():
    with pytest.raises(OSError) as e:
        at.index("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index_search():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    (index, num_matches, states_searched, pathways) = at.index_search(
        mol_block,
        timeout=None,
        canonize_str="tree-nauty",
        parallel_str="none",  # Make states_searched deterministic.
        memoize_str="none",
        kernel_str="none",
        bound_strs=["int", "matchable-edges"],
        max_pathways=1,
    )

    pathway_str = """digraph {
    0 [ label = "{14}" ]
    1 [ label = "{15}" ]
    2 [ label = "{14, 15}" ]
    3 [ label = "{7, 14, 15}" ]
    4 [ label = "{6, 7, 14, 15}" ]
    5 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
    6 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
    7 [ label = "{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}" ]
    0 -> 2 [ label = "{14}" ]
    1 -> 2 [ label = "{15}" ]
    0 -> 3 [ label = "{7}" ]
    2 -> 3 [ label = "{14, 15}" ]
    1 -> 4 [ label = "{6}" ]
    3 -> 4 [ label = "{7, 14, 15}" ]
    4 -> 5 [ label = "{6, 7, 14, 15}" ]
    3 -> 5 [ label = "{3, 4, 5}" ]
    2 -> 6 [ label = "{2, 13}" ]
    5 -> 6 [ label = "{3, 4, 5, 6, 7, 14, 15}" ]
    6 -> 7 [ label = "{2, 3, 4, 5, 6, 7, 13, 14, 15}" ]
    5 -> 7 [ label = "{0, 1, 8, 9, 10, 11, 12}" ]
}
"""

    assert (index, num_matches, states_searched, pathways) == (
        6,
        466,
        491,
        [pathway_str],
    )


def test_index_search_timeout():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    (_, _, states_searched, pathways) = at.index_search(
        mol_block,
        timeout=1,  # Limit search to 1 ms, which should time out.
        canonize_str="tree-nauty",
        parallel_str="none",
        memoize_str="none",
        kernel_str="none",
        bound_strs=[],
        max_pathways=0,
    )

    assert (states_searched is None) and (pathways == [])


def test_index_search_bad_molblock():
    with pytest.raises(OSError) as e:
        at.index_search("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index_search_bad_canonize():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, canonize_str="invalid-mode")

    assert e.type is ValueError and "Invalid canonization" in str(e.value)


def test_index_search_bad_parallel():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, parallel_str="invalid-mode")

    assert e.type is ValueError and "Invalid parallelization" in str(e.value)


def test_index_search_bad_memoize():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, memoize_str="invalid-mode")

    assert e.type is ValueError and "Invalid memoization" in str(e.value)


def test_index_search_bad_kernel():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, kernel_str="invalid-mode")

    assert e.type is ValueError and "Invalid kernelization" in str(e.value)


def test_index_search_bad_bound():
    with open(osp.join("data", "checks", "anthracene.mol")) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, bound_strs=["int", "invalid-bound"])

    assert e.type is ValueError and "Invalid bound" in str(e.value)
