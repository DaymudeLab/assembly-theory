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
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    info = at.mol_info(mol_block)
    num_atoms = info.count("label = \"Atom")
    num_single = info.count("label = \"Single\"")
    num_double = info.count("label = \"Double\"")
    num_triple = info.count("label = \"Triple\"")

    assert (num_atoms, num_single, num_double, num_triple) == (14, 9, 7, 0)


def test_mol_info_bad_molblock():
    with pytest.raises(OSError) as e:
        at.mol_info("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_depth():
    with open(osp.join('data', 'checks', 'benzene.mol')) as f:
        mol_block = f.read()

    assert at.depth(mol_block) == 3


def test_depth_bad_molblock():
    with pytest.raises(OSError) as e:
        at.depth("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    assert at.index(mol_block) == 6


def test_index_bad_molblock():
    with pytest.raises(OSError) as e:
        at.index("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index_search():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    (index, num_matches, states_searched) = at.index_search(
            mol_block,
            "grow-erode",
            "tree-nauty",
            "none",  # Disable parallelism for deterministic states_searched.
            "none",
            "none",
            ["int", "vec-simple", "vec-small-frags"])

    assert (index, num_matches, states_searched) == (6, 466, 2562)


def test_index_search_bad_molblock():
    with pytest.raises(OSError) as e:
        at.index_search("This string is not the contents of a .mol file.")

    assert e.type is OSError


def test_index_search_bad_enumerate():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, enumerate_str="invalid-mode")

    assert e.type is ValueError and "Invalid enumeration" in str(e.value)


def test_index_search_bad_canonize():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, canonize_str="invalid-mode")

    assert e.type is ValueError and "Invalid canonization" in str(e.value)


def test_index_search_bad_parallel():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, parallel_str="invalid-mode")

    assert e.type is ValueError and "Invalid parallelization" in str(e.value)


def test_index_search_bad_memoize():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, memoize_str="invalid-mode")

    assert e.type is ValueError and "Invalid memoization" in str(e.value)


def test_index_search_bad_kernel():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, kernel_str="invalid-mode")

    assert e.type is ValueError and "Invalid kernelization" in str(e.value)


def test_index_search_bad_bound():
    with open(osp.join('data', 'checks', 'anthracene.mol')) as f:
        mol_block = f.read()

    with pytest.raises(ValueError) as e:
        at.index_search(mol_block, bound_strs=["int", "invalid-bound"])

    assert e.type is ValueError and "Invalid bound" in str(e.value)
