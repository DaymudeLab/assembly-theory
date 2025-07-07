import assembly_theory as at
import rdkit.Chem


def test_aspirin_index():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly(mol) == 8


def test_aspirin_verbose():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly_verbose(mol) == {
        "duplicates": 20,
        "index": 8,
        "space": 35,
    }
    assert at.molecular_assembly_verbose(mol, no_bounds=True) == {
        "duplicates": 20,
        "index": 8,
        "space": 96,
    }

def test_anthracene_info():
    anthra_smi = "c1ccc2cc3ccccc3cc2c1"
    mol = rdkit.Chem.MolFromSmiles(anthra_smi)

    anthra_info = """graph {
    0 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    1 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    2 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    3 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    4 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    5 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    6 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    7 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    8 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    9 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    10 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    11 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    12 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    13 [ label = "Atom { element: Carbon, capacity: 0 }" ]
    0 -- 1 [ label = "Double" ]
    1 -- 2 [ label = "Single" ]
    2 -- 3 [ label = "Double" ]
    3 -- 4 [ label = "Single" ]
    4 -- 5 [ label = "Double" ]
    5 -- 6 [ label = "Single" ]
    6 -- 7 [ label = "Double" ]
    7 -- 8 [ label = "Single" ]
    8 -- 9 [ label = "Double" ]
    9 -- 10 [ label = "Single" ]
    10 -- 11 [ label = "Double" ]
    11 -- 12 [ label = "Single" ]
    12 -- 13 [ label = "Double" ]
    13 -- 0 [ label = "Single" ]
    12 -- 3 [ label = "Single" ]
    10 -- 5 [ label = "Single" ]
}"""
    assert at.molecule_info(mol) == anthra_info
