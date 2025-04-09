import assembly_theory as at
import rdkit.Chem

def test_aspirin_index():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly(mol) == 8

def test_aspirin_verbose():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly_verbose(mol) == {'duplicates': 20, 'index': 8, 'space': 36}

    assert at.molecular_assembly_verbose(mol, no_bounds=True) == {'duplicates': 20, 'index': 8, 'space': 96}

def test_anthracene_info():
    anthra_smi = "c1ccc2cc3ccccc3cc2c1"
    mol = rdkit.Chem.MolFromSmiles(anthra_smi)

    anthra_info = f"""0: Carbon
1: Carbon
2: Carbon
3: Carbon
4: Carbon
5: Carbon
6: Carbon
7: Carbon
8: Carbon
9: Carbon
10: Carbon
11: Carbon
12: Carbon
13: Carbon

0: Double, (0, 1), (Carbon, Carbon)
1: Single, (1, 2), (Carbon, Carbon)
2: Double, (2, 3), (Carbon, Carbon)
3: Single, (3, 4), (Carbon, Carbon)
4: Double, (4, 5), (Carbon, Carbon)
5: Single, (5, 6), (Carbon, Carbon)
6: Double, (6, 7), (Carbon, Carbon)
7: Single, (7, 8), (Carbon, Carbon)
8: Double, (8, 9), (Carbon, Carbon)
9: Single, (9, 10), (Carbon, Carbon)
10: Double, (10, 11), (Carbon, Carbon)
11: Single, (11, 12), (Carbon, Carbon)
12: Double, (12, 13), (Carbon, Carbon)
13: Single, (13, 0), (Carbon, Carbon)
14: Single, (12, 3), (Carbon, Carbon)
15: Single, (10, 5), (Carbon, Carbon)
"""
    assert at.molecule_info(mol) == anthra_info

