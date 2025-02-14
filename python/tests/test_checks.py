import orca
import rdkit.Chem

def test_aspirin_index():
    aspirin_fname = "data/checks/aspirin.mol"
    mol = rdkit.Chem.MolFromMolFile(aspirin_fname)
    assert orca.compute_ma(mol) == 8

def test_aspirin_verbose():
    aspirin_fname = "data/checks/aspirin.mol"
    mol = rdkit.Chem.MolFromMolFile(aspirin_fname)
    assert orca.compute_ma_verbose(mol) == {'duplicates': 20, 'index': 8, 'space': 96}