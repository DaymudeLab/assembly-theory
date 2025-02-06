import orca
import rdkit.Chem

def test_aspirin_index():
    aspirin_fname = "tests/inputs/aspirin.mol"
    mol = rdkit.Chem.MolFromMolFile(aspirin_fname)
    assert orca.compute_ma(mol) == 8