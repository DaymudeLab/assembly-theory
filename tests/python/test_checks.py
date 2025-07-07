import assembly_theory as at
import rdkit.Chem


def test_aspirin_index():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly(mol) == 8
