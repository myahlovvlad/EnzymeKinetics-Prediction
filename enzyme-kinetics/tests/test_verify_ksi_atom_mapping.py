import importlib.util
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

spec = importlib.util.spec_from_file_location("verify_ksi_atom_mapping", Path(__file__).resolve().parents[1] / "scripts/verify_ksi_atom_mapping.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_stereo_smiles_is_reproducible():
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)F"))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    coords = mol.GetConformer().GetPositions()
    assert module.stereo_smiles(mol, coords) == module.stereo_smiles(mol, coords)
