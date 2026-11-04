"""Small integration checks for the native chemistry dependencies."""

import numpy as np
import pytest


def test_rdkit_conformer_to_pyscf():
    chem = pytest.importorskip("rdkit.Chem")
    all_chem = pytest.importorskip("rdkit.Chem.AllChem")
    gto = pytest.importorskip("pyscf.gto")
    molecule = chem.AddHs(chem.MolFromSmiles("O"))
    assert all_chem.EmbedMolecule(molecule, randomSeed=42) == 0
    positions = molecule.GetConformer().GetPositions()
    assert np.isfinite(positions).all()
    atoms = [
        (atom.GetSymbol(), position)
        for atom, position in zip(molecule.GetAtoms(), positions, strict=True)
    ]
    quantum_molecule = gto.M(atom=atoms, basis="sto-3g", unit="Angstrom", verbose=0)
    assert quantum_molecule.natm == 3
    assert quantum_molecule.nelectron == 10


def test_pyscf_native_scf():
    gto = pytest.importorskip("pyscf.gto")
    scf = pytest.importorskip("pyscf.scf")
    lib = pytest.importorskip("pyscf.lib")
    previous_threads = lib.num_threads()
    try:
        lib.num_threads(1)
        molecule = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)
        method = scf.RHF(molecule)
        energy = method.kernel()
        assert method.converged
        assert energy == pytest.approx(-1.1167593, abs=1e-6)
    finally:
        lib.num_threads(previous_threads)
