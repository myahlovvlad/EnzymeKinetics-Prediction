"""Integration checks for the preserved experimental modules."""

import numpy as np
from rdkit import Chem

from enzymekinetics_prediction.experiments import (
    core_chemistry,
    core_physics,
    core_structure,
)


def test_retrosynthetic_network():
    """Проверяем, что RDKit находит 3 сложноэфирные связи и не падает."""
    smiles = "CCC(=O)OCC(COC(=O)CC)OC(=O)CC"
    reactions = core_chemistry.get_retrosynthetic_network(smiles)

    # Триацетин должен дать 3 пути
    assert len(reactions) == 3, f"Ожидалось 3 пути, получено {len(reactions)}"
    for rxn in reactions:
        assert "site" in rxn
        assert "final_products" in rxn
        assert "." in rxn["final_products"]  # Должны быть разделенные продукты
        fragments = rxn["final_products"].split(".")
        assert len(fragments) == 2
        assert all(Chem.MolFromSmiles(fragment) is not None for fragment in fragments)


def test_protein_loading(tmp_path):
    """Use a small PDB fixture so a missing dataset cannot silently pass."""
    pdb_path = tmp_path / "protein.pdb"
    lines = []
    for serial, residue in enumerate([700, 701, 725], start=1):
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA A{residue:4d}    "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C  "
        )
    pdb_path.write_text("\n".join(lines) + "\nTER\nEND\n")
    coords, colors, lid, lid_idx = core_structure.load_protein_and_lid(pdb_path)
    assert coords.shape == (3, 3)
    assert len(colors) == 3
    assert lid.shape == (2, 3)
    assert lid_idx == [0, 1]


def test_time_warping():
    """Проверяем математику деформации времени."""
    t = core_physics.get_time_warping(100)
    assert len(t) == 100, "Неверная длина массива времени"
    assert t[0] == 0.0, "Время должно начинаться с 0"
    assert abs(t[-1] - 1.0) < 1e-5, "Время должно заканчиваться на 1"


def test_evb_energy():
    """Проверяем, что энергия возвращает число (float)."""
    e = core_physics.evb_energy(
        x=2.0, y=2.0, temp=298.15, pH=7.0, allosteric_factor=1.0, interfacial_factor=1.0
    )
    assert isinstance(e, (float, np.floating)), "Энергия должна быть числом"


def test_diffusion_association_uses_litres_per_mole_per_second():
    from enzymekinetics_prediction.experiments.kinetics_estimator import estimate_kon
    import pytest

    # Two 1-nm radii; relative diffusivity 1e-9 m²/s.
    # The diffusion-controlled value is ~1.5135e10 L/(mol*s), not ~1.5135e4.
    assert estimate_kon(1.0, 1.0, 5e-10, 5e-10) == pytest.approx(1.5135e10, rel=1e-4)
