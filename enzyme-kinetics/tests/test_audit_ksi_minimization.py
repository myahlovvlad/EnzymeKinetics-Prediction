import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("audit_ksi_minimization", Path(__file__).resolve().parents[1] / "scripts/audit_ksi_minimization.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_force_stats_detect_nonfinite():
    assert module.force_stats([[3, 4, 0]])["max_atom_norm_kj_mol_nm"] == 5
    with pytest.raises(ValueError):
        module.force_stats([[np.nan, 0, 0]])


def test_reaction_geometry():
    xyz = np.array([[0., 0., 0.], [.1, 0., 0.], [.3, 0., 0.]])
    result = module.reaction_geometry(xyz, 0, 1, [("O", 2)])[0]
    assert result["c4_o_A"] == pytest.approx(3.)
    assert result["h_o_A"] == pytest.approx(2.)
    assert result["angle_deg"] == pytest.approx(180.)
