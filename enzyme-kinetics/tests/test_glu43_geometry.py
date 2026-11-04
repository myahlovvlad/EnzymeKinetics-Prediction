import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "glu43_geometry", Path(__file__).resolve().parents[1] / "scripts/glu43_geometry.py"
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_backbone_and_sidechain_distinct_and_periodic():
    atoms = {
        "C": [0.1, 0, 0],
        "O": [1.98, 0, 0],
        "CG": [1, 0, 0],
        "CD": [1, 1, 0],
        "OE1": [0, 1, 0],
        "OE2": [1, 1.5, 0],
    }
    result = m.geometry(atoms, [2, 2, 2])
    assert result["backbone_C_O_nm"] == pytest.approx(0.12)
    assert result["sidechain_CD_OE1_nm"] == pytest.approx(1.0)
    assert result["sidechain_CD_OE2_nm"] == pytest.approx(0.5)
    assert result["sidechain_OE1_CD_OE2_deg"] == pytest.approx(90.0)
    assert np.isfinite(result["sidechain_CG_CD_OE1_OE2_improper_deg"])


def test_reject_missing_and_degenerate_atoms():
    with pytest.raises(ValueError):
        m.geometry({"C": [0, 0, 0]}, [2, 2, 2])
