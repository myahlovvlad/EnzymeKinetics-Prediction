import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "assign_ksi_beta_h",
    Path(__file__).resolve().parents[1] / "scripts/assign_ksi_beta_h.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_beta_assignment_reflection_and_rotation():
    xyz = np.array(
        [[0.0, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, -1]]
    )
    names = {"C3": 0, "C4": 1, "C5": 2, "C10": 3, "C19": 4}
    assert module.identify_face(xyz, names, [5, 6])["beta_h_index"] == 5
    xyz[:, 2] *= -1
    assert module.identify_face(xyz, names, [5, 6])["beta_h_index"] == 5


def test_ambiguous_plane_rejected():
    with pytest.raises(ValueError):
        module.identify_face(
            np.zeros((7, 3)), {"C3": 0, "C4": 1, "C5": 2, "C10": 3, "C19": 4}, [5, 6]
        )
