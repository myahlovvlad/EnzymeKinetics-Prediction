import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("screen_ksi_poses", Path(__file__).resolve().parents[1] / "scripts/screen_ksi_poses.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_angle_and_candidate_geometry():
    xyz = np.array([[0., 0, 0], [1., 0, 0], [2.5, 0, 0]])
    result = module.assess(xyz, 0, [1], [2])[0]
    assert result["c4_oxygen_A"] == pytest.approx(2.5)
    assert result["h_oxygen_A"] == pytest.approx(1.5)
    assert result["angle_c4_h_oxygen_deg"] == pytest.approx(180)
