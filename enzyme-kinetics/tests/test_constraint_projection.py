import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("constraint_projection", Path(__file__).resolve().parents[1] / "scripts/constraint_projection.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_analytic_pair_projection():
    x = [[0, 0, 0], [1, 0, 0]]
    f = [[3, 2, 0], [-1, 4, 5]]
    result, info = m.project_forces(x, f, [[0, 1]])
    np.testing.assert_allclose(result, [[1, 2, 0], [1, 4, 5]], atol=1e-10)
    assert info["tangent_residual_max"] < 1e-10


def test_dependent_and_degenerate_constraints():
    result, _ = m.project_forces([[0, 0, 0], [1, 0, 0]], [[3, 0, 0], [-1, 0, 0]], [[0, 1], [1, 0]])
    np.testing.assert_allclose(result[:, 0], [1, 1], atol=1e-10)
    with pytest.raises(ValueError):
        m.project_forces([[0, 0, 0], [0, 0, 0]], [[0]*3]*2, [[0, 1]])


def test_three_consecutive_and_energy_threshold():
    b = {"iterations": 100, "stop_reason": "converged", "energy_change_kJ_mol": 0.01, "projected_force_max": 2.0, "max_displacement_angstrom": 0.001}
    assert not m.convergence_confirmed([b])
    assert m.convergence_confirmed([b, b, b])
    assert not m.convergence_confirmed([b, dict(b, energy_change_kJ_mol=0.1), b])
    assert not m.convergence_confirmed([b, b, dict(b, stop_reason="max_iterations")])
