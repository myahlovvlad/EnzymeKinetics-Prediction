import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("mdtraj")


spec = importlib.util.spec_from_file_location(
    "diagnose_dimer_rmsd", Path(__file__).resolve().parents[1] / "scripts/diagnose_dimer_rmsd.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_periodic_jump_does_not_move_dimer_in_reconstruction():
    xyz = np.array([[[0., 0., 0.], [1., 0., 0.]], [[0., 0., 0.], [11., 0., 0.]]])
    fixed, shift = module.reconstruct_pair(xyz, np.full((2, 3), 10.), [0], [1])
    np.testing.assert_allclose(fixed[0], fixed[1])
    np.testing.assert_allclose(shift[1], [10, 0, 0])
    assert xyz[1, 1, 0] == 11


def test_missing_box_rejected():
    with pytest.raises(ValueError):
        module.reconstruct_pair(np.zeros((2, 2, 3)), np.zeros((2, 3)), [0], [1])
