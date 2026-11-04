import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("analyze_ksi_complex", Path(__file__).resolve().parents[1] / "scripts/analyze_ksi_complex.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_coulomb_only_and_overlap():
    one = np.array([1.0])
    zero = np.array([0.0])
    assert module.pairwise_energy(np.array([[1.0]]), one, -one, one, one, zero, zero) == pytest.approx(-138.935456)
    with pytest.raises(ValueError):
        module.pairwise_energy(np.array([[0.0]]), one, one, one, one, zero, zero)
