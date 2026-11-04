import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "pbc_clash_audit",
    Path(__file__).resolve().parents[1] / "scripts/pbc_clash_audit.py",
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_minimum_image_across_boundary_and_shift():
    box = np.array([2.0, 2.0, 2.0])
    first, shift = m.minimum_image(np.array([1.9, 0, 0]), box)
    second, _ = m.minimum_image(np.array([3.9, 0, 0]), box)
    np.testing.assert_allclose(first, [-0.1, 0, 0])
    np.testing.assert_allclose(second, first)
    np.testing.assert_equal(shift, [-1, 0, 0])


def test_exclusion_water_and_interprotein():
    assert m.classify_pair("protein", "protein", "C", "N", False, True) is None
    assert m.classify_pair("water", "water", "O", "H", True, False) is None
    assert m.classify_pair("water", "water", "H", "H", False, False) is None
    assert m.classify_pair("water", "water", "O", "H", False, False) == "water-water"
    assert (
        m.classify_pair("protein", "protein", "C", "N", False, False)
        == "protein-protein"
    )
    assert 1.4 < 1.5 and 1.7 < 2.0 and not 1.7 < 1.5
