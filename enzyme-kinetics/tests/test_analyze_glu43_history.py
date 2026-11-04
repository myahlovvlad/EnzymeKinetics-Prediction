import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("openmm")
scripts = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(scripts))
try:
    spec = importlib.util.spec_from_file_location(
        "analyze_glu43_history", scripts / "analyze_glu43_history.py"
    )
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
finally:
    sys.path.remove(str(scripts))


def test_torsion_invariant_to_periodic_image():
    p = np.array([[0.1, 0.2, 0.3], [0.2, 0.2, 0.3], [0.2, 0.3, 0.3], [0.2, 0.3, 0.4]])
    result = m.dihedral(*p, np.array([2.0, 2.0, 2.0]))
    p[3] += np.array([2.0, 0.0, 0.0])
    assert m.dihedral(*p, np.array([2.0, 2.0, 2.0])) == pytest.approx(result)


def test_degenerate_torsion_rejected():
    with pytest.raises(ValueError):
        m.dihedral([0, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [3, 3, 3])
