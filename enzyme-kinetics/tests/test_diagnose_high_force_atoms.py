import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("openmm")
scripts = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(scripts))
try:
    spec = importlib.util.spec_from_file_location("diagnose_high_force_atoms", scripts / "diagnose_high_force_atoms.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
finally:
    sys.path.remove(str(scripts))


def test_summarize_categories():
    result = module.summarize(np.array([10., 60., 70., 1.]), np.array(["protein", "protein", "water", "ligand"]))
    assert result["protein"]["above_50"] == 1
    assert result["protein"]["fraction_of_all_above_50"] == .5
    assert result["water"]["max"] == 70
    assert result["ion"]["max"] is None
