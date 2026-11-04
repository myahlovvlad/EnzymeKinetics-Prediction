import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openmm")
spec = importlib.util.spec_from_file_location(
    "validate_opencl_gpu",
    Path(__file__).resolve().parents[1] / "scripts/validate_opencl_gpu.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_force_comparison_and_bad_values():
    result = module.compare_forces([[1.0, 0.0, 0.0]], [[2.0, 0.0, 0.0]])
    assert result["max_component_error_kj_mol_nm"] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        module.compare_forces([[1.0, 0.0, 0.0]], [[np.nan, 0.0, 0.0]])
