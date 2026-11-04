import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("pyscf")
spec = importlib.util.spec_from_file_location(
    "check_qmmm_embedding",
    Path(__file__).resolve().parents[1] / "scripts/check_qmmm_embedding.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_embedding_converges_and_changes_energy(monkeypatch):
    from pyscf import lib

    lib.num_threads(1)
    result = module.run()
    assert result["scf_converged"]
    assert abs(result["energy_difference_hartree"]) > 1e-5
