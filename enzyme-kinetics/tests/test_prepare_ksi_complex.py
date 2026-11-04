import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

spec = importlib.util.spec_from_file_location(
    "prepare_ksi_complex",
    Path(__file__).resolve().parents[1] / "scripts/prepare_ksi_complex.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_variants_require_both_chains():
    residues = [
        SimpleNamespace(name="ASP", id="99", chain=SimpleNamespace(id=c)) for c in "AB"
    ]
    top = SimpleNamespace(residues=lambda: iter(residues))
    assert module.select_ash_variants(top) == ["ASH", "ASH"]
    residues.pop()
    with pytest.raises(ValueError):
        module.select_ash_variants(top)
