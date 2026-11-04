import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "parameterize_ksi_ligand",
    Path(__file__).resolve().parents[1] / "scripts/parameterize_ksi_ligand.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_xyz_rejects_wrong_element_order(tmp_path):
    f = tmp_path / "ligand.xyz"
    f.write_text("2\ncomment\nC 0 0 0\nO 1 0 0\n")
    assert len(module.read_xyz(f, ["C", "O"])) == 2
    with pytest.raises(ValueError):
        module.read_xyz(f, ["O", "C"])
