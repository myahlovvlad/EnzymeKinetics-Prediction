"""Integration checks for relocated P2Rank-to-feature data contracts."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


def load_extractor():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts/structure/extract_struct_features_advanced.py"
    )
    spec = importlib.util.spec_from_file_location("pocket_extractor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_extractor_resolves_versioned_prediction_to_normalized_pdb(tmp_path):
    module = load_extractor()
    module.PDB_DIR = str(tmp_path)
    lines = []
    for serial, residue, name in [(1, 700, "LYS"), (2, 701, "ASP")]:
        lines.append(
            f"ATOM  {serial:5d}  CA  {name} A{residue:4d}    "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C  "
        )
    (tmp_path / "P12345.pdb").write_text("\n".join(lines) + "\nEND\n")
    predictions = tmp_path / "AF-P12345-F1-model_v6.pdb_predictions.csv"
    pd.DataFrame({"rank": [1], "score": [3.3], "residue_ids": ["A_700 A_701"]}).to_csv(
        predictions, index=False
    )
    result = module.process_single_pdb(str(predictions))
    assert result["freq_LYS"] == 0.5
    assert result["freq_ASP"] == 0.5
    assert result["pocket_charge"] == 0.0


def test_missing_pocket_structure_is_nan_instead_of_physical_zero(tmp_path):
    module = load_extractor()
    module.PDB_DIR = str(tmp_path)
    predictions = tmp_path / "AF-P12345-F1-model_v6.pdb_predictions.csv"
    pd.DataFrame({"rank": [1], "score": [3.3], "residue_ids": ["A_700 A_701"]}).to_csv(
        predictions, index=False
    )
    result = module.process_single_pdb(str(predictions))
    assert result["pocket_size"] == 2
    assert np.isnan(result["freq_LYS"])
    assert np.isnan(result["pocket_charge"])
