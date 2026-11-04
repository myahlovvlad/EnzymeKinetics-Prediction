import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "check_p00947_chemistry",
    Path(__file__).resolve().parents[1] / "scripts/check_p00947_chemistry.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

S = "C[C@]12CCC(=O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"
P = "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"


def test_reaction_site_and_indices():
    result = module.check_reaction_site(S, P)
    assert result["local_zero_based_indices"]["C4"] == 6
    assert result["local_zero_based_indices"]["C6"] == 8
    assert result["transferable_beta_hydrogen"] == "UNKNOWN"


def test_reversed_or_wrong_molecules_rejected():
    with pytest.raises(ValueError):
        module.check_reaction_site(P, S)


def test_enantiomer_or_epimer_not_exact_match():
    _, info = module.molecular_info(S)
    reference = {"smiles": S, "inchikey": info["inchikey"], "formula": info["formula"]}
    assert module.compare(S, reference)["same_inchikey"]
    assert not module.compare(S.replace("[C@]", "[C@@]", 1), reference)["same_inchikey"]
