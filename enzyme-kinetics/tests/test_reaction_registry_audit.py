"""Checks for the reaction registry's conservative chemistry accounting."""

import importlib.util
from pathlib import Path


path = Path(__file__).resolve().parents[1] / "scripts" / "audit_reaction_registry.py"
spec = importlib.util.spec_from_file_location("audit_reaction_registry", path)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_balanced_and_fragment_order():
    first = audit.inspect_reaction("O.CC>>CC.O")
    second = audit.inspect_reaction("CC.O>>O.CC")
    assert first["canonical_reaction_smiles"] == second["canonical_reaction_smiles"]
    assert first["atom_balanced"] and first["charge_balanced"]


def test_hydrogen_and_charge_accounted_separately():
    result = audit.inspect_reaction("O>>[OH-]")
    assert result["heavy_element_balanced"]
    assert result["element_delta_products_minus_reactants"] == {"H": -1}
    assert result["charge_delta_products_minus_reactants"] == -1


def test_invalid_reaction_is_not_balanced():
    assert audit.inspect_reaction("C>>")["parse_status"] == "incomplete"
    assert audit.inspect_reaction("C>>")["reactants"]["elements"] == {"C": 1, "H": 4}
    assert audit.inspect_reaction("C>O>C")['parse_status'] == "error"


def test_stoichiometry_is_not_deduplicated():
    result = audit.inspect_reaction("O.O>>O")
    assert result["element_delta_products_minus_reactants"] == {"H": -2, "O": -1}
