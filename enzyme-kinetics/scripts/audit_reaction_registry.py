"""Audit the local kcat reaction records; no external data or invented chemistry.

Output cards represent assay rows, not independently verified reaction equations.
Each dot-separated SMILES component has implicit coefficient one in this audit.
"""

import argparse
from collections import Counter
import csv
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase


@lru_cache(maxsize=None)
def inspect_side(text):
    if not text:
        raise ValueError("empty side")
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise ValueError("invalid SMILES")
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        raise ValueError("empty side")
    canonical = sorted(Chem.MolToSmiles(m, isomericSmiles=True) for m in fragments)
    elements = Counter()
    charge = 0
    maps = []
    unspecified_stereo = 0
    for fragment in fragments:
        unspecified_stereo += sum(
            label == "?" for _, label in Chem.FindMolChiralCenters(
                fragment, includeUnassigned=True, includeCIP=True
            )
        )
        for atom in Chem.AddHs(fragment).GetAtoms():
            elements[atom.GetSymbol()] += 1
            charge += atom.GetFormalCharge()
        maps.extend(a.GetAtomMapNum() for a in fragment.GetAtoms() if a.GetAtomMapNum())
    return {
        "components": canonical,
        "elements": dict(elements),
        "charge": charge,
        "mapped_atoms": len(maps),
        "heavy_atoms": sum(m.GetNumHeavyAtoms() for m in fragments),
        "unspecified_tetrahedral_centers": unspecified_stereo,
    }


def inspect_reaction(text):
    try:
        parts = text.split(">>")
        if len(parts) != 2:
            raise ValueError("expected reactants>>products")
        if not all(parts):
            return {"parse_status": "incomplete", "error": "empty side",
                    "reactants": inspect_side(parts[0]) if parts[0] else None,
                    "products": inspect_side(parts[1]) if parts[1] else None}
        left, right = map(inspect_side, parts)
        keys = set(left["elements"]) | set(right["elements"])
        delta = {k: right["elements"].get(k, 0) - left["elements"].get(k, 0)
                 for k in sorted(keys)}
        delta = {k: v for k, v in delta.items() if v}
        charge_delta = right["charge"] - left["charge"]
        key = ".".join(left["components"]) + ">>" + ".".join(right["components"])
        return {"parse_status": "ok", "canonical_reaction_smiles": key,
                "reactants": left, "products": right,
                "element_delta_products_minus_reactants": delta,
                "charge_delta_products_minus_reactants": charge_delta,
                "atom_balanced": not delta, "heavy_element_balanced": not any(k != "H" for k in delta),
                "charge_balanced": charge_delta == 0}
    except (ValueError, RuntimeError) as exc:
        return {"parse_status": "error", "error": str(exc)}


def present(value):
    return bool(value and value.strip() and value.strip().lower() not in {"nan", "none", "null"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {"rdkit_version": rdBase.rdkitVersion, "files": {}, "counts": {},
               "limitations": ["Counts are assay rows, not independent reactions or experiments.",
                   "Canonical SMILES are representation-normalized only; not tautomer- or protonation-normalized.",
                   "Component multiplicities are as encoded; absent species cannot be reconstructed.",
                   "Balance uses implicit hydrogens and formal charges; isotopes are not separately counted.",
                   "Atom-map counts do not establish correct or complete atom mapping.",
                   "Missing metadata are unknown, not zero."]}
    totals = Counter()
    unique_raw, unique_canonical = set(), set()
    with (args.output_dir / "assay-cards.jsonl").open("w", encoding="utf-8") as out:
        for split in ("trainval", "test"):
            path = args.input_dir / f"kcat-random_{split}.csv"
            summary["files"][path.name] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            with path.open(encoding="utf-8-sig", newline="") as source:
                for row_number, row in enumerate(csv.DictReader(source), 1):
                    totals["rows"] += 1
                    totals[f"rows_{split}"] += 1
                    for field in ("reaction_smiles", "ph", "temperature", "uniprot", "ec", "sequence"):
                        if present(row.get(field)):
                            totals[f"present_{field}"] += 1
                    if present(row.get("ph")) and present(row.get("temperature")):
                        totals["present_both_ph_temperature"] += 1
                    raw = row.get("reaction_smiles", "")
                    if present(raw):
                        unique_raw.add(raw)
                    result = inspect_reaction(raw)
                    if result["parse_status"] == "ok":
                        totals["parsed"] += 1
                        unique_canonical.add(result["canonical_reaction_smiles"])
                        for field in ("atom_balanced", "heavy_element_balanced", "charge_balanced"):
                            totals[field] += int(result[field])
                        totals["atom_and_charge_balanced"] += int(result["atom_balanced"] and result["charge_balanced"])
                        totals["hydrogen_only_imbalance"] += int(result["heavy_element_balanced"] and not result["atom_balanced"])
                        totals["any_atom_map"] += int(bool(result["reactants"]["mapped_atoms"] or result["products"]["mapped_atoms"]))
                        totals["unspecified_tetrahedral_stereo"] += int(bool(result["reactants"]["unspecified_tetrahedral_centers"] or result["products"]["unspecified_tetrahedral_centers"]))
                    else:
                        totals["parse_errors"] += 1
                    card = {"source_file": path.name, "source_row_1_based": row_number,
                            "reaction_smiles": raw, "reaction_audit": result,
                            "enzyme": {k: row.get(k) or None for k in ("uniprot", "sequence", "sequence_source", "ec", "taxonomy_id")},
                            "assay": {k: row.get(k) or None for k in ("value", "log10_value", "ph", "temperature")},
                            "unknown_not_in_source_columns": ["buffer", "ionic_strength", "solvent", "concentrations", "activities", "oligomer", "variant_annotation", "assay_primary_reference", "kcat_unit", "km", "reverse_rate", "equilibrium", "measurement_error"],
                            "direction_note": "SMILES left-to-right encoding; physiological direction and reversibility unknown",
                            "mechanism_status": "not inferred", "experimental_validation": "not verified"}
                    out.write(json.dumps(card, ensure_ascii=False) + "\n")
    totals["unique_raw_reaction_smiles"] = len(unique_raw)
    totals["unique_canonical_reaction_smiles"] = len(unique_canonical)
    summary["counts"] = dict(totals)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
