"""Compare the selected assay record with saved reference structures.

Steroid atom numbers below apply only to the exact input record; indices are
checked against the expected molecular graph, not inferred for arbitrary SMILES.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import inchi, rdMolDescriptors


def molecular_info(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return mol, {"isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
                 "inchikey": inchi.MolToInchiKey(mol),
                 "formula": rdMolDescriptors.CalcMolFormula(mol),
                 "formal_charge": Chem.GetFormalCharge(mol),
                 "stereocenters": [{"index": i, "cip": c} for i, c in centers]}


def compare(actual, reference):
    _, a = molecular_info(actual)
    _, b = molecular_info(reference["smiles"])
    if b["inchikey"] != reference["inchikey"] or b["formula"] != reference["formula"]:
        raise ValueError("Reference transcription does not match stated PubChem properties")
    return {"local": a, "reference": b,
            "same_isomeric_smiles": a["isomeric_smiles"] == b["isomeric_smiles"],
            "same_inchikey": a["inchikey"] == b["inchikey"]}


def check_reaction_site(substrate, product):
    s = Chem.MolFromSmiles(substrate)
    p = Chem.MolFromSmiles(product)
    if s is None or p is None or s.GetNumAtoms() != 21 or p.GetNumAtoms() != 21:
        raise ValueError("Unexpected molecules")
    # The selected CSV entry retains identical atom order across both sides.
    s_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in s.GetBonds()}
    p_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in p.GetBonds()}
    if s_edges != p_edges or any(s.GetAtomWithIdx(i).GetSymbol() != p.GetAtomWithIdx(i).GetSymbol() for i in range(21)):
        raise ValueError("Atom indices cannot be paired by input order")
    labels = {"C3": 4, "C4": 6, "C5": 7, "C6": 8, "O3": 5}
    if any(s.GetAtomWithIdx(idx).GetSymbol() != ("O" if name == "O3" else "C") for name, idx in labels.items()):
        raise ValueError("Unexpected atom identities")
    pairs = [(4, 5), (4, 6), (6, 7), (7, 8)]
    expected_s = [2, 1, 1, 2]
    expected_p = [2, 1, 2, 1]
    def orders(m):
        return [m.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() if m.GetBondBetweenAtoms(a, b) else None for a, b in pairs]
    if orders(s) != expected_s or orders(p) != expected_p:
        raise ValueError("Unexpected reaction-site bond orders")
    return {"local_zero_based_indices": labels, "substrate_bond_orders": orders(s),
            "product_bond_orders": orders(p), "bond_pairs": pairs,
            "mapping_status": "input-order graph correspondence; not validated stereospecific hydrogen atom mapping",
            "transferable_beta_hydrogen": "UNKNOWN"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--row", type=int, required=True, help="1-based data row")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.row < 1:
        parser.error("row must be >= 1")
    with args.csv.open(newline="", encoding="utf-8-sig") as f:
        row = next((r for i, r in enumerate(csv.DictReader(f), 1) if i == args.row), None)
    if row is None or row.get("uniprot") != "P00947":
        raise ValueError("Missing row or wrong UniProt")
    s, p = row["reaction_smiles"].split(">>")
    reference = json.loads(args.reference.read_text())
    result = {"rdkit_version": rdBase.rdkitVersion, "source_row_1_based": args.row,
              "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (args.csv, args.reference)},
              "substrate_comparison": compare(s, reference["substrate"]),
              "product_comparison": compare(p, reference["product"]),
              "reaction_site": check_reaction_site(s, p),
              "limitations": ["PubChem identity does not confirm stereochemical identity of the precise assay lot.",
                              "Index correspondence is not a validated hydrogen-transfer mechanism."]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError("Choose a new output path; existing results are not overwritten")
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
