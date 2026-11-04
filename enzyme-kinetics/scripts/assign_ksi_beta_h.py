"""Assign C4 H faces relative to angular C19 in a modelled 5-AND conformation.

The beta face is defined relative to the C19 angular methyl (steroid convention).
This geometry-dependent assignment is not experimental validation of a pose.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

SMILES = "C[C@]12CCC(=O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"


def identify_face(xyz, names, hydrogens):
    c3, c4, c5, c10, c19 = (xyz[names[n]] for n in ("C3", "C4", "C5", "C10", "C19"))
    normal = np.cross(c5 - c3, c10 - c3)
    if np.linalg.norm(normal) < 1e-6:
        raise ValueError("Degenerate steroid ring plane")
    normal /= np.linalg.norm(normal)
    methyl = float(normal @ (c19 - c10))
    projections = {int(h): float(normal @ (xyz[h] - c4)) for h in hydrogens}
    if abs(methyl) < 0.15 or len(hydrogens) != 2 or any(abs(v) < 0.15 for v in projections.values()):
        raise ValueError("Ambiguous face geometry")
    beta = [h for h, value in projections.items() if value * methyl > 0]
    if len(beta) != 1:
        raise ValueError("Expected exactly one hydrogen on beta face")
    return {"beta_h_index": beta[0], "alpha_h_index": next(h for h in projections if h != beta[0]),
            "methyl_normal_projection_A": methyl, "h_normal_projections_A": projections}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", type=Path, required=True)
    p.add_argument("--pose-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    lines = args.template.read_text().splitlines()
    block = "\n".join(line for line in lines if line.startswith("HETATM") and line[17:20] == "ASD" and line[21] == "A") + "\nEND\n"
    mol = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(SMILES), Chem.MolFromPDBBlock(block, removeHs=True))
    Chem.AssignStereochemistryFrom3D(mol)
    if Chem.MolToSmiles(mol) != Chem.MolToSmiles(Chem.MolFromSmiles(SMILES)):
        raise ValueError("Template-generated ligand stereo mismatch")
    names = {a.GetPDBResidueInfo().GetName().strip(): a.GetIdx() for a in mol.GetAtoms()}
    mol = Chem.AddHs(mol, addCoords=True)
    hs = [a.GetIdx() for a in mol.GetAtomWithIdx(names["C4"]).GetNeighbors() if a.GetAtomicNum() == 1]
    files = sorted(args.pose_dir.glob("candidate-*.xyz"))
    if not files:
        raise ValueError("No candidate XYZ files")
    results = []
    for path in files:
        lines = path.read_text().splitlines()
        entries = [line.split() for line in lines[2:]]
        if int(lines[0]) != len(entries) or len(entries) != mol.GetNumAtoms():
            raise ValueError("XYZ does not match template-derived atom count")
        if [r[0] for r in entries] != [a.GetSymbol() for a in mol.GetAtoms()]:
            raise ValueError("XYZ element sequence differs from template-derived ligand")
        xyz = np.array([[float(x) for x in r[1:]] for r in entries])
        result = identify_face(xyz, names, hs)
        results.append({"file": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), **result})
    report = {"template_sha256": hashlib.sha256(args.template.read_bytes()).hexdigest(),
              "rdkit_version": rdBase.rdkitVersion, "atom_names_zero_based": {k: names[k] for k in ("C3", "C4", "C5", "C10", "C19")},
              "method": "C4-H and C10-C19 projections on normal to plane C3/C5/C10; same side as C19 = beta",
              "limitations": ["Geometry-based model assignment; local ring plane is an approximation.",
                              "Atom order verified by element sequence, not by coordinates or full graph of XYZ.",
                              "Requires independent graphical or stereochemical verification before QM/MM."], "poses": results}
    if args.output.exists():
        raise FileExistsError("Choose new output")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
