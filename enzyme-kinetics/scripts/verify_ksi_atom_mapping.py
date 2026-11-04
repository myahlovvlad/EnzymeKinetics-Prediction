"""Cross-check 5-AND atom order and stereochemistry before/after minimization."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase


def stereo_smiles(mol, coords):
    clone = Chem.Mol(mol)
    for idx, point in enumerate(coords):
        clone.GetConformer().SetAtomPosition(idx, point)
    Chem.AssignStereochemistryFrom3D(clone, replaceExistingTags=True)
    return Chem.MolToSmiles(Chem.RemoveHs(clone), isomericSmiles=True)


def main():
    from openmm import app, unit

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sdf", type=Path, required=True)
    p.add_argument("--xyz", type=Path, required=True)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--assignment", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose new output")
    mol = next((m for m in Chem.SDMolSupplier(str(args.sdf), removeHs=False) if m), None)
    if mol is None or mol.GetNumAtoms() != 47:
        raise ValueError("Unexpected SDF")
    rows = [x.split() for x in args.xyz.read_text().splitlines()[2:]]
    if [r[0] for r in rows] != [a.GetSymbol() for a in mol.GetAtoms()]:
        raise ValueError("XYZ atom sequence mismatch")
    xyz = np.array([[float(v) for v in r[1:]] for r in rows])
    sdf_xyz = np.array(mol.GetConformer().GetPositions())
    pdb = app.PDBFile(str(args.pdb))
    lig = [a for a in pdb.topology.atoms() if a.residue.name == "UNK" and a.residue.chain.id == "C"]
    if [a.element.symbol for a in lig] != [a.GetSymbol() for a in mol.GetAtoms()]:
        raise ValueError("PDB ligand atom sequence mismatch")
    pos = np.asarray(pdb.positions.value_in_unit(unit.angstrom))[[a.index for a in lig]]
    assignment = json.loads(args.assignment.read_text())
    idx = assignment["poses"][0]["beta_h_index"]
    c4 = assignment["atom_names_zero_based"]["C4"]
    if mol.GetAtomWithIdx(idx).GetSymbol() != "H" or mol.GetBondBetweenAtoms(idx, c4) is None:
        raise ValueError("Assigned H is not bonded to C4 in SDF graph")
    input_s = stereo_smiles(mol, sdf_xyz)
    output_s = stereo_smiles(mol, pos)
    report = {"input_sha256": {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in (args.sdf, args.xyz, args.pdb, args.assignment)},
              "rdkit_version": rdBase.rdkitVersion, "xyz_sdf_max_coordinate_difference_A": float(np.max(np.linalg.norm(xyz-sdf_xyz, axis=1))),
              "c4_index": c4, "beta_h_index": idx, "beta_h_bonded_to_c4_in_sdf": True,
              "initial_C4_H_A": float(np.linalg.norm(sdf_xyz[c4]-sdf_xyz[idx])),
              "minimized_C4_H_A": float(np.linalg.norm(pos[c4]-pos[idx])),
              "sdf_3d_stereo_smiles": input_s, "minimized_3d_stereo_smiles": output_s,
              "stereo_preserved": input_s == output_s,
              "limitations": ["PDB and XYZ atom correspondence also relies on preserved writer order; identical elements are indistinguishable by sequence alone", "not an experimental confirmation of binding"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
