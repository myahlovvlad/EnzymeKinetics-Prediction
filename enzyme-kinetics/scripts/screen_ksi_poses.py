"""Exploratory rigid-body geometry screen, not docking or MM interaction energy."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

from Bio.PDB import PDBParser
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enzymekinetics_prediction.research.proton_surface import build_cluster, SMILES  # noqa: E402


def angle_degrees(a, h, b):
    x, y = a - h, b - h
    return float(
        np.degrees(
            np.arccos(np.clip(x @ y / (np.linalg.norm(x) * np.linalg.norm(y)), -1, 1))
        )
    )


def assess(xyz, c4, hydrogens, oxygens):
    combos = [
        {
            "hydrogen_index": int(h),
            "oxygen_index": int(o),
            "c4_oxygen_A": float(np.linalg.norm(xyz[c4] - xyz[o])),
            "h_oxygen_A": float(np.linalg.norm(xyz[h] - xyz[o])),
            "angle_c4_h_oxygen_deg": angle_degrees(xyz[c4], xyz[h], xyz[o]),
        }
        for h in hydrogens
        for o in oxygens
    ]
    return sorted(combos, key=lambda x: x["h_oxygen_A"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wt", type=Path, required=True)
    p.add_argument("--template", type=Path, required=True)
    p.add_argument(
        "--water-source",
        type=Path,
        help="8CHO assembly in same coordinate frame; water atoms from all models",
    )
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=[101, 202, 303])
    args = p.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Output must be empty; existing results are never overwritten")
    if len(set(args.seeds)) != len(args.seeds):
        raise ValueError("Seeds must be distinct")
    z, xyz, _, _, _, _, _, info = build_cluster(args.wt, args.template)
    n = Chem.AddHs(Chem.MolFromSmiles(SMILES)).GetNumAtoms()
    ligand = xyz[:n].copy()
    protein = PDBParser(QUIET=True).get_structure("wt", args.wt)[0]
    wt = protein["A"]
    oxygens = np.array([wt[38][name].coord.astype(float) for name in ("OD1", "OD2")])
    water_atoms = [
        atom
        for atom in protein.get_atoms()
        if atom.get_parent().id[0] == "W" and atom.element not in {"H", "D"}
    ]
    if args.water_source:
        if water_atoms:
            raise ValueError(
                "WT already has waters; refusing ambiguous duplicate water sources"
            )
        water_atoms = [
            atom
            for atom in PDBParser(QUIET=True)
            .get_structure("water", args.water_source)
            .get_atoms()
            if atom.get_parent().id[0] == "W" and atom.element not in {"H", "D"}
        ]
    receptor = np.array(
        [
            atom.coord
            for atom in protein.get_atoms()
            if atom.element not in {"H", "D"} and atom.get_parent().id[0] == " "
        ]
        + [atom.coord for atom in water_atoms],
        float,
    )
    receptor = receptor[
        np.linalg.norm(receptor - ligand[:21].mean(axis=0), axis=1) < 12
    ]
    # Reconstruct explicit ligand bonds in template ordering to identify both C4 hydrogens.
    lines = args.template.read_text().splitlines()
    block = (
        "\n".join(
            line
            for line in lines
            if line.startswith("HETATM") and line[17:20] == "ASD" and line[21] == "A"
        )
        + "\nEND\n"
    )
    mol = AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(SMILES), Chem.MolFromPDBBlock(block, removeHs=True)
    )
    mol = Chem.AddHs(mol, addCoords=True)
    c4 = info["C4_index"]
    hs = [
        a.GetIdx()
        for a in mol.GetAtomWithIdx(c4).GetNeighbors()
        if a.GetAtomicNum() == 1
    ]
    if len(hs) != 2 or len(ligand) != mol.GetNumAtoms():
        raise ValueError("Unexpected ligand atom mapping or C4 hydrogen count")
    center = ligand.mean(axis=0)
    reference = ligand - center
    heavy = np.arange(21)
    # O atoms need smaller generic repulsive contact radii than carbon atoms.
    contact = np.where(np.asarray(z[:21]) == 8, 2.7, 3.05)

    def transformed(x):
        return Rotation.from_rotvec(x[3:]).apply(reference) + center + x[:3]

    def objective(x):
        q = transformed(x)
        d = np.linalg.norm(q[heavy, None, :] - receptor[None, :, :], axis=2)
        overlap = np.maximum(contact[:, None] - d, 0)
        donor = np.linalg.norm(q[c4] - oxygens, axis=1).min()
        return float(
            20 * (overlap**2).sum() + 8 * max(donor - 3.5, 0) ** 2 + 0.05 * (x @ x)
        )

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for seed in args.seeds:
        rng = np.random.default_rng(seed)
        start = np.r_[rng.normal(0, 0.7, 3), rng.normal(0, 0.25, 3)]
        fit = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=[(-3, 3)] * 3 + [(-1.5, 1.5)] * 3,
            options={"maxiter": 300},
        )
        q = transformed(fit.x)
        d = np.linalg.norm(q[heavy, None, :] - receptor[None, :, :], axis=2)
        allxyz = np.vstack([q, oxygens])
        geometries = assess(allxyz, c4, hs, [n, n + 1])
        chosen = geometries[0]
        result = {
            "seed": seed,
            "optimizer_success": bool(fit.success),
            "message": str(fit.message),
            "geometry_objective_arbitrary": float(fit.fun),
            "protein_heavy_min_distance_A": float(d.min()),
            "c4_h_donor_pairs": geometries,
            "screen_distance_pass": min(g["c4_oxygen_A"] for g in geometries) <= 4,
            "best_h_is_beta": "UNKNOWN",
            "mm_interaction_energy": "UNKNOWN",
            "closest_h_angle_deg": chosen["angle_c4_h_oxygen_deg"],
            "candidate_status": "unvalidated; restraint-biased; no verified beta-H",
        }
        rows.append(result)
        with (args.output / f"candidate-{seed}.xyz").open("w") as f:
            f.write(
                f"{len(q)}\n5-AND ligand only; rigid-body restrained geometry; not WT-ES validation\n"
            )
            for atom, pos in zip(z[:n], q, strict=True):
                f.write(
                    Chem.GetPeriodicTable().GetElementSymbol(int(atom))
                    + " "
                    + " ".join(f"{v:.7f}" for v in pos)
                    + "\n"
                )
    inputs = [args.wt, args.template] + (
        [args.water_source] if args.water_source else []
    )
    summary = {
        "inputs_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
        "rdkit_version": rdBase.rdkitVersion,
        "reference": info,
        "receptor_heavy_atoms_within_12A": len(receptor),
        "crystallographic_water_heavy_atoms_in_input": len(water_atoms),
        "experimental_water_included": bool(water_atoms),
        "method": "independent seeded rigid-body starts; geometric overlap and Asp38 distance restraint; no physical energy",
        "note": "Clashes scored against WT and any waters actually present, which are not relaxed. The supplied water source must share the WT coordinate frame; catalytic roles are not established. β-H stereochemistry unknown. Not a docking score, binding energy, or MD ensemble.",
        "candidates": rows,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
