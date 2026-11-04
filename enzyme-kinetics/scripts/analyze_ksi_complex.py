"""Check minimized candidate geometry and direct pairwise force-field interaction.

Direct nonperiodic protein–ligand Coulomb + LJ is not a binding free energy,
not the PME decomposition and not a solvent-screened interaction energy.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def pairwise_energy(r, q1, q2, s1, s2, e1, e2):
    if np.any(r <= 0):
        raise ValueError("Overlapping atom positions")
    sigma = (s1[:, None] + s2[None, :]) / 2
    epsilon = np.sqrt(e1[:, None] * e2[None, :])
    ratio = sigma / r
    return float(
        np.sum(
            138.935456 * q1[:, None] * q2[None, :] / r
            + 4 * epsilon * (ratio**12 - ratio**6)
        )
    )


def main():
    import openmm
    from openmm import app, unit

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--system", type=Path, required=True)
    p.add_argument("--beta-assignment", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose a new output")
    pdb = app.PDBFile(str(args.pdb))
    system = openmm.XmlSerializer.deserialize(args.system.read_text())
    atoms = list(pdb.topology.atoms())
    if system.getNumParticles() != len(atoms):
        raise ValueError("Topology and force field particle counts differ")
    ligand = [
        a.index for a in atoms if a.residue.name == "UNK" and a.residue.chain.id == "C"
    ]
    protein = [
        a.index
        for a in atoms
        if a.residue.chain.id in ("A", "B")
        and a.residue.name not in ("HOH", "NA", "CL")
    ]
    if len(ligand) != 47 or len(protein) < 3000:
        raise ValueError("Unexpected protein/ligand topology")
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(nb) != 1:
        raise ValueError("Expected exactly one NonbondedForce")
    params = [nb[0].getParticleParameters(i) for i in range(len(atoms))]
    charges = np.array([v[0].value_in_unit(unit.elementary_charge) for v in params])
    sigma = np.array([v[1].value_in_unit(unit.nanometer) for v in params])
    epsilon = np.array([v[2].value_in_unit(unit.kilojoule_per_mole) for v in params])
    xyz = pdb.positions.value_in_unit(unit.nanometer)
    xyz = np.asarray(xyz)
    distances = np.linalg.norm(xyz[ligand, None, :] - xyz[None, protein, :], axis=2)
    direct = pairwise_energy(
        distances,
        charges[ligand],
        charges[protein],
        sigma[ligand],
        sigma[protein],
        epsilon[ligand],
        epsilon[protein],
    )
    assignment = json.loads(args.beta_assignment.read_text())
    beta = assignment["poses"][0]["beta_h_index"]
    c4 = assignment["atom_names_zero_based"]["C4"]
    oxygens = [
        a
        for a in atoms
        if a.residue.chain.id == "A"
        and a.residue.id == "38"
        and a.residue.name == "ASP"
        and a.name in ("OD1", "OD2")
    ]
    if (
        len(oxygens) != 2
        or atoms[ligand[beta]].element.symbol != "H"
        or atoms[ligand[c4]].element.symbol != "C"
    ):
        raise ValueError("Unexpected atom mapping for beta-H geometry")
    rows = []
    for o in oxygens:
        c, h, v = xyz[ligand[c4]], xyz[ligand[beta]], xyz[o.index]
        u, w = c - h, v - h
        angle = float(
            np.degrees(
                np.arccos(
                    np.clip(u @ w / (np.linalg.norm(u) * np.linalg.norm(w)), -1, 1)
                )
            )
        )
        rows.append(
            {
                "oxygen": o.name,
                "C4_O_A": float(np.linalg.norm(c - v) * 10),
                "Hbeta_O_A": float(np.linalg.norm(h - v) * 10),
                "C4_Hbeta_O_deg": angle,
            }
        )
    report = {
        "input_sha256": {
            str(x): hashlib.sha256(x.read_bytes()).hexdigest()
            for x in (args.pdb, args.system, args.beta_assignment)
        },
        "ligand_atom_count": len(ligand),
        "protein_atom_count": len(protein),
        "beta_h_local_index": beta,
        "geometry": rows,
        "direct_protein_ligand_coulomb_lj_kj_mol": direct,
        "method": "untruncated nonperiodic Coulomb + Lorentz-Berthelot LJ using force-field charges; protein–ligand only",
        "limitations": [
            "not PME energy, not solvent screened, not binding free energy",
            "PDB coordinates rounded to 0.001 A",
            "one unvalidated pose and protonation state",
            "beta-H assignment is geometry-based and order-dependent",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
