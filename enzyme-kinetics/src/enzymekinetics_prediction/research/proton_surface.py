"""Exploratory electronic proton-position surface of an explicitly capped KSI cluster.

A substrate pose is modelled from the 1OHS saturated analog by bond-order
assignment, stereochemistry validation and CA alignment to WT 8CHO. The protein
is NOT an experimentally observed ES. Only the transferred proton moves in a
2D affine plane; all other atoms are frozen. Stationary points in this plane
are NOT full enzyme transition states or free-energy barriers.
"""

from pathlib import Path

from Bio.PDB import PDBParser, Superimposer
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import root
from threadpoolctl import threadpool_limits

from enzymekinetics_prediction.physics.artifacts import write_json, digest_file
from enzymekinetics_prediction.physics.quantum import BOHR_A, HARTREE_EV
from enzymekinetics_prediction.physics.mechanism import (
    KSI_MECHANISM,
    classify_stationary,
    numerical_hessian,
)
from .journal import Journal

SMILES = "C[C@]12CCC(=O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"


def fragment(smiles, atom_names, residue):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom in mol.GetAtoms():
        conf.SetAtomPosition(
            atom.GetIdx(), residue[atom_names[atom.GetAtomMapNum()]].coord.astype(float)
        )
    mol.AddConformer(conf)
    return Chem.AddHs(mol, addCoords=True)


def build_cluster(reference, analog):
    parser = PDBParser(QUIET=True)
    wt = parser.get_structure("wt", reference)[0]["A"]
    donor = parser.get_structure("analog", analog)[0]["A"]
    sup = Superimposer()
    sup.set_atoms(
        [wt[i]["CA"] for i in range(1, 126)], [donor[i]["CA"] for i in range(1, 126)]
    )
    lines = Path(analog).read_text().splitlines()
    ligand_codes = {
        line[17:20] for line in lines if line.startswith("HETATM") and line[21] == "A"
    }
    ligand_code = (
        "5SD" if "5SD" in ligand_codes else "ASD" if "ASD" in ligand_codes else None
    )
    if ligand_code is None:
        raise ValueError("Template must contain 5SD or ASD in chain A")
    block = (
        "\n".join(
            line
            for line in lines
            if line.startswith("HETATM")
            and line[17:20] == ligand_code
            and line[21] == "A"
        )
        + "\nEND\n"
    )
    ligand = AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(SMILES), Chem.MolFromPDBBlock(block, removeHs=True)
    )
    Chem.AssignStereochemistryFrom3D(ligand)
    if Chem.MolToSmiles(ligand) != Chem.MolToSmiles(Chem.MolFromSmiles(SMILES)):
        raise ValueError("Transferred pose does not match target stereo-SMILES")
    names = {
        a.GetPDBResidueInfo().GetName().strip(): a.GetIdx() for a in ligand.GetAtoms()
    }
    rot, trans = sup.rotran
    conf = ligand.GetConformer()
    for i, xyz in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(i, xyz @ rot + trans)
    ligand = Chem.AddHs(ligand, addCoords=True)
    # Modest positional restraint regularizes the saturated analog coordinates.
    properties = AllChem.MMFFGetMoleculeProperties(ligand)
    ff = AllChem.MMFFGetMoleculeForceField(ligand, properties)
    for i in range(21):
        ff.MMFFAddPositionConstraint(i, 0.2, 100)
    status = ff.Minimize(maxIts=1000)
    if status:
        raise ValueError("Restrained substrate pose relaxation did not converge")
    tyrosine = fragment(
        "[CH3:1][c:2]1[cH:3][cH:4][c:5]([OH:8])[cH:6][cH:7]1",
        {1: "CB", 2: "CG", 3: "CD1", 4: "CE1", 5: "CZ", 6: "CE2", 7: "CD2", 8: "OH"},
        wt[14],
    )
    base = fragment(
        "[CH3:1][C:2](=[O:3])[O-:4]", {1: "CB", 2: "CG", 3: "OD1", 4: "OD2"}, wt[38]
    )
    acid = fragment(
        "[CH3:1][C:2](=[O:3])[OH:4]", {1: "CB", 2: "CG", 3: "OD1", 4: "OD2"}, wt[99]
    )
    pieces = [ligand, tyrosine, base, acid]
    positions = np.vstack([m.GetConformer().GetPositions() for m in pieces])
    numbers = np.array(
        [a.GetAtomicNum() for m in pieces for a in m.GetAtoms()], np.int32
    )
    offset = ligand.GetNumAtoms() + tyrosine.GetNumAtoms()
    oxygens = [offset + a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() == 8]
    hydrogens = [
        a.GetIdx()
        for a in ligand.GetAtomWithIdx(names["C4"]).GetNeighbors()
        if a.GetAtomicNum() == 1
    ]
    distance = [
        (np.linalg.norm(positions[h] - positions[o]), h, o)
        for h in hydrogens
        for o in oxygens
    ]
    _, proton, oxygen = min(distance)
    donor_distance = min(
        np.linalg.norm(positions[names["C4"]] - positions[o]) for o in oxygens
    )
    if donor_distance > 4.0:
        raise ValueError(
            f"Unproductive transferred analog pose: C4--Asp38 oxygen distance {donor_distance:.3f} A > 4.0 A; do not compute a reaction barrier from this pose"
        )
    c4, c6 = positions[names["C4"]], positions[names["C6"]]
    o = positions[oxygen]
    anchors = np.array(
        [
            positions[proton],
            o + (c4 - o) / np.linalg.norm(c4 - o),
            c6 + 1.09 * (o - c6) / np.linalg.norm(o - c6),
        ]
    )
    e1 = anchors[1] - anchors[0]
    e1 /= np.linalg.norm(e1)
    e2 = anchors[2] - anchors[0]
    e2 -= e1 * (e2 @ e1)
    e2 /= np.linalg.norm(e2)
    basis = np.column_stack([e1, e2])
    coordinates = (anchors - anchors[0]) @ basis
    charge = sum(Chem.GetFormalCharge(m) for m in pieces)
    return (
        numbers,
        positions,
        proton,
        anchors[0],
        basis,
        coordinates,
        charge,
        {
            "CA_alignment_rmsd_A": float(sup.rms),
            "source_analog": f"{ligand_code} in {Path(analog).name}; template ligand is not the target Delta5 substrate",
            "target_SMILES": SMILES,
            "C4_index": names["C4"],
            "C6_index": names["C6"],
            "Asp38_oxygen_index": oxygen,
            "transferred_proton_index": proton,
            "charge_e": charge,
            "fragments": [
                "substrate",
                "p-cresol Tyr14 cap",
                "acetate Asp38 cap",
                "acetic acid Asp99 cap",
            ],
        },
    )


def run(reference, analog, output, grid=13):
    from tblite.interface import Calculator

    output = Path(output)
    if (output / "surface.csv").exists():
        raise ValueError("Surface already exists; use a fresh experiment directory")
    output.mkdir(parents=True, exist_ok=True)
    journal = Journal(output)
    with journal.run(
        "H5-projected-proton-surface",
        {
            "grid": grid,
            "method": "GFN2-xTB",
            "solvent": "ALPB water",
            "mobile_dof": 2,
            "threads": 1,
        },
        {
            "reference_sha256": digest_file(reference),
            "analog_sha256": digest_file(analog),
        },
    ):
        z, xyz, proton, origin, basis, anchors, charge, preparation = build_cluster(
            reference, analog
        )
        if (int(z.sum()) - charge) % 2:
            raise ValueError("Unexpected open-shell cluster")
        write_json(
            output / "cluster.json",
            {
                "atomic_numbers": z.tolist(),
                "positions_A": xyz.tolist(),
                "origin_A": origin.tolist(),
                "basis": basis.tolist(),
                "anchors_A": anchors.tolist(),
                **preparation,
            },
        )
        with (output / "cluster.xyz").open("w") as f:
            f.write(f"{len(z)}\nModelled capped WT cluster; not an experimental ES\n")
            for number, position in zip(z, xyz, strict=True):
                f.write(
                    Chem.GetPeriodicTable().GetElementSymbol(int(number))
                    + " "
                    + " ".join(map(str, position))
                    + "\n"
                )
        evaluations = 0

        def calculate(q):
            nonlocal evaluations
            coordinates = xyz.copy()
            coordinates[proton] = origin + basis @ np.asarray(q)
            distances = np.linalg.norm(coordinates - coordinates[proton], axis=1)
            distances[proton] = 100
            if distances.min() < 0.35:
                raise ValueError("Proton grid point too close to a nucleus")
            calc = Calculator(
                "GFN2-xTB", z, coordinates / BOHR_A, charge=float(charge), uhf=0
            )
            calc.set("verbosity", 0)
            calc.set("max-iter", 250)
            calc.add("alpb-solvation", "water")
            with threadpool_limits(1):
                r = calc.singlepoint()
            energy = float(r.get("energy")) * HARTREE_EV
            gradient = np.asarray(r.get("gradient")) * HARTREE_EV / BOHR_A
            if not np.isfinite(energy) or not np.isfinite(gradient).all():
                raise ValueError("Nonfinite quantum result")
            evaluations += 1
            return energy, gradient[proton] @ basis, float(np.max(np.abs(gradient)))

        lo = anchors.min(axis=0) - 0.25
        hi = anchors.max(axis=0) + 0.25
        rows = []
        for x in np.linspace(lo[0], hi[0], grid):
            for y in np.linspace(lo[1], hi[1], grid):
                try:
                    e, g, full = calculate([x, y])
                    rows.append(
                        {"x_A": x, "y_A": y, "energy_eV": e, "status": "complete"}
                    )
                except (ValueError, RuntimeError) as exc:
                    rows.append(
                        {"x_A": x, "y_A": y, "energy_eV": None, "status": str(exc)}
                    )
        frame = pd.DataFrame(rows)
        frame["delta_energy_kj_mol"] = (
            frame.energy_eV - frame.energy_eV.min()
        ) * 96.4853321233
        frame.to_csv(output / "surface.csv", index=False)
        stationary = []
        for name, start in [
            ("ES_EI_midpoint", (anchors[0] + anchors[1]) / 2),
            ("EI_EP_midpoint", (anchors[1] + anchors[2]) / 2),
        ]:
            try:

                def gradient(q):
                    if np.any(q < lo - 0.5) or np.any(q > hi + 0.5):
                        raise ValueError(
                            "Stationary search left declared surface domain"
                        )
                    return calculate(q)[1]

                solution = root(gradient, start, options={"maxfev": 50})
                e, g, full = calculate(solution.x)
                h = numerical_hessian(gradient, solution.x, step=0.005)
                classification = classify_stationary(g, h, gradient_tolerance=0.01)
                stationary.append(
                    {
                        "candidate": name,
                        "coordinates_A": solution.x.tolist(),
                        "energy_eV": e,
                        "solver_success": bool(solution.success),
                        "full_frozen_cluster_gradient_max_eV_per_A": full,
                        **classification,
                    }
                )
            except (ValueError, RuntimeError) as exc:
                stationary.append(
                    {"candidate": name, "status": "unresolved", "error": str(exc)}
                )
        write_json(output / "mechanism.json", KSI_MECHANISM)
        journal.result(
            "H5-projected-proton-surface",
            {
                "preparation": preparation,
                "evaluations": evaluations,
                "successful_grid_points": int(frame.energy_eV.notna().sum()),
                "stationary_candidates": stationary,
                "full_TS_validated": False,
                "kcat_computed": False,
                "limitations": [
                    "Modelled substrate placement from bound ligand template; binding pose not experimentally validated",
                    "Capped cluster with ALPB water, no full protein QM/MM embedding or ensemble",
                    "All coordinates except two proton-position coordinates frozen; subspace saddles are not full 3N transition states",
                    "Electronic potential surface is not a potential of mean force or activation free energy",
                    "No KIE, kcat or chemical barrier accepted from this exploratory surface",
                ],
            },
        )
