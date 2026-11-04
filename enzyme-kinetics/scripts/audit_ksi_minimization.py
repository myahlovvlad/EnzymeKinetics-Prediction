"""Audit and continue a candidate KSI minimization without altering its inputs.

Force statistics are raw Cartesian forces, not constraint-projected gradients;
two finite iteration budgets do not by themselves prove convergence.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def force_stats(forces):
    values = np.asarray(forces, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3 or not np.isfinite(values).all():
        raise ValueError("Expected finite N x 3 forces")
    return {"rms_component_kj_mol_nm": float(np.sqrt(np.mean(values**2))),
            "max_atom_norm_kj_mol_nm": float(np.linalg.norm(values, axis=1).max())}


def reaction_geometry(xyz, c4, h, oxygens):
    rows = []
    for name, index in oxygens:
        v1 = xyz[c4] - xyz[h]
        v2 = xyz[index] - xyz[h]
        norm = np.linalg.norm(v1) * np.linalg.norm(v2)
        if norm <= 0:
            raise ValueError("Degenerate reaction geometry")
        rows.append({"oxygen": name, "c4_o_A": float(np.linalg.norm(xyz[c4] - xyz[index]) * 10),
                     "h_o_A": float(np.linalg.norm(v2) * 10),
                     "angle_deg": float(np.degrees(np.arccos(np.clip(v1 @ v2 / norm, -1, 1))))})
    return rows


def main():
    import openmm
    from openmm import app, unit

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--system", type=Path, required=True)
    p.add_argument("--assignment", type=Path, required=True)
    p.add_argument("--positions-npy", type=Path, help="exact N x 3 coordinates in nm from a previous audit of the same system")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--iterations", type=int, default=300)
    p.add_argument("--tolerance", type=float, default=1.0, help="OpenMM kJ/mol/nm")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    if args.iterations <= 0 or args.tolerance <= 0 or args.threads <= 0:
        p.error("Iterations, tolerance and threads must be positive")
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError("Choose a new output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    pdb = app.PDBFile(str(args.pdb))
    system = openmm.XmlSerializer.deserialize(args.system.read_text())
    atoms = list(pdb.topology.atoms())
    if system.getNumParticles() != len(atoms):
        raise ValueError("System/topology atom count mismatch")
    ligand = [a.index for a in atoms if a.residue.name == "UNK" and a.residue.chain.id == "C"]
    protein = [a.index for a in atoms if a.residue.chain.id in ("A", "B") and a.element is not None and a.element.symbol not in ("H", "D")]
    ligand_heavy = [i for i in ligand if atoms[i].element.symbol not in ("H", "D")]
    waters = [a.index for a in atoms if a.residue.name == "HOH" and a.name == "O"]
    crystal = [a.index for a in atoms if a.residue.name == "HOH" and a.residue.chain.id == "D" and a.name == "O"]
    if len(ligand) != 47 or len(crystal) != 48 or not protein or not waters:
        raise ValueError("Unexpected ligand, protein or crystal-water atom counts")
    assignment = json.loads(args.assignment.read_text())
    c4 = ligand[assignment["atom_names_zero_based"]["C4"]]
    h = ligand[assignment["poses"][0]["beta_h_index"]]
    oxygens = [(a.name, a.index) for a in atoms if a.residue.chain.id == "A" and a.residue.id == "38" and a.residue.name == "ASP" and a.name in ("OD1", "OD2")]
    if len(oxygens) != 2 or atoms[h].element.symbol != "H":
        raise ValueError("Unexpected reaction atom mapping")
    nonbonded = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(nonbonded) != 1:
        raise ValueError("Expected one NonbondedForce")
    params = [nonbonded[0].getParticleParameters(i) for i in range(len(atoms))]
    q = np.array([v[0].value_in_unit(unit.elementary_charge) for v in params])
    sigma = np.array([v[1].value_in_unit(unit.nanometer) for v in params])
    epsilon = np.array([v[2].value_in_unit(unit.kilojoule_per_mole) for v in params])
    if not np.isfinite(q).all() or not np.isfinite(sigma).all() or not np.isfinite(epsilon).all():
        raise ValueError("Nonfinite nonbonded parameters")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"), {"Threads": str(args.threads)})
    if args.positions_npy:
        coords = np.load(args.positions_npy, allow_pickle=False)
        if coords.shape != (len(atoms), 3) or not np.isfinite(coords).all():
            raise ValueError("Invalid exact restart positions")
        context.setPositions(coords * unit.nanometer)
    else:
        context.setPositions(pdb.positions)
    def snapshot(label):
        state = context.getState(getEnergy=True, getForces=True, getPositions=True)
        coords = np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        forces = np.array(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        pl = cKDTree(coords[protein]).query(coords[ligand_heavy], k=1)[0]
        wl = cKDTree(coords[waters]).query(coords[ligand_heavy], k=1)[0]
        max_constraint_error = max((abs(np.linalg.norm(coords[int(system.getConstraintParameters(i)[0])] - coords[int(system.getConstraintParameters(i)[1])]) - system.getConstraintParameters(i)[2].value_in_unit(unit.nanometer)) for i in range(system.getNumConstraints())), default=0.)
        data = {"label": label, "potential_kj_mol": float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),
                "force_all": force_stats(forces), "force_ligand": force_stats(forces[ligand]),
                "reaction_geometry": reaction_geometry(coords, c4, h, oxygens),
                "min_ligand_protein_heavy_A": float(pl.min() * 10),
                "min_ligand_water_oxygen_A": float(wl.min() * 10),
                "ligand_heavy_atoms_with_protein_within_2A": int(np.sum(pl < .2)),
                "max_constraint_distance_error_A": float(max_constraint_error * 10)}
        if not np.isfinite(data["potential_kj_mol"]):
            raise ValueError("Nonfinite potential energy")
        np.save(args.output / f"{label}-positions-nm.npy", coords)
        with (args.output / f"{label}.pdb").open("w") as f:
            app.PDBFile.writeFile(pdb.topology, state.getPositions(), f, keepIds=False)
        return data, coords
    snapshots = []
    initial, ref = snapshot("initial")
    snapshots.append(initial)
    for stage in (1, 2):
        started = time.monotonic()
        openmm.LocalEnergyMinimizer.minimize(context, tolerance=args.tolerance * unit.kilojoule_per_mole / unit.nanometer, maxIterations=args.iterations)
        current, coords = snapshot(f"stage-{stage}")
        current["wall_seconds"] = time.monotonic() - started
        displacement = np.linalg.norm(coords[crystal] - ref[crystal], axis=1) * 10
        current["crystal_water_displacement_A_median"] = float(np.median(displacement))
        current["crystal_water_displacement_A_max"] = float(displacement.max())
        snapshots.append(current)
    inputs = [args.pdb, args.system, args.assignment] + ([args.positions_npy] if args.positions_npy else [])
    report = {"inputs_sha256": {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
              "openmm_version": openmm.__version__, "platform": "CPU", "threads": args.threads,
              "iterations_per_stage": args.iterations, "tolerance_kj_mol_nm": args.tolerance,
              "particle_count": len(atoms), "constraints": system.getNumConstraints(),
              "total_charge_e": float(q.sum()), "ligand_charge_e": float(q[ligand].sum()),
              "nonfinite_parameters": 0, "negative_sigma_or_epsilon_count": int(np.count_nonzero((sigma < 0) | (epsilon < 0))),
              "crystal_waters": len(crystal), "snapshots": snapshots,
              "limitations": ["Exact NPY restart; caller must ensure it belongs to this system" if args.positions_npy else "PDB coordinates rounded to 0.001 Å at restart", "Raw force statistics include constraint forces; not minimizer convergence criterion", "Finite iterations and similar energies do not establish equilibrium", "No proof of productive ES, correct protonation or thermodynamic stability", "Clash distances are nonperiodic; this is not a complete all-atom clash audit"]}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
