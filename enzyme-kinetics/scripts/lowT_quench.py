"""Short NVT diagnostic from exact coordinates; not equilibrated production MD."""

import argparse
import csv
import hashlib
import json
import math
import traceback
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit
from scipy.spatial import cKDTree

from constraint_projection import project_forces
from glu43_geometry import geometry, image


METRICS = [
    "time_ps",
    "potential_energy_kJ_mol",
    "kinetic_energy_kJ_mol",
    "total_energy_kJ_mol",
    "temperature_K",
    "constraint_violation_max_nm",
    "constraint_violation_relative_max",
    "ligand_rmsd_nm_vs_start",
    "protein_backbone_rmsd_nm_vs_start",
    "radius_of_gyration_nm",
    "num_ligand_protein_contacts_under_4A",
    "num_ligand_water_hbonds",
    "num_ligand_protein_hbonds",
    "projected_force_max_kJ_mol_nm",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def constraint_errors(x, pairs, targets):
    errors = np.abs(np.linalg.norm(x[pairs[:, 0]] - x[pairs[:, 1]], axis=1) - targets)
    return float(errors.max(initial=0)), float((errors / targets).max(initial=0))


def assert_valid_frame(
    x,
    forces,
    potential,
    kinetic,
    rel_error,
    tolerance,
    ligand_rmsd,
    backbone_rmsd,
    box,
    initial_box,
):
    if not all(
        np.isfinite(z).all()
        for z in (
            x,
            forces,
            np.asarray([potential, kinetic, rel_error, ligand_rmsd, backbone_rmsd]),
            box,
        )
    ):
        raise ValueError("nonfinite coordinate, force, energy or metric")
    if rel_error > tolerance:
        raise ValueError(f"relative constraint violation {rel_error:g} > {tolerance:g}")
    if ligand_rmsd >= 0.3 or backbone_rmsd >= 0.5:
        raise ValueError(
            f"RMSD threshold exceeded: ligand={ligand_rmsd:g}, backbone={backbone_rmsd:g} nm"
        )
    if not np.allclose(box, initial_box, atol=1e-12, rtol=0):
        raise ValueError("NVT box changed")


def end_checks(rows, target):
    if len(rows) < 5:
        return {
            "energy_stable": False,
            "temperature_stable": False,
            "reason": "fewer than 5 saved frames",
        }
    last = rows[-5:]
    energies = np.array([r["potential_energy_kJ_mol"] for r in last])
    temps = np.array([r["temperature_K"] for r in last])
    monotone = bool(np.all(np.diff(energies) < 0) and energies[0] - energies[-1] > 500)
    return {
        "energy_stable": not monotone,
        "energy_last5_kj_mol": energies.tolist(),
        "temperature_stable": bool(
            abs(float(temps.mean()) - target) <= 20
            and np.max(np.abs(temps - temps.mean())) <= 20
        ),
        "temperature_last5_K": temps.tolist(),
        "temperature_last5_mean_K": float(temps.mean()),
    }


def checked_blocks(blocks, metrics, target):
    """Yield one block at a time; inspect its saved metrics before the next step."""
    for block in range(blocks + 1):
        yield block
        if len(metrics) >= 5:
            checks = end_checks(metrics, target)
            # Temperature's existing gate excludes the initial sampled frame;
            # the energy safety rule includes it as a consecutive saved frame.
            temperature_failed = (
                len(metrics[1:]) >= 5
                and not end_checks(metrics[1:], target)["temperature_stable"]
            )
            if not checks["energy_stable"] or temperature_failed:
                raise ValueError(
                    "rolling energy/temperature criterion failed: " + json.dumps(checks)
                )


def aligned_coordinates(x, ref, backbone, ligand, lengths):
    # Reconstruct each atom in the image closest to its initial position; this
    # short-run convention fails loudly on the prescribed large RMSD limits.
    moved = ref + image(x - ref, lengths)
    ref_center = ref[backbone].mean(axis=0)
    cur_center = moved[backbone].mean(axis=0)
    u, _, vt = np.linalg.svd(
        (moved[backbone] - cur_center).T @ (ref[backbone] - ref_center)
    )
    fix = np.eye(3)
    fix[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ fix @ vt
    aligned = (moved - cur_center) @ rotation + ref_center
    backbone_rmsd = float(
        np.sqrt(np.mean(np.sum((aligned[backbone] - ref[backbone]) ** 2, axis=1)))
    )
    ligand_rmsd = float(
        np.sqrt(np.mean(np.sum((aligned[ligand] - ref[ligand]) ** 2, axis=1)))
    )
    return aligned, backbone_rmsd, ligand_rmsd


def hbond_count(x, lengths, ligand_acceptors, water_triplets, ligand_donors):
    # Geometric candidate: donor-acceptor <3.5 A, H-acceptor <2.5 A and
    # D-H-A >=150 degrees. Ligand N/O acceptor typing is provisional.
    count = 0
    waters = np.array([w[0] for w in water_triplets], dtype=int)
    tree = cKDTree(np.mod(x[waters], lengths), boxsize=lengths)
    for acceptor in ligand_acceptors:
        for wi in tree.query_ball_point(np.mod(x[acceptor], lengths), 0.35):
            o, *hs = water_triplets[wi]
            for h in hs:
                dh = image(x[o] - x[h], lengths)
                ah = image(x[acceptor] - x[h], lengths)
                if np.linalg.norm(ah) < 0.25 and np.dot(dh, ah) / (
                    np.linalg.norm(dh) * np.linalg.norm(ah)
                ) <= math.cos(math.radians(150)):
                    count += 1
    for donor, h in ligand_donors:
        for wi in tree.query_ball_point(np.mod(x[donor], lengths), 0.35):
            o = waters[wi]
            dh = image(x[donor] - x[h], lengths)
            ah = image(x[o] - x[h], lengths)
            if np.linalg.norm(ah) < 0.25 and np.dot(dh, ah) / (
                np.linalg.norm(dh) * np.linalg.norm(ah)
            ) <= math.cos(math.radians(150)):
                count += 1
    return count


def protein_ligand_hbonds(x, lengths, donors, acceptors, categories):
    """Count geometric N/O donor-H to other-group N/O candidates, angle >120°."""
    count = 0
    for donor, hydrogen in donors:
        for acceptor in acceptors:
            if (
                categories[donor] == categories[acceptor]
                or categories[donor] not in ("protein", "ligand")
                or categories[acceptor] not in ("protein", "ligand")
            ):
                continue
            da = image(x[acceptor] - x[donor], lengths)
            if np.linalg.norm(da) >= 0.35:
                continue
            dh = image(x[donor] - x[hydrogen], lengths)
            ah = image(x[acceptor] - x[hydrogen], lengths)
            denom = np.linalg.norm(dh) * np.linalg.norm(ah)
            if denom > 0 and np.dot(dh, ah) / denom < -0.5:
                count += 1
    return count


def add_backbone_restraint(system, indices, reference_nm, k_kcal_mol_a2):
    """Add periodic k*r² bias and return its force index and provenance."""
    x = np.asarray(reference_nm, dtype=float)
    ids = np.asarray(indices, dtype=int)
    if x.shape != (system.getNumParticles(), 3) or not np.isfinite(x).all():
        raise ValueError("invalid restraint reference")
    if (
        ids.ndim != 1
        or len(ids) == 0
        or len(np.unique(ids)) != len(ids)
        or np.any(ids < 0)
        or np.any(ids >= len(x))
    ):
        raise ValueError("invalid restrained indices")
    if not math.isfinite(k_kcal_mol_a2) or k_kcal_mol_a2 <= 0:
        raise ValueError("restraint strength must be positive and finite")
    force = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", k_kcal_mol_a2 * 418.4)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    for i in ids:
        force.addParticle(int(i), list(map(float, x[i])))
    index = system.addForce(force)
    return index, {
        "type": "backbone_position_periodic",
        "expression": "k*periodicdistance^2",
        "k_kcal_mol_A2": k_kcal_mol_a2,
        "atom_indices": ids.tolist(),
        "reference": "input NPY",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in (
        "positions-npy",
        "state-xml",
        "system-xml",
        "box-vectors",
        "restraints-json",
        "topology-pdb",
        "outdir",
    ):
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument("--temperature", type=float, required=True)
    p.add_argument("--friction", type=float, default=10.0)
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--duration-ps", type=float, default=20.0)
    p.add_argument("--save-interval-ps", type=float, default=1.0)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument(
        "--backbone-restraints-kcal",
        type=float,
        default=0.0,
        help="k in k*r^2, kcal/(mol A^2); zero disables",
    )
    a = p.parse_args()
    if openmm.__version__ != "8.6.1" or a.outdir.exists():
        raise ValueError("OpenMM 8.6.1 and a new output directory required")
    if (
        min(
            a.temperature,
            a.friction,
            a.timestep_fs,
            a.duration_ps,
            a.save_interval_ps,
            a.seed,
        )
        <= 0
        or not math.isfinite(a.backbone_restraints_kcal)
        or a.backbone_restraints_kcal < 0
    ):
        raise ValueError("parameters must be positive")
    steps_per = round(a.save_interval_ps * 1000 / a.timestep_fs)
    blocks = round(a.duration_ps / a.save_interval_ps)
    if not math.isclose(
        steps_per * a.timestep_fs, a.save_interval_ps * 1000
    ) or not math.isclose(blocks * a.save_interval_ps, a.duration_ps):
        raise ValueError("duration and interval must contain integral steps and blocks")
    system = openmm.XmlSerializer.deserialize(a.system_xml.read_text())
    if any("Barostat" in type(force).__name__ for force in system.getForces()):
        raise ValueError("barostat not allowed in NVT")
    restraints = json.loads(a.restraints_json.read_text())
    if restraints["soft_pose_restraints"]:
        raise ValueError("unsupported soft restraints")
    state_in = openmm.XmlSerializer.deserialize(a.state_xml.read_text())
    x = np.load(a.positions_npy, allow_pickle=False)
    box = np.load(a.box_vectors, allow_pickle=False)
    state_x = np.asarray(
        state_in.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    )
    state_box = np.asarray(
        state_in.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    )
    if (
        x.shape != (system.getNumParticles(), 3)
        or not np.isfinite(x).all()
        or box.shape != (3, 3)
        or not np.allclose(box, np.diag(np.diag(box)))
        or not np.allclose(x, state_x, rtol=0, atol=1e-12)
        or not np.allclose(box, state_box, rtol=0, atol=1e-12)
    ):
        raise ValueError("NPY, State or box mismatch")
    lengths = np.diag(box)
    if np.any(lengths <= 0):
        raise ValueError("invalid box")
    top = app.PDBFile(str(a.topology_pdb)).topology
    atoms = list(top.atoms())
    if len(atoms) != len(x):
        raise ValueError("topology size mismatch")
    ligand = np.array(
        [
            z.index
            for z in atoms
            if z.residue.chain.id == "C"
            and z.residue.name == "UNK"
            and z.element.symbol != "H"
        ]
    )
    protein = np.array(
        [
            z.index
            for z in atoms
            if z.residue.chain.id in ("A", "B") and z.element.symbol != "H"
        ]
    )
    backbone = np.array(
        [
            z.index
            for z in atoms
            if z.residue.chain.id in ("A", "B") and z.name in ("N", "CA", "C", "O")
        ]
    )
    glu = {
        z.name: z.index
        for z in atoms
        if z.residue.chain.id == "A"
        and z.residue.id == "43"
        and z.residue.name == "GLU"
    }
    if (
        len(ligand) != 21
        or not len(backbone)
        or glu.get("C") != 601
        or not {"C", "O", "CG", "CD", "OE1", "OE2"} <= glu.keys()
    ):
        raise ValueError("unexpected topology or GLU43 mapping")
    if a.backbone_restraints_kcal:
        _, info = add_backbone_restraint(
            system, backbone, x, a.backbone_restraints_kcal
        )
        restraints["soft_pose_restraints"].append(info)
    residues = list(top.residues())
    water_triplets = []
    for residue in residues:
        if residue.name == "HOH":
            water = {z.name: z.index for z in residue.atoms()}
            if not {"O", "H1", "H2"} <= water.keys():
                raise ValueError("unrecognised water atom names")
            water_triplets.append((water["O"], water["H1"], water["H2"]))
    acceptors = [i for i in ligand if atoms[i].element.symbol in ("N", "O")]
    lig_donors = []
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for i in range(force.getNumBonds()):
                ai, aj, _, _ = force.getBondParameters(i)
                ai, aj = int(ai), int(aj)
                if ai in acceptors and atoms[aj].element.symbol == "H":
                    lig_donors.append((ai, aj))
                if aj in acceptors and atoms[ai].element.symbol == "H":
                    lig_donors.append((aj, ai))
    pairs, targets = [], []
    for i in range(system.getNumConstraints()):
        j, k, d = system.getConstraintParameters(i)
        pairs.append((int(j), int(k)))
        targets.append(d.value_in_unit(unit.nanometer))
    pairs, targets = np.asarray(pairs, dtype=int), np.asarray(targets)
    if len(pairs) != restraints["hard_distance_constraints"]["count"]:
        raise ValueError("constraint count differs from provenance")
    categories = [
        "protein"
        if z.residue.chain.id in ("A", "B")
        else "ligand"
        if z.residue.chain.id == "C" and z.residue.name == "UNK"
        else "other"
        for z in atoms
    ]
    acceptor_candidates = [
        z.index
        for z in atoms
        if categories[z.index] in ("protein", "ligand")
        and z.element.symbol in ("N", "O")
    ]
    donor_pairs = set()
    all_bonds = [tuple(pair) for pair in pairs]
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            all_bonds.extend(
                (int(force.getBondParameters(i)[0]), int(force.getBondParameters(i)[1]))
                for i in range(force.getNumBonds())
            )
    for i, j in all_bonds:
        for donor, hydrogen in ((i, j), (j, i)):
            if (
                categories[donor] in ("protein", "ligand")
                and atoms[donor].element.symbol in ("N", "O")
                and atoms[hydrogen].element.symbol == "H"
            ):
                donor_pairs.add((donor, hydrogen))
    integrator = openmm.LangevinMiddleIntegrator(
        a.temperature * unit.kelvin,
        a.friction / unit.picosecond,
        a.timestep_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(a.seed)
    tolerance = integrator.getConstraintTolerance()
    absolute, relative = constraint_errors(x, pairs, targets)
    if relative > tolerance:
        raise ValueError(
            f"initial relative constraint violation {relative} > {tolerance}"
        )
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"}
    )
    context.setPeriodicBoxVectors(*[openmm.Vec3(*v) * unit.nanometer for v in box])
    context.setPositions(x * unit.nanometer)
    context.setVelocitiesToTemperature(a.temperature * unit.kelvin, a.seed)
    context.applyVelocityConstraints(tolerance)
    masses = np.array(
        [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(len(x))]
    )
    dof = (
        3 * int(np.count_nonzero(masses))
        - len(pairs)
        - (
            3
            if any(isinstance(f, openmm.CMMotionRemover) for f in system.getForces())
            else 0
        )
    )
    inputs = [
        a.positions_npy,
        a.state_xml,
        a.system_xml,
        a.box_vectors,
        a.restraints_json,
        a.topology_pdb,
        Path(__file__),
        Path(__file__).with_name("constraint_projection.py"),
        Path(__file__).with_name("glu43_geometry.py"),
    ]
    config = {
        "label": "short NVT MD test at 298.15 K"
        if a.temperature == 298.15
        else "short NVT MD diagnostic",
        "classification": "restart from coordinates, not bit-exact continuation; velocities newly sampled",
        "temperature_K": a.temperature,
        "friction_per_ps": a.friction,
        "timestep_fs": a.timestep_fs,
        "duration_ps": a.duration_ps,
        "interval_ps": a.save_interval_ps,
        "seed": a.seed,
        "platform": "CPU",
        "threads": 4,
        "openmm": openmm.__version__,
        "constraints": "as serialized; HBonds preparation",
        "constraints_count": len(pairs),
        "constraint_tolerance_relative": tolerance,
        "restraints": restraints,
        "input_sha256": {str(path): sha(path) for path in inputs},
        "hbond_definition": "geometric candidates, N/O ligand acceptors provisional; D-A <3.5 A, H-A <2.5 A, D-H-A >=150 deg",
        "rmsd_definition": "minimum-image to initial, protein backbone Kabsch; ligand measured in protein-aligned frame",
        "energy_check": "last five frames strictly monotonic decline >500 kJ/mol fails",
        "temperature_check": "last five frames mean within 20 K of target and each within 20 K of mean",
    }
    a.outdir.mkdir(parents=True)
    # JSON is a valid YAML 1.2 document.
    (a.outdir / "run-config.yaml").write_text(json.dumps(config, indent=2) + "\n")
    (a.outdir / "restraints.json").write_text(json.dumps(restraints, indent=2) + "\n")
    metrics, glu_rows = [], []
    error = None
    try:
        with (
            (a.outdir / "trajectory.dcd").open("wb") as traj,
            (a.outdir / "metrics.csv").open("w", newline="") as mf,
            (a.outdir / "glu43-tracking.csv").open("w", newline="") as gf,
        ):
            dcd = app.DCDFile(
                traj,
                top,
                a.timestep_fs * unit.femtoseconds,
                firstStep=steps_per,
                interval=steps_per,
            )
            mw = csv.DictWriter(mf, fieldnames=METRICS)
            mw.writeheader()
            gw = csv.DictWriter(
                gf,
                fieldnames=[
                    "time_ps",
                    "backbone_C_projected_force_kJ_mol_nm",
                    "backbone_C_O_nm",
                    "sidechain_CD_OE1_nm",
                    "sidechain_CD_OE2_nm",
                    "sidechain_OE1_CD_OE2_deg",
                    "sidechain_CG_CD_OE1_OE2_improper_deg",
                    "backbone_C_neighbors_under_3_5A",
                ],
            )
            gw.writeheader()
            for block in checked_blocks(blocks, metrics, a.temperature):
                if block:
                    integrator.step(steps_per)
                state = context.getState(
                    getEnergy=True,
                    getForces=True,
                    getPositions=True,
                    getVelocities=True,
                    enforcePeriodicBox=False,
                )
                coords = np.asarray(
                    state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                )
                forces = np.asarray(
                    state.getForces(asNumpy=True).value_in_unit(
                        unit.kilojoule_per_mole / unit.nanometer
                    )
                )
                current_box = np.asarray(
                    state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
                        unit.nanometer
                    )
                )
                pot = float(
                    state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                )
                kin = float(
                    state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                )
                # Persist the state before any validation so failures remain inspectable.
                label = f"{block:02d}"
                (a.outdir / f"state-block-{label}.xml").write_text(
                    openmm.XmlSerializer.serialize(state)
                )
                np.save(a.outdir / f"box-block-{label}.npy", current_box)
                np.save(a.outdir / f"positions-block-{label}.npy", coords)
                (a.outdir / f"restraints-block-{label}.json").write_text(
                    json.dumps(restraints, indent=2) + "\n"
                )
                absolute, relative = constraint_errors(coords, pairs, targets)
                aligned, brms, lrms = aligned_coordinates(
                    coords, x, backbone, ligand, lengths
                )
                assert_valid_frame(
                    coords,
                    forces,
                    pot,
                    kin,
                    relative,
                    tolerance,
                    lrms,
                    brms,
                    current_box,
                    box,
                )
                if absolute >= 1e-4:
                    raise ValueError(
                        f"absolute constraint error {absolute:g} nm exceeds 1e-4 nm"
                    )
                proj, _ = project_forces(coords, forces, pairs)
                norms = np.linalg.norm(proj, axis=1)
                tree = cKDTree(np.mod(coords[protein], lengths), boxsize=lengths)
                contacts = sum(
                    len(v)
                    for v in tree.query_ball_point(np.mod(coords[ligand], lengths), 0.4)
                )
                rg = float(
                    np.sqrt(
                        np.mean(
                            np.sum(
                                (aligned[backbone] - aligned[backbone].mean(axis=0))
                                ** 2,
                                axis=1,
                            )
                        )
                    )
                )
                row = {
                    "time_ps": block * a.save_interval_ps,
                    "potential_energy_kJ_mol": pot,
                    "kinetic_energy_kJ_mol": kin,
                    "total_energy_kJ_mol": pot + kin,
                    "temperature_K": 2 * kin / (dof * 0.00831446261815324),
                    "constraint_violation_max_nm": absolute,
                    "constraint_violation_relative_max": relative,
                    "ligand_rmsd_nm_vs_start": lrms,
                    "protein_backbone_rmsd_nm_vs_start": brms,
                    "radius_of_gyration_nm": rg,
                    "num_ligand_protein_contacts_under_4A": contacts,
                    "num_ligand_water_hbonds": hbond_count(
                        coords, lengths, acceptors, water_triplets, lig_donors
                    ),
                    "num_ligand_protein_hbonds": protein_ligand_hbonds(
                        coords, lengths, donor_pairs, acceptor_candidates, categories
                    ),
                    "projected_force_max_kJ_mol_nm": float(norms.max()),
                }
                g = geometry({name: coords[i] for name, i in glu.items()}, lengths)
                displacement = image(coords - coords[glu["C"]], lengths)
                nearest = np.where(
                    (np.linalg.norm(displacement, axis=1) < 0.35)
                    & (np.arange(len(coords)) != glu["C"])
                )[0]
                g.update(
                    {
                        "time_ps": row["time_ps"],
                        "backbone_C_projected_force_kJ_mol_nm": float(norms[glu["C"]]),
                        "backbone_C_neighbors_under_3_5A": json.dumps(
                            [
                                {
                                    "index": int(i),
                                    "distance_A": float(
                                        np.linalg.norm(displacement[i]) * 10
                                    ),
                                }
                                for i in nearest
                            ]
                        ),
                    }
                )
                mw.writerow(row)
                mf.flush()
                gw.writerow(g)
                gf.flush()
                metrics.append(row)
                glu_rows.append(g)
                if block:
                    dcd.writeModel(
                        state.getPositions(),
                        periodicBoxVectors=state.getPeriodicBoxVectors(),
                    )
                    traj.flush()
            checks = end_checks(metrics[1:], a.temperature)
            if not checks["energy_stable"] or not checks["temperature_stable"]:
                raise ValueError(
                    "final energy/temperature criterion failed: " + json.dumps(checks)
                )
    except Exception as exc:
        error = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    report = {
        "input_sha256": config["input_sha256"],
        "run_config_sha256": sha(a.outdir / "run-config.yaml"),
        "status": "failed" if error else "completed",
        "error": error,
        "saved_metric_frames": len(metrics),
        "end_checks": end_checks(metrics[1:], a.temperature),
        "initial": metrics[0] if metrics else None,
        "final": metrics[-1] if metrics else None,
        "glu43_initial": glu_rows[0] if glu_rows else None,
        "glu43_final": glu_rows[-1] if glu_rows else None,
        "limitations": [
            "short diagnostic is not equilibrated production MD",
            "force projection not a convergence criterion",
            "PDB only supplies topology; atom ordering assumed",
            "fixed-coordinate image reconstruction assumes short displacements",
            "ligand acceptor types and hbond counts provisional",
            "no kinetic or binding free-energy interpretation",
        ],
    }
    (a.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
