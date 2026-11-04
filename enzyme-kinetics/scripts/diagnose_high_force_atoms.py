"""Audit projected forces by topology; coordinates are always loaded from NPY."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit

from constraint_projection import project_forces


def category(atom):
    if atom.residue.name in {"HOH", "WAT"}:
        return "water"
    if atom.residue.chain.id == "C" and atom.residue.name == "UNK":
        return "ligand"
    if atom.residue.chain.id in {"A", "B"}:
        return "protein"
    return "ion"


def summarize(norms, categories):
    result = {}
    high = norms > 50
    for name in ("protein", "ligand", "water", "ion"):
        mask = categories == name
        values = norms[mask]
        result[name] = {
            "count": int(mask.sum()),
            "max": float(values.max()) if len(values) else None,
            "rms_atom_norm": float(np.sqrt(np.mean(values**2)))
            if len(values)
            else None,
            "above_50": int((mask & high).sum()),
            "fraction_of_all_above_50": float((mask & high).sum() / high.sum())
            if high.any()
            else None,
        }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "input-npy",
        "system-xml",
        "topology-pdb",
        "restraints-json",
        "outdir",
    ):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    if openmm.__version__ != "8.6.1":
        raise RuntimeError("OpenMM 8.6.1 required")
    restraints = json.loads(a.restraints_json.read_text())
    if restraints["soft_pose_restraints"]:
        raise ValueError("Unrecognised active pose restraints")
    topology = app.PDBFile(str(a.topology_pdb)).topology
    atoms = list(topology.atoms())
    system = openmm.XmlSerializer.deserialize(a.system_xml.read_text())
    x = np.load(a.input_npy, allow_pickle=False)
    if (
        x.shape != (system.getNumParticles(), 3)
        or len(atoms) != len(x)
        or not np.isfinite(x).all()
    ):
        raise ValueError("coordinates, topology and system are inconsistent")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"}
    )
    context.setPositions(x * unit.nanometer)
    state = context.getState(getForces=True, getEnergy=True, getPositions=True)
    f = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
    )
    pairs = np.asarray(
        [
            (
                int(system.getConstraintParameters(i)[0]),
                int(system.getConstraintParameters(i)[1]),
            )
            for i in range(system.getNumConstraints())
        ]
    )
    projected, diag = project_forces(x, f, pairs)
    norms = np.linalg.norm(projected, axis=1)
    cats = np.array([category(atom) for atom in atoms])
    rows = [
        {
            "atom_index": atom.index,
            "atom_name": atom.name,
            "residue": atom.residue.name + str(atom.residue.id),
            "chain": atom.residue.chain.id,
            "category": cats[atom.index],
            "fproj_norm": float(norms[atom.index]),
        }
        for atom in atoms
    ]
    a.outdir.mkdir(parents=True)
    (a.outdir / "state.xml").write_text(openmm.XmlSerializer.serialize(state))
    np.save(
        a.outdir / "box-vectors-nm.npy",
        np.asarray(
            state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        ),
    )
    np.save(a.outdir / "positions-nm.npy", x)
    (a.outdir / "restraints.json").write_text(a.restraints_json.read_text())
    for name, data in (
        ("per-atom.csv", rows),
        ("top100.csv", sorted(rows, key=lambda r: -r["fproj_norm"])[:100]),
    ):
        with (a.outdir / name).open("x", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(data)
    inputs = (
        a.input_npy,
        a.system_xml,
        a.topology_pdb,
        a.restraints_json,
        Path(__file__),
        Path(__file__).with_name("constraint_projection.py"),
    )
    report = {
        "input_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs
        },
        "classification": "fixed-coordinate diagnostic, topology PDB coordinates not used",
        "openmm_version": openmm.__version__,
        "threads": 4,
        "projection": diag,
        "categories": summarize(norms, cats),
    }
    (a.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
