"""Conditional side-chain potential scans with explicit coupling double differences."""

from pathlib import Path
import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import digest_file
from .journal import Journal


def rotating_atoms(topology, chain, resid):
    residues = [
        r for r in topology.residues() if r.chain.id == chain and r.id == str(resid)
    ]
    if len(residues) != 1:
        raise ValueError("Ambiguous residue selection")
    residue = residues[0]
    names = {a.name: a for a in residue.atoms()}
    if not {"CA", "CB"} <= names.keys():
        raise ValueError("chi1 axis unavailable")
    a, b = names["CA"].index, names["CB"].index
    graph = {x.index: set() for x in topology.atoms()}
    for left, right in topology.bonds():
        if {left.index, right.index} == {a, b}:
            continue
        graph[left.index].add(right.index)
        graph[right.index].add(left.index)
    visited = set()
    pending = [b]
    while pending:
        node = pending.pop()
        if node in visited:
            continue
        visited.add(node)
        pending.extend(graph[node] - visited)
    if a in visited:
        raise ValueError("Ring or crosslink prevents a simple side-chain rotation")
    if any(x.index in visited and x.residue != residue for x in topology.atoms()):
        raise ValueError("Rotation would move atoms outside chosen residue")
    return a, b, sorted(visited)


def rotate(xyz, selection, degrees):
    a, b, indices = selection
    out = np.asarray(xyz, float).copy()
    axis = out[b] - out[a]
    norm = np.linalg.norm(axis)
    if norm < 1e-10:
        raise ValueError("Zero rotation axis")
    axis /= norm
    angle = np.radians(degrees)
    x = out[indices] - out[b]
    out[indices] = (
        x * np.cos(angle)
        + np.cross(axis, x) * np.sin(angle)
        + np.outer(x @ axis, axis) * (1 - np.cos(angle))
        + out[b]
    )
    return out


def run(prepared, output, threads=24):
    import openmm as mm
    from openmm import app, unit

    prepared, output = Path(prepared), Path(output)
    journal = Journal(output)
    pdb = app.PDBFile(str(prepared / "topology.pdb"))
    system = mm.XmlSerializer.deserialize((prepared / "system.xml").read_text())
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("CPU"),
        {"Threads": str(threads)},
    )

    def energy(positions):
        context.setPositions(positions)
        value = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        if not np.isfinite(value):
            raise ValueError("Nonfinite conditional potential")
        return value

    with journal.run(
        "H4-conditional-torsions",
        {
            "threads": threads,
            "angles_deg": list(range(-180, 181, 30)),
            "frozen_environment": True,
        },
        {"system_sha256": digest_file(prepared / "system.xml")},
    ):
        reference = energy(xyz)
        rows = []
        for resid in [14, 38, 54, 99, 116]:
            selection = rotating_atoms(pdb.topology, "A", resid)
            for angle in range(-180, 181, 30):
                rows.append(
                    {
                        "resid": resid,
                        "rotation_deg": angle,
                        "delta_potential_kj_mol": energy(rotate(xyz, selection, angle))
                        - reference,
                    }
                )
        pd.DataFrame(rows).to_csv(output / "torsions.csv", index=False)
        a = rotating_atoms(pdb.topology, "A", 14)
        b = rotating_atoms(pdb.topology, "A", 99)
        couplings = []
        for x in [-30, 0, 30]:
            for y in [-30, 0, 30]:
                ea = energy(rotate(xyz, a, x))
                eb = energy(rotate(xyz, b, y))
                eab = energy(rotate(rotate(xyz, a, x), b, y))
                couplings.append(
                    {
                        "Tyr14_rotation_deg": x,
                        "Asp99_rotation_deg": y,
                        "coupling_kj_mol": eab - ea - eb + reference,
                        "delta_potential_kj_mol": eab - reference,
                    }
                )
        pd.DataFrame(couplings).to_csv(output / "torsion-coupling.csv", index=False)
        journal.result(
            "H4-conditional-torsions",
            {
                "max_abs_coupling_kj_mol": max(
                    abs(r["coupling_kj_mol"]) for r in couplings
                ),
                "reference_energy_kj_mol": reference,
                "limitations": [
                    "Conditional fixed-environment potential; no relaxation, not PMF or free energy",
                    "Large clash energies are retained and indicate inaccessible frozen configurations, not physical barriers",
                    "Only side-chain chi1 scanned; biochemical mutations and constant-pH coupling are not simulated",
                ],
            },
        )
