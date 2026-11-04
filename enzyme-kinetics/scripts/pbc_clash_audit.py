"""Orthorhombic minimum-image intermolecular contacts (not a full clash score)."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def minimum_image(delta, box):
    delta = np.asarray(delta, dtype=float)
    box = np.asarray(box, dtype=float)
    if box.shape != (3,) or not np.isfinite(box).all() or (box <= 0).any():
        raise ValueError("positive orthorhombic lengths required")
    shift = -np.rint(delta / box).astype(int)
    return delta + shift * box, shift


def classify_pair(a, b, ei, ej, same_molecule, bonded):
    if same_molecule or bonded:
        return None
    groups = frozenset((a, b))
    if groups == {"ligand", "protein"}:
        return "ligand-protein"
    if groups == {"ligand", "water"}:
        return "ligand-water"
    if a == b == "protein":
        return "protein-protein"
    if a == b == "water" and not (ei == "H" and ej == "H") and (ei == "O" or ej == "O"):
        return "water-water"
    return None


def components(n, bonds):
    parent = np.arange(n)
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in bonds:
        x, y = find(i), find(j)
        if x != y:
            parent[y] = x
    return np.array([find(i) for i in range(n)])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    import openmm
    from openmm import app, unit
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("input-npy", "system-xml", "topology-pdb", "state-xml", "box-npy", "restraints-json", "outdir"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    if openmm.__version__ != "8.6.1":
        raise ValueError("requires OpenMM 8.6.1")
    system = openmm.XmlSerializer.deserialize(a.system_xml.read_text())
    state = openmm.XmlSerializer.deserialize(a.state_xml.read_text())
    box = np.load(a.box_npy, allow_pickle=False)
    x = np.load(a.input_npy, allow_pickle=False)
    statebox = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    statex = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-12) or not np.allclose(box, statebox, atol=1e-12, rtol=0):
        raise ValueError("box mismatch or unsupported nonorthorhombic box")
    length = np.diag(box)
    if not np.isfinite(x).all() or not np.isfinite(box).all() or (length <= 0).any() or x.shape != (system.getNumParticles(), 3) or not np.allclose(x, statex, rtol=0, atol=1e-12):
        raise ValueError("positions/box/state mismatch")
    atoms = list(app.PDBFile(str(a.topology_pdb)).topology.atoms())
    if len(atoms) != len(x):
        raise ValueError("topology atom count mismatch")
    bonds = set()
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for k in range(force.getNumBonds()):
                i, j, _, _ = force.getBondParameters(k)
                bonds.add(tuple(sorted((int(i), int(j)))))
    for k in range(system.getNumConstraints()):
        i, j, _ = system.getConstraintParameters(k)
        bonds.add(tuple(sorted((int(i), int(j)))))
    mol = components(len(atoms), bonds)
    cats = ["water" if z.residue.name in ("HOH", "WAT") else "ligand" if z.residue.name == "UNK" and z.residue.chain.id == "C" else "protein" if z.residue.chain.id in ("A", "B") else "ion" for z in atoms]
    elems = [z.element.symbol if z.element else "?" for z in atoms]
    if cats.count("ligand") != 47 or cats.count("protein") != 3760 or "?" in elems:
        raise ValueError("unexpected atom classification")
    wrapped = np.mod(x, length)
    candidates = cKDTree(wrapped, boxsize=length).query_pairs(0.2, output_type="ndarray")
    rows = []
    for i, j in candidates:
        i, j = int(i), int(j)
        category = classify_pair(cats[i], cats[j], elems[i], elems[j], mol[i] == mol[j], (i, j) in bonds)
        if category is None:
            continue
        delta, shift = minimum_image(x[j] - x[i], length)
        distance = float(np.linalg.norm(delta) * 10)
        if distance >= 2.0:
            continue
        ai, aj = atoms[i], atoms[j]
        hydrogen_class = "H-H" if elems[i] == elems[j] == "H" else "H-heavy" if "H" in (elems[i], elems[j]) else "heavy-heavy"
        rows.append({"category": category, "atom_i": i, "name_i": ai.name, "res_i": ai.residue.name + ai.residue.id, "chain_i": ai.residue.chain.id, "atom_j": j, "name_j": aj.name, "res_j": aj.residue.name + aj.residue.id, "chain_j": aj.residue.chain.id, "distance_angstrom": distance, "image_shift": json.dumps(shift.tolist()), "is_hard": distance < 1.5, "hydrogen_class": hydrogen_class})
    rows.sort(key=lambda row: row["distance_angstrom"])
    fields = ["category", "atom_i", "name_i", "res_i", "chain_i", "atom_j", "name_j", "res_j", "chain_j", "distance_angstrom", "image_shift", "is_hard", "hydrogen_class"]
    a.outdir.mkdir(parents=True)
    for filename, selected in (("pairs.csv", rows), ("top50.csv", rows[:50])):
        with (a.outdir / filename).open("x", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(selected)
    crossings = set()
    for i, j in bonds:
        if np.any(np.floor(x[i] / length) != np.floor(x[j] / length)):
            crossings.add(int(mol[i]))
    categories = {}
    for key in ("ligand-protein", "ligand-water", "water-water", "protein-protein"):
        selected = [r for r in rows if r["category"] == key]
        categories[key] = {"under_2A": len(selected), "under_1_5A": sum(r["is_hard"] for r in selected), "min_distance_A": min((r["distance_angstrom"] for r in selected), default=None), "max_reported_distance_A": max((r["distance_angstrom"] for r in selected), default=None), "by_hydrogen_class": {h: sum(r["hydrogen_class"] == h for r in selected) for h in ("heavy-heavy", "H-heavy", "H-H")}}
    inputs = (a.input_npy, a.system_xml, a.topology_pdb, a.state_xml, a.box_npy, a.restraints_json, Path(__file__))
    used = {"sha256": {str(z): sha(z) for z in inputs}, "openmm_version": openmm.__version__, "box_nm": box.tolist(), "box_angstrom": (box * 10).tolist(), "note": "PDB used only for labels; coordinates from NPY verified against State. Orthorhombic minimum image; no new molecular copies."}
    (a.outdir / "state-used.json").write_text(json.dumps(used, indent=2) + "\n")
    report = {"categories": categories, "total_under_2A": len(rows), "total_under_1_5A": sum(r["is_hard"] for r in rows), "molecules_crossing_wrapped_boundary": len(crossings), "limitations": ["only specified intermolecular categories", "distance thresholds are flags, not chemically validated clash definitions", "max/min only for reported pairs under 2 A", "molecular identity inferred from harmonic bonds and distance constraints", "topology atom ordering assumed to match system"]}
    (a.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
