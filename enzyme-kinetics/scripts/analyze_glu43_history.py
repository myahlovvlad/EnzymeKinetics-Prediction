"""Audit saved KSI NVT states; frames after 6 ps violated the stop protocol."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit
from scipy.spatial import cKDTree

from constraint_projection import project_forces
from glu43_geometry import image


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dihedral(p0, p1, p2, p3, box):
    """Signed torsion in degrees, using minimum-image consecutive bonds."""
    b0 = -image(np.asarray(p1) - p0, box)
    b1 = image(np.asarray(p2) - p1, box)
    b2 = image(np.asarray(p3) - p2, box)
    n = np.linalg.norm(b1)
    if n < 1e-12:
        raise ValueError("degenerate torsion")
    b1 = b1 / n
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-12:
        raise ValueError("degenerate torsion planes")
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def selected_indices(top):
    mapping = {}
    for res in top.residues():
        if res.chain.id in ("A", "B") and res.id in ("42", "43", "44"):
            key = (res.chain.id, int(res.id))
            if key in mapping:
                raise ValueError("non-unique residue identifier")
            mapping[key] = {atom.name: atom.index for atom in res.atoms()}
    for chain in ("A", "B"):
        for resid in (42, 43, 44):
            if not {"N", "CA", "C", "O"} <= mapping.get((chain, resid), {}).keys():
                raise ValueError(f"missing backbone {chain}:{resid}")
    return mapping


def frame_geometry(x, box, mapping, chain):
    prev, cur, nxt = (mapping[(chain, i)] for i in (42, 43, 44))
    xyz = lambda group, name: x[group[name]]
    torsion = lambda *args: dihedral(*args, box)
    dist = lambda a, b: float(np.linalg.norm(image(a-b, box))*10)
    return {"omega_SER42_GLU43_deg": torsion(xyz(prev,"CA"),xyz(prev,"C"),xyz(cur,"N"),xyz(cur,"CA")),
            "omega_GLU43_PRO44_deg": torsion(xyz(cur,"CA"),xyz(cur,"C"),xyz(nxt,"N"),xyz(nxt,"CA")),
            "phi_GLU43_deg": torsion(xyz(prev,"C"),xyz(cur,"N"),xyz(cur,"CA"),xyz(cur,"C")),
            "psi_GLU43_deg": torsion(xyz(cur,"N"),xyz(cur,"CA"),xyz(cur,"C"),xyz(nxt,"N")),
            "C43_O43_A": dist(xyz(cur,"C"),xyz(cur,"O")),
            "CA43_C43_A": dist(xyz(cur,"CA"),xyz(cur,"C")),
            "C43_N44_A": dist(xyz(cur,"C"),xyz(nxt,"N")),
            "C42_N43_A": dist(xyz(prev,"C"),xyz(cur,"N"))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("frames", "system-xml", "topology-pdb", "outdir"):
        p.add_argument("--" + name, required=True, type=Path)
    a = p.parse_args()
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    if openmm.__version__ != "8.6.1":
        raise RuntimeError("OpenMM 8.6.1 required")
    system = openmm.XmlSerializer.deserialize(a.system_xml.read_text())
    topology = app.PDBFile(str(a.topology_pdb)).topology
    atoms = list(topology.atoms())
    mapping = selected_indices(topology)
    if len(atoms) != system.getNumParticles() or mapping[("A", 43)]["C"] != 601:
        raise ValueError("atom count or mapping mismatch")
    constraints = [system.getConstraintParameters(i) for i in range(system.getNumConstraints())]
    pairs = np.array([(int(i), int(j)) for i, j, _ in constraints])
    # Constraints include water H–H; using them as 1–2 exclusions is deliberate.
    bonds = {tuple(sorted(pair)) for pair in pairs.tolist()}
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            bonds.update(tuple(sorted(map(int, force.getBondParameters(i)[:2]))) for i in range(force.getNumBonds()))
    adjacency = [set() for _ in atoms]
    for i, j in bonds:
        adjacency[i].add(j); adjacency[j].add(i)
    geometry_rows, environment_rows, local_rows, force_rows = [], [], [], []
    inputs = [a.system_xml, a.topology_pdb, Path(__file__), Path(__file__).with_name("constraint_projection.py"), Path(__file__).with_name("glu43_geometry.py")]
    for frame in range(21):
        state_path = a.frames / f"state-block-{frame:02d}.xml"
        coord_path = a.frames / f"positions-block-{frame:02d}.npy"
        box_path = a.frames / f"box-block-{frame:02d}.npy"
        inputs.extend((state_path, coord_path, box_path))
        state = openmm.XmlSerializer.deserialize(state_path.read_text())
        x = np.load(coord_path, allow_pickle=False)
        box_matrix = np.load(box_path, allow_pickle=False)
        sx = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        sb = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
        f = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
        if x.shape != (len(atoms), 3) or not np.isfinite(x).all() or not np.allclose(x, sx, rtol=0, atol=1e-12) or box_matrix.shape != (3, 3) or not np.allclose(box_matrix, sb, rtol=0, atol=1e-12) or not np.allclose(box_matrix, np.diag(np.diag(box_matrix)), rtol=0, atol=1e-12):
            raise ValueError(f"invalid frame {frame}")
        box = np.diag(box_matrix)
        if np.any(box <= 0) or not np.isfinite(f).all():
            raise ValueError(f"invalid box or forces: {frame}")
        d = x[pairs[:, 0]] - x[pairs[:, 1]]
        if not np.allclose(d, image(d, box), rtol=0, atol=1e-10):
            raise ValueError(f"constraint spans periodic image in frame {frame}; unwrap before projection")
        projected, diagnostic = project_forces(x, f, pairs)
        common = {"frame": frame, "time_ps": float(frame), "post_violation": frame > 6}
        for chain in ("A", "B"):
            geometry_rows.append({**common, "chain": chain, **frame_geometry(x, box, mapping, chain)})
            atom = mapping[(chain, 43)]["C"]
            force_rows.append({**common, "chain": chain, "atom_index": atom, "projected_force_kj_mol_nm": float(np.linalg.norm(projected[atom])), "tangent_residual": diagnostic["tangent_residual_max"]})
        center = mapping[("A", 43)]["C"]
        tree = cKDTree(np.mod(x, box), boxsize=box)
        neighbors = tree.query_ball_point(np.mod(x[center], box), .4)
        local = tree.query_ball_point(np.mod(x[center], box), .8)
        for i in neighbors:
            if i == center:
                continue
            dist = float(np.linalg.norm(image(x[i]-x[center], box))*10)
            indirect = i in set().union(*(adjacency[z] for z in adjacency[center])) if adjacency[center] else False
            atom = atoms[i]
            environment_rows.append({**common, "index": int(i), "name": atom.name, "residue": atom.residue.name + atom.residue.id, "chain": atom.residue.chain.id, "element": atom.element.symbol if atom.element else "?", "distance_A": dist, "direct_bond": i in adjacency[center], "two_bonds": indirect})
        region = np.asarray(local, dtype=int)
        region_tree = cKDTree(np.mod(x[region], box), boxsize=box)
        for li, lj in region_tree.query_pairs(.2, output_type="ndarray"):
            i, j = int(region[li]), int(region[lj])
            if j in adjacency[i] or any(j in adjacency[k] for k in adjacency[i]):
                continue
            delta = image(x[j] - x[i], box)
            local_rows.append({**common, "atom_i": i, "atom_j": j, "distance_A": float(np.linalg.norm(delta)*10), "residue_i": atoms[i].residue.name + atoms[i].residue.id, "chain_i": atoms[i].residue.chain.id, "residue_j": atoms[j].residue.name + atoms[j].residue.id, "chain_j": atoms[j].residue.chain.id})
    a.outdir.mkdir(parents=True)
    for name, rows in (("geometry.csv", geometry_rows), ("environment.csv", environment_rows), ("local-pbc-pairs.csv", local_rows), ("projected-force.csv", force_rows)):
        with (a.outdir / name).open("x", newline="") as f:
            if rows:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
    local_mapping = [{"system_index": atom.index, "chain": atom.residue.chain.id, "resid": atom.residue.id, "resname": atom.residue.name, "atom": atom.name, "element": atom.element.symbol if atom.element else None} for atom in atoms if atom.residue.chain.id in ("A", "B") and atom.residue.id in ("42", "43", "44")]
    with (a.outdir / "atom-map.csv").open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(local_mapping[0])); writer.writeheader(); writer.writerows(local_mapping)
    report = {"inputs_sha256": {str(p): sha(p) for p in inputs}, "frames": 21, "post_violation_frames": list(range(7, 21)), "starting_force_A": force_rows[0], "force_at_6ps_A": force_rows[12], "local_pairs_under_2A": len(local_rows), "limitations": ["historical 7-20 ps frames follow violated stop rule; not valid baseline", "stored Cartesian forces reprojected, not independently recomputed", "0-6 ps is a single short nonequilibrium trajectory", "atom identities depend on topology/system order", "local distance flags are not chemical clash or H-bond classification", "no crystal geometry or force component decomposition", "not a comparative sweep"]}
    (a.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "inputs_sha256"}, indent=2))


if __name__ == "__main__":
    main()
