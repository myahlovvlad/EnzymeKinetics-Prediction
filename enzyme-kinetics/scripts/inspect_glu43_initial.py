"""Initial GLU43 geometry/provenance; crystal is not an NVT box."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm
from Bio.PDB import PDBParser
from openmm import app, unit

from glu43_geometry import geometry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "positions-npy",
        "topology-pdb",
        "system-xml",
        "crystal-pdb",
        "outdir",
    ):
        p.add_argument("--" + name, required=True, type=Path)
    a = p.parse_args()
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    if openmm.__version__ != "8.6.1":
        raise ValueError("requires OpenMM 8.6.1")
    system = openmm.XmlSerializer.deserialize(a.system_xml.read_text())
    top = app.PDBFile(str(a.topology_pdb)).topology
    atoms = list(top.atoms())
    coords = np.load(a.positions_npy, allow_pickle=False)
    if (
        len(atoms) != system.getNumParticles()
        or coords.shape != (len(atoms), 3)
        or not np.isfinite(coords).all()
    ):
        raise ValueError("incompatible coordinates, topology or system")
    box = np.asarray(
        [v.value_in_unit(unit.nanometer) for v in system.getDefaultPeriodicBoxVectors()]
    )
    if not np.allclose(box, np.diag(np.diag(box)), rtol=0, atol=1e-12):
        raise ValueError("only orthorhombic simulation box supported")
    lengths = np.diag(box)
    crystal = PDBParser(QUIET=True).get_structure("8CHO", str(a.crystal_pdb))
    nonbonded = next(
        f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)
    )
    result = {}
    for model, chain in ((0, "A"), (1, "B")):
        residues = [r for r in top.residues() if r.chain.id == chain and r.id == "43"]
        if len(residues) != 1 or residues[0].name != "GLU":
            raise ValueError(f"Expected GLU43 in chain {chain}")
        residue = residues[0]
        mapped = {atom.name: atom.index for atom in residue.atoms()}
        if not {"C", "O", "CG", "CD", "OE1", "OE2"} <= mapped.keys():
            raise ValueError("missing residue atoms")
        original = crystal[model]["A"][43]
        orig_xyz = {
            atom.name: np.asarray(atom.coord, dtype=float) / 10
            for atom in original.get_atoms()
        }
        # Original crystal coordinates are local, not periodic in the simulation
        # box. Use direct vectors, with a safely oversized box for this residue.
        if not {"C", "O", "CG", "CD", "OE1", "OE2"} <= orig_xyz.keys():
            raise ValueError("missing crystal atoms")
        charges = {
            name: nonbonded.getParticleParameters(mapped[name])[0].value_in_unit(
                unit.elementary_charge
            )
            for name in ("CD", "OE1", "OE2")
        }
        result[chain] = {
            "residue_name_in_topology": residue.name,
            "model_atom_indices": mapped,
            "sidechain_OE_hydrogen_names": [
                name for name in mapped if name in ("HE1", "HE2")
            ],
            "sidechain_charges_e": charges,
            "model_geometry": geometry(
                {name: coords[i] for name, i in mapped.items()}, lengths
            ),
            "crystal_geometry": geometry(orig_xyz, np.array([1e6] * 3)),
            "interpretation": "GLU topology and absence of OE hydrogens support modeled deprotonated glutamate, not experimentally established protonation",
        }
    inputs = (
        a.positions_npy,
        a.topology_pdb,
        a.system_xml,
        a.crystal_pdb,
        Path(__file__),
        Path(__file__).with_name("glu43_geometry.py"),
    )
    report = {
        "input_sha256": {
            str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs
        },
        "openmm_version": openmm.__version__,
        "chains": result,
        "limitations": [
            "crystal hydrogen/protonation not independently established",
            "coordinates of model use periodic minimum image; crystal geometry uses direct local coordinates",
            "atom ordering assumed",
            "no time evolution, local clash analysis, or force decomposition",
        ],
    }
    a.outdir.mkdir(parents=True)
    (a.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
