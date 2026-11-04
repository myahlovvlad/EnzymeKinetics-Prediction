"""All-atom trajectory distances and projected electric fields along an explicit pair.

The field model is deliberately the nonperiodic, homogeneous screened Coulomb
approximation. It is not the PME force field used to generate the trajectory,
APBS, a bond-breaking trajectory, or a quantum electric field.
"""

import argparse
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .artifacts import artifact, digest_file, write_json
from .electrostatics import screened_coulomb
from .geometry import make_whole, orthorhombic_lengths, structural_metrics


def load_dcd(path, topology):
    import mdtraj as md

    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(prefix="enzyme-trace-") as temporary:
        # Native MDTraj readers reject Unicode filenames, unlike Python I/O.
        dcd, pdb = Path(temporary) / "trajectory.dcd", Path(temporary) / "topology.pdb"
        dcd.symlink_to(path)
        pdb.symlink_to(Path(topology).resolve())
        return md.load_dcd(str(dcd), top=str(pdb))


def protein_structure(traj, *, periodic=None):
    """Return whole protein coordinates and canonical CA structural metrics.

    Reconstruct full bonded molecules before selecting CA. Every protein chain
    must lie in one connected bond graph; incomplete topologies fail explicitly.
    periodic=False is reserved for a known nonperiodic System, since DCD may
    contain default box vectors even for such systems.
    """
    protein = traj.atom_slice(traj.topology.select("protein"))
    if not protein.n_atoms:
        raise ValueError("No protein atoms for structural analysis")
    boxes = protein.unitcell_vectors
    if periodic is None:
        periodic = boxes is not None
    if periodic:
        if boxes is None:
            raise ValueError(
                "Periodic structural analysis requires per-frame box vectors"
            )
        orthorhombic_lengths(boxes)
        molecules = list(protein.topology.find_molecules())
        membership = {a.index: i for i, group in enumerate(molecules) for a in group}
        for chain in protein.topology.chains:
            components = {membership[a.index] for a in chain.atoms}
            if len(components) > 1:
                raise ValueError("Incomplete protein bond topology: disconnected chain")
        bonds = [(a.index, b.index) for a, b in protein.topology.bonds]
        whole, _, _ = make_whole(protein.xyz, boxes, bonds)
        # MDTraj storage is float32; metrics retain the double reconstruction.
    else:
        whole = protein.xyz.astype(float)
    ca = protein.topology.select("name CA")
    if not len(ca):
        raise ValueError("No protein CA atoms for structural analysis")
    metrics = structural_metrics(whole[:, ca])
    protein.xyz = whole
    return protein, metrics


def projected_fields(
    coordinates_nm,
    charges_e,
    pair,
    source_indices,
    *,
    temperature_K=298.15,
    ionic_strength_M=0.15,
    dielectric=78.5,
):
    xyz = np.asarray(coordinates_nm, float)
    charges = np.asarray(charges_e, float)
    if (
        xyz.ndim != 3
        or xyz.shape[2] != 3
        or len(xyz) == 0
        or charges.shape != (xyz.shape[1],)
        or not np.isfinite(xyz).all()
        or not np.isfinite(charges).all()
    ):
        raise ValueError(
            "Expected finite frames x atoms x 3 coordinates and matching charges"
        )
    a, b = pair
    if a == b or min(a, b) < 0 or max(a, b) >= xyz.shape[1]:
        raise ValueError("Invalid explicit reaction atom pair")
    sources = np.asarray(sorted(set(source_indices)), dtype=int)
    if (
        sources.size == 0
        or sources.min() < 0
        or sources.max() >= xyz.shape[1]
        or a in sources
        or b in sources
    ):
        raise ValueError(
            "Sources must be nonempty, valid, and exclude both probe atoms"
        )
    rows = []
    for i, frame in enumerate(xyz):
        axis = frame[b] - frame[a]
        distance = np.linalg.norm(axis)
        if distance <= 1e-10:
            raise ValueError("Probe atom positions coincide")
        midpoint = (frame[a] + frame[b]) * 5  # nm -> A, then midpoint
        phi, field = screened_coulomb(
            frame[sources] * 10,
            charges[sources],
            [midpoint],
            temperature_K=temperature_K,
            ionic_strength_M=ionic_strength_M,
            dielectric=dielectric,
        )
        rows.append(
            {
                "frame": i,
                "distance_nm": distance,
                "potential_V": phi[0],
                "field_projection_V_per_A": field[0] @ (axis / distance),
                "field_norm_V_per_A": np.linalg.norm(field[0]),
            }
        )
    return pd.DataFrame(rows)


def analyze_trajectory(
    trajectory,
    topology,
    system_xml,
    pair,
    source_selection,
    label,
    output,
    *,
    prepared_nonperiodic=False,
    **conditions,
):
    import openmm as mm
    from openmm import unit

    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose a new trace output directory")
    traj = load_dcd(trajectory, topology)
    system = mm.XmlSerializer.deserialize(Path(system_xml).read_text())
    if system.usesPeriodicBoundaryConditions() and not prepared_nonperiodic:
        raise ValueError(
            "Periodic input requires explicitly prepared nonperiodic coordinates; "
            "use prepared_nonperiodic=True only after image preparation. "
            "The field model is not a periodic sum."
        )
    if system.getNumParticles() != traj.n_atoms:
        raise ValueError(
            "System and trajectory atom counts differ; original atom order is required"
        )
    forces = [
        force for force in system.getForces() if isinstance(force, mm.NonbondedForce)
    ]
    if len(forces) != 1:
        raise ValueError(
            "Exactly one standard NonbondedForce is required for charge extraction"
        )
    charges = np.array(
        [
            forces[0].getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
            for i in range(traj.n_atoms)
        ]
    )
    sources = sorted(set(traj.topology.select(source_selection)) - set(pair))
    frame = projected_fields(traj.xyz, charges, pair, sources, **conditions)
    frame.insert(1, "time_ps", traj.time)
    result = artifact(
        "trajectory-screened-coulomb",
        {
            **{
                name + "_sha256": digest_file(p)
                for name, p in [
                    ("trajectory", trajectory),
                    ("topology", topology),
                    ("system", system_xml),
                ]
            },
            "label": label,
        },
        {
            "atom_pair_zero_based": list(pair),
            "source_selection": source_selection,
            "excluded_probe_atoms": list(pair),
            "source_atom_count": len(sources),
            "prepared_nonperiodic": prepared_nonperiodic,
            **conditions,
        },
        {
            "frames": len(frame),
            "distance_mean_nm": float(frame.distance_nm.mean()),
            "field_projection_mean_V_per_A": float(
                frame.field_projection_V_per_A.mean()
            ),
            "field_projection_std_V_per_A": float(
                frame.field_projection_V_per_A.std(ddof=0)
            ),
        },
        {
            "frames": "count",
            "distance_mean_nm": "nm",
            "field_projection_mean_V_per_A": "Volt/Angstrom",
            "field_projection_std_V_per_A": "Volt/Angstrom",
        },
        limitations=[
            "Homogeneous nonperiodic screened Coulomb from fixed partial charges; not PME, APBS or QM",
            "Source selection and probe orientation must be justified for the reaction; no automatic TS assignment",
            "No periodic minimum-image summation; boundary crossings require separately prepared unwrapped/aligned trajectories",
            "No convergence or causal rate-acceleration claim follows from these moments",
        ],
    )
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / "trace.csv", index=False)
    write_json(output / "trace.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["trajectory", "topology", "system", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--pair", nargs=2, type=int, required=True)
    parser.add_argument("--source-selection", default="protein")
    parser.add_argument("--label", required=True)
    parser.add_argument(
        "--prepared-nonperiodic",
        action="store_true",
        help="Assert that periodic images were prepared explicitly for the nonperiodic field model",
    )
    parser.add_argument("--temperature-K", type=float, default=298.15)
    parser.add_argument("--ionic-strength-M", type=float, required=True)
    parser.add_argument("--dielectric", type=float, required=True)
    args = parser.parse_args()
    result = analyze_trajectory(
        args.trajectory,
        args.topology,
        args.system,
        args.pair,
        args.source_selection,
        args.label,
        args.output,
        prepared_nonperiodic=args.prepared_nonperiodic,
        temperature_K=args.temperature_K,
        ionic_strength_M=args.ionic_strength_M,
        dielectric=args.dielectric,
    )
    print(result["values"])


if __name__ == "__main__":
    main()
