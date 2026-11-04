"""Evaluate the fixed-geometry constraint force projection, without minimization."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm
from openmm import unit

from constraint_projection import project_forces


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-npy", type=Path, required=True)
    parser.add_argument("--system-xml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if openmm.__version__ != "8.6.1":
        raise RuntimeError(f"OpenMM 8.6.1 required, found {openmm.__version__}")
    system = openmm.XmlSerializer.deserialize(args.system_xml.read_text())
    positions = np.load(args.input_npy, allow_pickle=False)
    if (
        positions.shape != (system.getNumParticles(), 3)
        or not np.isfinite(positions).all()
    ):
        raise ValueError("invalid coordinates")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"}
    )
    context.setPositions(positions * unit.nanometer)
    state = context.getState(getForces=True, getEnergy=True)
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
    )
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not np.isfinite(energy) or not np.isfinite(forces).all():
        raise ValueError("nonfinite energy or forces")
    constraints = [
        system.getConstraintParameters(i) for i in range(system.getNumConstraints())
    ]
    pairs = np.array([(int(i), int(j)) for i, j, _ in constraints], dtype=int)
    targets = np.array([d.value_in_unit(unit.nanometer) for _, _, d in constraints])
    errors = np.abs(
        np.linalg.norm(positions[pairs[:, 0]] - positions[pairs[:, 1]], axis=1)
        - targets
    )
    relative = errors / targets
    projected, diagnostic = project_forces(positions, forces, pairs)
    norms = np.linalg.norm(projected, axis=1)
    tolerance = integrator.getConstraintTolerance()
    report = {
        "input_sha256": {
            str(p): sha(p)
            for p in (
                args.input_npy,
                args.system_xml,
                Path(__file__),
                Path(__file__).with_name("constraint_projection.py"),
            )
        },
        "openmm_version": openmm.__version__,
        "platform": "CPU",
        "threads": 4,
        "classification": "fixed-coordinate force projection diagnostic, not minimization",
        "constraints": len(pairs),
        "energy_kj_mol": energy,
        "constraint_tolerance_relative": tolerance,
        "constraint_error_max_nm": float(errors.max()),
        "constraint_error_rms_nm": float(np.sqrt(np.mean(errors**2))),
        "constraint_error_relative_max": float(relative.max()),
        "constraints_above_relative_tolerance": int(
            np.count_nonzero(relative > tolerance)
        ),
        "projected_force_max_kj_mol_nm": float(norms.max()),
        "projected_force_rms_component_kj_mol_nm": float(
            np.sqrt(np.mean(projected**2))
        ),
        "projection_diagnostics": diagnostic,
        "note": "Projection changes forces only, not positions: constraint distance residuals are identical before and after by construction. Tangency residual checks J Fproj, not minimizer convergence.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
