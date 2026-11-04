"""Compare OpenMM CPU and OpenCL energies/forces at identical coordinates."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit


def compare_forces(cpu, gpu):
    cpu, gpu = np.asarray(cpu), np.asarray(gpu)
    if cpu.shape != gpu.shape or not np.isfinite(cpu).all() or not np.isfinite(gpu).all():
        raise ValueError("Non-finite or mismatched forces")
    diff = gpu - cpu
    return {"rms_component_error_kj_mol_nm": float(np.sqrt(np.mean(diff**2))),
            "max_component_error_kj_mol_nm": float(np.max(np.abs(diff))),
            "cpu_rms_component_kj_mol_nm": float(np.sqrt(np.mean(cpu**2)))}


def evaluate(name, system_xml, positions, threads, steps, precision, cpu_pme):
    system = openmm.XmlSerializer.deserialize(system_xml)
    integ = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    props = {"Threads": str(threads)} if name == "CPU" else {"Precision": precision}
    if name == "OpenCL" and cpu_pme:
        props["UseCpuPme"] = "true"
    start = time.perf_counter()
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName(name), props)
    context.setPositions(positions)
    state = context.getState(getEnergy=True, getForces=True)
    creation_and_first_s = time.perf_counter() - start
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    start = time.perf_counter()
    if steps:
        integ.step(steps)
    after = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    elapsed = time.perf_counter() - start
    if not np.isfinite(energy) or not np.isfinite(after):
        raise ValueError(f"{name} returned non-finite energy")
    result = {"energy_kj_mol": float(energy), "energy_after_steps_kj_mol": float(after),
              "context_and_first_state_seconds": creation_and_first_s,
              "steps": steps, "steps_seconds": elapsed, "properties": props}
    del context, integ
    return result, forces


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--system", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--steps", type=int, default=2)
    p.add_argument("--precision", choices=("single", "mixed", "double"), default="single")
    p.add_argument("--cpu-pme", action="store_true")
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose a new output path")
    if args.steps < 0 or args.threads < 1:
        raise ValueError("Invalid step or thread count")
    pdb = app.PDBFile(str(args.pdb))
    system_xml = args.system.read_text()
    cpu, cpu_f = evaluate("CPU", system_xml, pdb.positions, args.threads, args.steps, args.precision, args.cpu_pme)
    gpu, gpu_f = evaluate("OpenCL", system_xml, pdb.positions, args.threads, args.steps, args.precision, args.cpu_pme)
    report = {"input_sha256": {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in (args.pdb, args.system)},
              "openmm_version": openmm.__version__, "atoms": len(pdb.positions),
              "cpu": cpu, "opencl": gpu, "energy_difference_kj_mol": gpu["energy_kj_mol"] - cpu["energy_kj_mol"],
              "force_comparison": compare_forces(cpu_f, gpu_f),
              "limitations": ["OpenMM warns Rusticl implementation unsupported", "Single snapshot and a few steps are not a full validation or benchmark", "No claims of stable trajectory or scientific accuracy"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
