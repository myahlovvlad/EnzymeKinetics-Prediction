"""Benchmark prepared OpenMM systems; choose measured throughput over thread count."""

import argparse
import os
from pathlib import Path
import platform
import statistics
import time

import numpy as np

from enzymekinetics_prediction.physics.artifacts import digest_file, write_json


def hardware():
    affinity = (
        sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else list(range(os.cpu_count() or 1))
    )
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "logical_cpus": os.cpu_count(),
        "affinity_cpus": affinity,
        "memory_total_bytes": os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"),
    }


def benchmark(
    system_xml, topology, thread_counts=(1, 2, 4, 8, 12, 24), steps=50, repeats=3
):
    import openmm as mm
    from openmm import app, unit

    if (
        steps < 10
        or repeats < 2
        or any(t < 1 or t > len(os.sched_getaffinity(0)) for t in thread_counts)
    ):
        raise ValueError("Invalid benchmark limits")
    system = mm.XmlSerializer.deserialize(Path(system_xml).read_text())
    pdb = app.PDBFile(str(topology))
    if system.getNumParticles() != pdb.topology.getNumAtoms():
        raise ValueError("Benchmark System/topology atom count mismatch")
    records = []
    # Fixed order is recorded; warm-up and medians reduce initialization noise.
    for threads in thread_counts:
        integrator = mm.LangevinMiddleIntegrator(298.15, 1, 0.002)
        integrator.setRandomNumberSeed(42)
        context = mm.Context(
            system,
            integrator,
            mm.Platform.getPlatformByName("CPU"),
            {"Threads": str(threads)},
        )
        context.setPositions(pdb.positions)
        energy = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        context.setVelocitiesToTemperature(298.15, 42)
        integrator.step(10)
        durations = []
        for _ in range(repeats):
            start = time.perf_counter()
            integrator.step(steps)
            context.getState(getEnergy=True)  # synchronize native work
            durations.append(time.perf_counter() - start)
        seconds = statistics.median(durations)
        records.append(
            {
                "threads": threads,
                "seconds_per_repeat": durations,
                "steps_per_s": steps / seconds,
                "ns_per_day": steps * 0.002 / 1000 * 86400 / seconds,
                "initial_energy_kj_mol": energy,
            }
        )
        del context, integrator
    reference = records[0]["initial_energy_kj_mol"]
    for record in records:
        record["energy_difference_from_1thread_kj_mol"] = (
            record["initial_energy_kj_mol"] - reference
        )
        if not np.isclose(
            record["initial_energy_kj_mol"], reference, rtol=1e-6, atol=0.01
        ):
            raise RuntimeError(
                "Thread configuration changes reference energy beyond tolerance"
            )
    best = max(records, key=lambda row: row["steps_per_s"])
    return {
        "hardware": hardware(),
        "inputs": {
            "system_sha256": digest_file(system_xml),
            "topology_sha256": digest_file(topology),
        },
        "steps_per_repeat": steps,
        "repeats": repeats,
        "step_fs": 2,
        "records": records,
        "best_threads": best["threads"],
        "speedup_over_1thread": best["steps_per_s"] / records[0]["steps_per_s"]
        if records[0]["threads"] == 1
        else None,
        "limitations": [
            "Short single-system throughput benchmark; concurrent replicas require a separate throughput check",
            "PME/FFT, memory bandwidth and SMT can limit scaling; maximum threads need not be fastest",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4, 8, 12, 24])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--opencl", action="store_true")
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    write_json(
        args.output,
        check_opencl(args.system, args.topology)
        if args.opencl
        else benchmark(
            args.system, args.topology, args.threads, args.steps, args.repeats
        ),
    )
    print(args.output)


def check_opencl(system_xml, topology):
    """Validate a real system, not just platform enumeration; no silent GPU fallback."""
    import openmm as mm
    from openmm import app, unit

    result = {
        "backend": "OpenCL",
        "scope": "Experimental single-geometry energy/force and short stability check",
    }
    try:
        system = mm.XmlSerializer.deserialize(Path(system_xml).read_text())
        pdb = app.PDBFile(str(topology))
        integrator = mm.VerletIntegrator(0.002)
        cpu = mm.Context(
            system, integrator, mm.Platform.getPlatformByName("CPU"), {"Threads": "1"}
        )
        cpu.setPositions(pdb.positions)
        state = cpu.getState(getForces=True, getEnergy=True)
        reference_energy = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        reference_forces = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        del cpu, integrator
        integrator = mm.LangevinMiddleIntegrator(298.15, 1, 0.002)
        gpu = mm.Context(
            system,
            integrator,
            mm.Platform.getPlatformByName("OpenCL"),
            {"Precision": "mixed", "UseCpuPme": "true"},
        )
        gpu.setPositions(pdb.positions)
        state = gpu.getState(getForces=True, getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        relative = float(
            np.linalg.norm(forces - reference_forces) / np.linalg.norm(reference_forces)
        )
        result.update(
            device=gpu.getPlatform().getPropertyValue(gpu, "DeviceName"),
            relative_force_error=relative,
            energy_difference_kj_mol=float(energy - reference_energy),
        )
        if (
            not np.isfinite(relative)
            or relative > 1e-4
            or abs(energy - reference_energy) > 1
        ):
            raise ValueError("CPU/GPU energy-force comparison failed")
        gpu.setVelocitiesToTemperature(298.15, 42)
        integrator.step(10)
        start = time.perf_counter()
        integrator.step(50)
        gpu.getState(getEnergy=True)
        result.update(
            status="experimental_check_passed",
            steps_per_s=50 / (time.perf_counter() - start),
        )
    except Exception as exc:
        result.update(status="failed", error_type=type(exc).__name__, error=str(exc))
    return result


if __name__ == "__main__":
    main()
