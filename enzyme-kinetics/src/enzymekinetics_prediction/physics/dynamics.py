"""Actual OpenMM trajectories and ensemble observables, independent of animations."""

from dataclasses import asdict, dataclass
from pathlib import Path
import random
import tempfile

import numpy as np

from .artifacts import artifact, digest_file, package_version, write_json
from .kinetics import positive
from .geometry import STRUCTURAL_METHOD, structural_metrics


@dataclass(frozen=True)
class MDSettings:
    temperature_K: float = 298.15
    ph: float = 7.0
    ionic_strength_M: float = 0.15
    step_fs: float = 2.0
    equilibration_steps: int = 100
    production_steps: int = 1000
    report_interval: int = 100
    seed: int = 42
    threads: int = 2
    padding_nm: float = 1.0
    protonation_variants: dict | None = None
    minimization_iterations: int = 200

    def __post_init__(self):
        positive(self.temperature_K, "temperature_K")
        positive(self.ionic_strength_M, "ionic_strength_M", allow_zero=True)
        if not 0 <= self.ph <= 14 or not 0 < self.step_fs <= 2:
            raise ValueError(
                "Invalid pH or timestep; default protocol supports at most 2 fs"
            )
        if (
            self.equilibration_steps < 0
            or self.report_interval < 1
            or self.production_steps < 2 * self.report_interval
            or self.threads < 1
        ):
            raise ValueError(
                "Need at least two trajectory frames and positive resource limits"
            )
        if self.padding_nm < 1:
            raise ValueError("Explicit-solvent padding must be at least 1 nm")
        if (
            not isinstance(self.minimization_iterations, int)
            or self.minimization_iterations < 0
        ):
            raise ValueError("minimization_iterations must be a nonnegative integer")


def ensemble_observables(coordinates_nm, distance_pairs=(), productive_upper_nm=None):
    """Summaries for whole, coherently imaged coordinates (no PBC guessing)."""
    xyz = np.asarray(coordinates_nm, dtype=float)
    if (
        xyz.ndim != 3
        or xyz.shape[0] < 2
        or xyz.shape[2] != 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected finite [frames>=2, atoms, 3] coordinates in nm")
    metrics = structural_metrics(xyz)
    rmsd = metrics["rmsd_nm"]
    rg = metrics["radius_gyration_nm"]
    dccm = metrics["dccm"]
    off_diagonal = dccm[~np.eye(xyz.shape[1], dtype=bool)]
    values = {
        "rmsd_mean_nm": float(rmsd.mean()),
        "rmsd_std_nm": float(rmsd.std()),
        "radius_gyration_mean_nm": float(rg.mean()),
        "dccm_abs_offdiagonal_mean": float(np.abs(off_diagonal).mean())
        if len(off_diagonal)
        else 0.0,
    }
    units = {
        "rmsd_mean_nm": "nm",
        "rmsd_std_nm": "nm",
        "radius_gyration_mean_nm": "nm",
        "dccm_abs_offdiagonal_mean": "dimensionless",
    }
    if productive_upper_nm is not None and len(productive_upper_nm) != len(
        distance_pairs
    ):
        raise ValueError(
            "One explicit productive threshold is required per distance pair"
        )
    productive = np.ones(len(xyz), dtype=bool)
    for i, pair in enumerate(distance_pairs):
        a, b = pair
        if min(a, b) < 0 or max(a, b) >= xyz.shape[1] or a == b:
            raise ValueError("Invalid trajectory atom indices")
        distances = np.linalg.norm(xyz[:, a] - xyz[:, b], axis=1)
        values[f"distance_{i}_mean_nm"] = float(distances.mean())
        units[f"distance_{i}_mean_nm"] = "nm"
        if productive_upper_nm is not None:
            positive(productive_upper_nm[i], "productive distance threshold")
            productive &= distances <= productive_upper_nm[i]
    if productive_upper_nm is not None and distance_pairs:
        values["productive_fraction"] = float(productive.mean())
        units["productive_fraction"] = "fraction"
    return values, units


def run_md(
    pdb_path,
    output,
    settings=MDSettings(),
    system_xml=None,
    distance_pairs=(),
    productive_upper_nm=None,
    restart_checkpoint=None,
):
    try:
        import openmm as mm
        from openmm import app, unit
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError("Install the physics extra for OpenMM and MDTraj") from exc
    pdb_path = Path(pdb_path).resolve()
    output = Path(output).resolve()
    if not pdb_path.is_file():
        raise FileNotFoundError(pdb_path)
    if output.exists() and any(output.iterdir()):
        raise ValueError("MD output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    pdb = app.PDBFile(str(pdb_path))
    inputs = {"pdb": str(pdb_path), "pdb_sha256": digest_file(pdb_path)}
    if system_xml:
        system = mm.XmlSerializer.deserialize(Path(system_xml).read_text())
        topology, positions = pdb.topology, pdb.positions
        if system.getNumParticles() != topology.getNumAtoms():
            raise ValueError("Prepared System and PDB atom counts differ")
        inputs["system_xml_sha256"] = digest_file(system_xml)
        parameterization = "user-supplied System XML; original atom order required"
    else:
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        modeller = app.Modeller(pdb.topology, pdb.positions)
        known_residues = {f"{r.chain.id}:{r.id}" for r in modeller.topology.residues()}
        if set(settings.protonation_variants or {}) - known_residues:
            raise ValueError("Protonation override refers to an unknown chain/residue")
        preparation_platform = mm.Platform.getPlatformByName("CPU")
        previous_threads = preparation_platform.getPropertyDefaultValue("Threads")
        rng_state = random.getstate()
        try:
            preparation_platform.setPropertyDefaultValue(
                "Threads", str(settings.threads)
            )
            random.seed(settings.seed)
            modeller.addHydrogens(
                forcefield,
                pH=settings.ph,
                platform=preparation_platform,
                variants=[
                    (settings.protonation_variants or {}).get(f"{r.chain.id}:{r.id}")
                    for r in modeller.topology.residues()
                ],
            )
        finally:
            preparation_platform.setPropertyDefaultValue("Threads", previous_threads)
            random.setstate(rng_state)
        modeller.addSolvent(
            forcefield,
            padding=settings.padding_nm * unit.nanometer,
            ionicStrength=settings.ionic_strength_M * unit.molar,
            neutralize=True,
        )
        topology, positions = modeller.topology, modeller.positions
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1 * unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
        )
        parameterization = (
            "amber14-all.xml + amber14/tip3pfb.xml; standard residues only"
        )
    integrator = mm.LangevinMiddleIntegrator(
        settings.temperature_K * unit.kelvin,
        1 / unit.picosecond,
        settings.step_fs * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(settings.seed)
    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(
        topology, system, integrator, platform, {"Threads": str(settings.threads)}
    )
    if restart_checkpoint:
        if not system_xml:
            raise ValueError(
                "Checkpoint continuation requires the original prepared System XML"
            )
        inputs["restart_checkpoint_sha256"] = digest_file(restart_checkpoint)
        simulation.context.loadCheckpoint(Path(restart_checkpoint).read_bytes())
    else:
        simulation.context.setPositions(positions)
        if settings.minimization_iterations:
            simulation.minimizeEnergy(maxIterations=settings.minimization_iterations)
        simulation.context.setVelocitiesToTemperature(
            settings.temperature_K * unit.kelvin, settings.seed
        )
    simulation.reporters.append(
        app.CheckpointReporter(str(output / "checkpoint.chk"), 1000)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(output / "progress.csv"),
            settings.report_interval,
            step=True,
            temperature=True,
            potentialEnergy=True,
            elapsedTime=True,
            speed=True,
        )
    )
    simulation.step(settings.equilibration_steps)
    initial = simulation.context.getState(getPositions=True)
    topology_path = output / "topology.pdb"
    with topology_path.open("w") as handle:
        app.PDBFile.writeFile(topology, initial.getPositions(), handle, keepIds=True)
    (output / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    trajectory_path = output / "trajectory.dcd"
    simulation.reporters.append(
        app.DCDReporter(str(trajectory_path), settings.report_interval)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(output / "state.csv"),
            settings.report_interval,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
        )
    )
    simulation.step(settings.production_steps)
    checkpoint = simulation.context.createCheckpoint()
    (output / "checkpoint.chk").write_bytes(checkpoint)
    # MDTraj's native DCD reader requires ASCII filenames. The workspace may
    # legitimately contain Cyrillic; an ASCII temporary symlink avoids copying
    # a potentially large trajectory and preserves the recorded input path.
    with tempfile.TemporaryDirectory(prefix="enzyme-md-") as temporary:
        link = Path(temporary) / "trajectory.dcd"
        link.symlink_to(trajectory_path)
        trajectory = md.load_dcd(str(link), top=md.Topology.from_openmm(topology))
    ca_indices = trajectory.topology.select("protein and name CA")
    if not len(ca_indices):
        raise ValueError("No protein CA atoms for ensemble observables")
    from .trajectory import protein_structure

    prepared, structure = protein_structure(
        trajectory, periodic=system.usesPeriodicBoundaryConditions()
    )
    # Preserve the existing CA-relative distance-pair indexing.
    values, units = ensemble_observables(
        structure["coordinates_nm"],
        distance_pairs,
        productive_upper_nm,
    )
    values.update(
        production_time_ps=settings.production_steps * settings.step_fs / 1000,
        frames=int(trajectory.n_frames),
    )
    units.update(production_time_ps="ps", frames="count")
    result = artifact(
        "OpenMM",
        inputs,
        {
            **asdict(settings),
            "parameterization": parameterization,
            "openmm_version": package_version("openmm"),
            "mdtraj_version": package_version("mdtraj"),
            "trajectory_sha256": digest_file(trajectory_path),
            "observable_selection": "protein CA",
            "structural_method": STRUCTURAL_METHOD,
            "periodic_analysis": system.usesPeriodicBoundaryConditions(),
            "distance_pairs_CA_indices": list(distance_pairs),
            "productive_upper_nm": productive_upper_nm,
        },
        values,
        units,
        limitations=[
            "Classical fixed-topology MD does not simulate the chemical reaction",
            "Default preparation is an apo standard-residue protein; ES ligands require a supplied parameterized System",
            "Hydrogen assignment at nominal pH is not constant-pH sampling or a validated catalytic protonation state",
            "Reported moments and DCCM do not establish convergence; short runs are backend checks",
        ],
    )
    write_json(output / "result.json", result)
    del simulation, integrator
    return result
