"""Replicate-level local flexibility and coupling of annotated catalytic residues."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import write_json, digest_file
from enzymekinetics_prediction.physics.dynamics import MDSettings, run_md
from enzymekinetics_prediction.physics.trajectory import load_dcd, protein_structure
from enzymekinetics_prediction.physics.geometry import STRUCTURAL_METHOD
from enzymekinetics_prediction.physics.sites import KSI_ROLES, KSI_BINDING_CANDIDATES
from .journal import Journal
from .statistics import (
    autocorrelation_diagnostics,
    paired_sign_flip,
    benjamini_hochberg,
)


def replica_sampling_statistics(rmsd_nm, state):
    """Sampling diagnostics, separate from structural coordinate processing."""
    rmsd = np.asarray(rmsd_nm, dtype=float)
    if rmsd.ndim != 1 or len(rmsd) < 2 or not np.isfinite(rmsd).all():
        raise ValueError("Need at least two finite RMSD observations")
    temperatures = state["Temperature (K)"]
    if not len(temperatures) or not np.isfinite(temperatures).all():
        raise ValueError("Need finite recorded temperatures")
    mid = len(rmsd) // 2
    return {
        "rmsd_autocorrelation": autocorrelation_diagnostics(rmsd),
        "temperature_mean_K": float(temperatures.mean()),
        "temperature_second_half_K": float(
            temperatures.iloc[len(temperatures) // 2 :].mean()
        ),
        "rmsd_first_half_mean_nm": float(rmsd[:mid].mean()),
        "rmsd_second_half_mean_nm": float(rmsd[mid:].mean()),
    }


def analyze_replica(directory, output=None):
    import mdtraj as md
    from scipy.stats import circmean, circstd

    directory = Path(directory)
    traj = load_dcd(directory / "trajectory.dcd", directory / "topology.pdb")
    original = directory
    targets = (
        "analysis.json",
        "residue-rmsf.csv",
        "sidechain-torsions.csv",
        "site-dccm.csv",
        "structure.npz",
    )
    if output is None:
        output = (
            directory / "analysis-c1"
            if any((directory / p).exists() for p in targets)
            else directory
        )
    output = Path(output)
    if output == directory:
        if any((output / name).exists() for name in targets):
            raise FileExistsError("Historical analysis must not be overwritten")
    elif output.exists() and any(output.iterdir()):
        raise FileExistsError(output)
    traj, metrics = protein_structure(traj)
    aligned = traj.atom_slice(traj.topology.select("name CA"))
    xyz = metrics["aligned_nm"]
    rmsf, dccm = metrics["rmsf_nm"], metrics["dccm"]
    rows = []
    for i, atom in enumerate(aligned.topology.atoms):
        r = atom.residue
        rows.append(
            {
                "chain": r.chain.index,
                "resid": r.resSeq,
                "resname": r.name,
                "rmsf_nm": rmsf[i],
                "role": KSI_ROLES.get(r.resSeq, ("", "scaffold"))[1],
                "catalytic": r.resSeq in [14, 38, 99],
                "positioning": r.resSeq in [39, 54, 116],
                "analog_binding_candidate": r.resSeq in KSI_BINDING_CANDIDATES,
            }
        )
    frame = pd.DataFrame(rows)
    torsions = []
    for label, function in [("chi1", md.compute_chi1), ("chi2", md.compute_chi2)]:
        indices, angles = function(traj)
        for column, atoms in enumerate(indices):
            residue = traj.topology.atom(int(atoms[1])).residue
            if residue.resSeq not in (set(KSI_ROLES) | KSI_BINDING_CANDIDATES):
                continue
            values = angles[:, column].astype(float)
            centers = np.array([-np.pi / 3, np.pi / 3, np.pi])
            delta = np.angle(np.exp(1j * (values[:, None] - centers[None, :])))
            rotamers = np.argmin(np.abs(delta), axis=1)
            torsions.append(
                {
                    "chain": residue.chain.index,
                    "resid": residue.resSeq,
                    "torsion": label,
                    "circular_mean_deg": float(
                        np.degrees(circmean(values, high=np.pi, low=-np.pi))
                    ),
                    "circular_std_deg": float(
                        np.degrees(circstd(values, high=np.pi, low=-np.pi))
                    ),
                    "rotamer_transitions_observed": int(
                        np.sum(rotamers[1:] != rotamers[:-1])
                    ),
                    **{
                        f"rotamer_{i}_fraction": float(np.mean(rotamers == i))
                        for i in range(3)
                    },
                }
            )
    pairs = []
    for i, a in enumerate(rows):
        if a["resid"] not in (set(KSI_ROLES) | KSI_BINDING_CANDIDATES):
            continue
        for j, b in enumerate(rows):
            if j <= i or b["resid"] not in (set(KSI_ROLES) | KSI_BINDING_CANDIDATES):
                continue
            pairs.append(
                {
                    "a_chain": a["chain"],
                    "a_resid": a["resid"],
                    "b_chain": b["chain"],
                    "b_resid": b["resid"],
                    "dccm": float(dccm[i, j]),
                }
            )
    state = pd.read_csv(original / "state.csv")
    result = {
        "frames": len(xyz),
        "structural_method": STRUCTURAL_METHOD,
        "analysis_output": str(output),
        "mdtraj_version": md.__version__,
        "numpy_version": np.__version__,
        "topology_sha256": digest_file(original / "topology.pdb"),
        "state_sha256": digest_file(original / "state.csv"),
        "analysis_code_sha256": {
            "ensemble.py": digest_file(__file__),
            "statistics.py": digest_file(Path(__file__).with_name("statistics.py")),
            **{
                name: digest_file(Path(__file__).parents[1] / "physics" / name)
                for name in ("geometry.py", "trajectory.py", "sites.py")
            },
        },
        **replica_sampling_statistics(metrics["rmsd_nm"], state),
        "catalytic_minus_scaffold_rmsf_nm": float(
            frame[frame.catalytic].rmsf_nm.mean()
            - frame[~frame.catalytic & ~frame.positioning].rmsf_nm.mean()
        ),
        "positioning_minus_scaffold_rmsf_nm": float(
            frame[frame.positioning].rmsf_nm.mean()
            - frame[~frame.catalytic & ~frame.positioning].rmsf_nm.mean()
        ),
        "trajectory_sha256": digest_file(original / "trajectory.dcd"),
        "limitations": [
            "Fixed-protonation apo WT dimer, not ES or constant-pH reaction dynamics",
            "Short exploratory sampling; equal initial structure, independent velocity/thermostat seeds",
            "DCCM is correlation, not causal influence or a force constant",
            "Bonded reconstruction and compact nearest-image assembly; orthorhombic cells only",
            "Replica-level significance requires a separate consistent reanalysis of the complete series",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / "residue-rmsf.csv", index=False)
    pd.DataFrame(torsions).to_csv(output / "sidechain-torsions.csv", index=False)
    pd.DataFrame(pairs).to_csv(output / "site-dccm.csv", index=False)
    np.savez_compressed(
        output / "structure.npz",
        **{k: metrics[k] for k in ("rmsd_nm", "rmsf_nm", "dccm", "radius_gyration_nm")},
    )
    write_json(output / "analysis.json", result)
    return result


def worker(arguments):
    source, output, parameters = arguments
    root = Path(source)
    run_md(root / "topology.pdb", output, MDSettings(**parameters), root / "system.xml")
    return analyze_replica(output)


def run_replicas(
    prepared, output, threads=8, equilibration_steps=5000, production_steps=10000
):
    output = Path(output)
    journal = Journal(output)
    seeds = [101, 202, 303]
    settings = [
        {
            "seed": seed,
            "threads": threads,
            "equilibration_steps": equilibration_steps,
            "production_steps": production_steps,
            "report_interval": 100,
            "protonation_variants": {"A:99": "ASH", "B:99": "ASH"},
            "minimization_iterations": 0,
        }
        for seed in seeds
    ]
    with journal.run(
        "H3-replicated-md",
        {
            "settings": settings,
            "concurrent_processes": 3,
            "starting_system": str(prepared),
        },
    ):
        with ProcessPoolExecutor(3) as pool:
            results = list(
                pool.map(
                    worker,
                    [
                        (str(prepared), str(output / f"replica-{s}"), p)
                        for s, p in zip(seeds, settings, strict=True)
                    ],
                )
            )
        tests = {
            key: paired_sign_flip([r[key] for r in results])
            for key in [
                "catalytic_minus_scaffold_rmsf_nm",
                "positioning_minus_scaffold_rmsf_nm",
            ]
        }
        corrected = benjamini_hochberg([r["p_two_sided_exact"] for r in tests.values()])
        for result, q in zip(tests.values(), corrected, strict=True):
            result["q_bh"] = float(q)
        journal.result(
            "H3-replicated-md",
            {
                "replicates": results,
                "tests": tests,
                "conclusion": "Exploratory pilot: n=3 cannot reject a two-sided exact sign-flip test at alpha=.05; no convergence claim",
            },
        )
