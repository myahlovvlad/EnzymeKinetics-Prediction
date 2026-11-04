"""Replicate-level local flexibility and coupling of annotated catalytic residues."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import write_json, digest_file
from enzymekinetics_prediction.physics.dynamics import MDSettings, run_md
from enzymekinetics_prediction.physics.trajectory import load_dcd
from enzymekinetics_prediction.physics.sites import KSI_ROLES, KSI_BINDING_CANDIDATES
from .journal import Journal
from .statistics import (
    autocorrelation_diagnostics,
    paired_sign_flip,
    benjamini_hochberg,
)


def analyze_replica(directory):
    import mdtraj as md
    from scipy.stats import circmean, circstd

    directory = Path(directory)
    traj = load_dcd(directory / "trajectory.dcd", directory / "topology.pdb")
    ca = traj.topology.select("protein and name CA")
    aligned = traj.atom_slice(ca)
    aligned.superpose(aligned, 0, parallel=False)
    xyz = aligned.xyz.astype(float)
    fluct = xyz - xyz.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum(fluct**2, axis=-1), axis=0))
    covariance = np.einsum("tix,tjx->ij", fluct, fluct) / len(xyz)
    denom = np.sqrt(np.diag(covariance)[:, None] * np.diag(covariance)[None, :])
    dccm = np.divide(
        covariance, denom, out=np.zeros_like(covariance), where=denom > 1e-15
    )
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
    frame.to_csv(directory / "residue-rmsf.csv", index=False)
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
    pd.DataFrame(torsions).to_csv(directory / "sidechain-torsions.csv", index=False)
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
    pd.DataFrame(pairs).to_csv(directory / "site-dccm.csv", index=False)
    rmsd = np.sqrt(np.mean(np.sum((xyz - xyz[0]) ** 2, axis=-1), axis=1))
    state = pd.read_csv(directory / "state.csv")
    result = {
        "frames": len(xyz),
        "rmsd_autocorrelation": autocorrelation_diagnostics(rmsd),
        "temperature_mean_K": float(state["Temperature (K)"].mean()),
        "temperature_second_half_K": float(
            state["Temperature (K)"].iloc[len(state) // 2 :].mean()
        ),
        "rmsd_first_half_mean_nm": float(rmsd[: len(rmsd) // 2].mean()),
        "rmsd_second_half_mean_nm": float(rmsd[len(rmsd) // 2 :].mean()),
        "catalytic_minus_scaffold_rmsf_nm": float(
            frame[frame.catalytic].rmsf_nm.mean()
            - frame[~frame.catalytic & ~frame.positioning].rmsf_nm.mean()
        ),
        "positioning_minus_scaffold_rmsf_nm": float(
            frame[frame.positioning].rmsf_nm.mean()
            - frame[~frame.catalytic & ~frame.positioning].rmsf_nm.mean()
        ),
        "trajectory_sha256": digest_file(directory / "trajectory.dcd"),
        "limitations": [
            "Fixed-protonation apo WT dimer, not ES or constant-pH reaction dynamics",
            "Short exploratory sampling; equal initial structure, independent velocity/thermostat seeds",
            "DCCM is correlation, not causal influence or a force constant",
        ],
    }
    write_json(directory / "analysis.json", result)
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
