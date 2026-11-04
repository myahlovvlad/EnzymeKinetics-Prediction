"""Replay saved P00947 H3 trajectories and C1 statistics/plots, without simulation."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import (
    digest_file,
    package_version,
    write_json,
)
from enzymekinetics_prediction.physics.geometry import (
    STRUCTURAL_METHOD,
    structural_metrics,
)
from enzymekinetics_prediction.physics.trajectory import load_dcd
from enzymekinetics_prediction.research import ensemble
from enzymekinetics_prediction.research.statistics import (
    autocorrelation_diagnostics,
    paired_sign_flip,
    benjamini_hochberg,
)

COMPARISONS = (
    "catalytic_minus_scaffold_rmsf_nm",
    "positioning_minus_scaffold_rmsf_nm",
)


def compare_replicas(results):
    """The registered H3 comparison: independent replica contrasts, two BH tests."""
    tests = {}
    for key in COMPARISONS:
        differences = [r[key] for r in results]
        tests[key] = {**paired_sign_flip(differences), "differences_nm": differences}
    adjusted = benjamini_hochberg([r["p_two_sided_exact"] for r in tests.values()])
    for result, q in zip(tests.values(), adjusted, strict=True):
        result["q_bh"] = float(q)
    return tests


def figures(output, records):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(records), 1, figsize=(9, 2.7 * len(records)), squeeze=False
    )
    for ax, record in zip(axes[:, 0], records, strict=True):
        ax.plot(record["raw"]["rmsd_nm"], label="Raw-image control")
        ax.plot(record["corrected"]["rmsd_nm"], label="C1 reconstructed")
        ax.set(
            xlabel="Saved frame (zero based)",
            ylabel="CA RMSD, nm",
            title=record["name"],
        )
        ax.legend()
    fig.tight_layout()
    fig.savefig(output / "rmsd.png", dpi=150)
    plt.close(fig)

    # Recreate the historical replicas.png from one consistent C1 result set.
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    for record in records:
        table = record["residues"]
        chain = table[table.chain == table.chain.min()]
        axes[0].plot(chain.resid, chain.rmsf_nm, label=record["name"])
        state = record["state"]
        axes[1].plot(state.iloc[:, 0], state["Temperature (K)"], label=record["name"])
    axes[0].set(
        xlabel="Author residue, first chain",
        ylabel="CA RMSF, nm",
        title="C1: short exploratory apo-dimer trajectories",
    )
    axes[1].set(xlabel="Recorded integrator step", ylabel="Temperature, K")
    for ax in axes:
        ax.legend()
    fig.tight_layout()
    fig.savefig(output / "replicas.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for record in records:
        axes[0].plot(record["sampling"]["acf"], label=record["name"])
    names = [r["name"] for r in records]
    axes[0].set(
        xlabel="Lag, saved frames",
        ylabel="RMSD autocorrelation",
        title="Biased ACF; short-series diagnostic",
    )
    axes[0].legend()
    axes[1].bar(names, [r["sampling"]["effective_samples"] for r in records])
    axes[1].set(
        ylabel="Exploratory effective samples", title="Not independent replicates"
    )
    fig.tight_layout()
    fig.savefig(output / "sampling.png", dpi=150)
    plt.close(fig)

    for record in records:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, key, title in zip(
            axes,
            ("raw", "corrected"),
            ("Raw-image control", "C1 reconstructed"),
            strict=True,
        ):
            im = ax.imshow(
                record[key]["dccm"], vmin=-1, vmax=1, cmap="coolwarm", origin="lower"
            )
            ax.set(title=title, xlabel="Ordered CA index", ylabel="Ordered CA index")
        fig.colorbar(im, ax=axes.tolist(), label="DCCM, dimensionless", shrink=0.8)
        fig.suptitle(record["name"] + ": correlation is not causation")
        fig.savefig(output / (record["name"] + "-dccm.png"), dpi=150)
        plt.close(fig)


def reassess(replicas, output):
    sources = [Path(p).resolve() for p in replicas]
    names = [p.name for p in sources]
    if not 2 <= len(sources) <= 20 or len(set(names)) != len(names):
        raise ValueError("Need 2..20 uniquely named independent replica directories")
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    required = ("trajectory.dcd", "topology.pdb", "state.csv")
    missing = [
        str(p / name) for p in sources for name in required if not (p / name).is_file()
    ]
    if missing:
        raise FileNotFoundError("BLOCKED: missing inputs: " + ", ".join(missing))
    output.mkdir(parents=True)
    source_hashes = {
        str(p / name): digest_file(p / name) for p in sources for name in required
    }
    records, results = [], []
    for source in sources:
        destination = output / source.name
        print("Reanalyzing " + source.name, flush=True)
        result = ensemble.analyze_replica(source, destination)
        results.append(result)
        with np.load(destination / "structure.npz") as saved:
            corrected = {key: saved[key].copy() for key in saved.files}
        traj = load_dcd(source / "trajectory.dcd", source / "topology.pdb")
        raw = structural_metrics(
            traj.xyz[:, traj.topology.select("protein and name CA")]
        )
        raw = {key: raw[key] for key in corrected}
        np.savez_compressed(destination / "raw-image-control.npz", **raw)
        sampling = autocorrelation_diagnostics(corrected["rmsd_nm"], include_acf=True)
        if sampling["acf"] is None:
            raise ValueError(
                "RMSD sampling is numerically unresolved; do not invent ESS"
            )
        write_json(destination / "sampling.json", sampling)
        records.append(
            {
                "name": source.name,
                "raw": raw,
                "corrected": corrected,
                "sampling": sampling,
                "residues": pd.read_csv(destination / "residue-rmsf.csv"),
                "state": pd.read_csv(source / "state.csv"),
            }
        )
    summary = {
        "method": STRUCTURAL_METHOD,
        "annotation_protocol": "P00947 H3; existing KSI roles, independent trajectories are experimental units",
        "replicas": {name: result for name, result in zip(names, results, strict=True)},
        "tests": compare_replicas(results),
        "limitations": [
            "n=3 in the registered study: minimum two-sided exact p is .25; no significance at .05",
            "Exploratory short-series ACF/ESS; shared starting structure; independence/exchangeability are assumptions",
            "No equilibrium, catalytic mechanism, causal coupling or absence-of-effect claim",
        ],
    }
    write_json(output / "summary.json", summary)
    figures(output, records)
    if any(digest_file(p) != h for p, h in source_hashes.items()):
        raise RuntimeError("Source changed during reanalysis")
    package = Path(ensemble.__file__).parents[1]
    code = [
        package / name
        for name in (
            "physics/geometry.py",
            "physics/trajectory.py",
            "physics/sites.py",
            "research/ensemble.py",
            "research/statistics.py",
        )
    ]
    manifest = {
        "schema_version": 1,
        "status": "completed",
        "scientific_status": "exploratory_not_validated",
        "analysis_only": True,
        "source_sha256": source_hashes,
        "code_sha256": {str(p.relative_to(package)): digest_file(p) for p in code},
        "driver_sha256": digest_file(__file__),
        "versions": {
            name: package_version(name)
            for name in ("numpy", "scipy", "pandas", "mdtraj", "matplotlib")
        },
        "seed": None,
        "seed_reason": "deterministic reanalysis; original trajectories use registered independent seeds",
        "output_sha256": {
            str(p.relative_to(output)): digest_file(p)
            for p in sorted(output.rglob("*"))
            if p.is_file()
        },
    }
    write_json(output / "manifest.json", manifest)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replica", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = reassess(args.replica, args.output)
    print(json.dumps(summary["tests"], indent=2))


if __name__ == "__main__":
    main()
