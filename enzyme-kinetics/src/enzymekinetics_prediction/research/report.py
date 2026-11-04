"""Portable evidence and readable diary for registered active-site experiments."""

import argparse
from datetime import datetime
from pathlib import Path
import json
import shutil
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import write_json, digest_file


def export(source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for p in [
        *source.glob("H*.json"),
        source / "registration.json",
        source / "events.jsonl",
        source / "assembly.json",
        source / "conditions-factorial.csv",
        source / "torsions.csv",
        source / "torsion-coupling.csv",
        source / "conditions-factorial-focused.csv",
        source / "performance-WT-throughput.json",
    ]:
        if p.exists():
            shutil.copyfile(p, output / p.name)
    evidence = {}
    for p in source.glob("results/*.json"):
        dest = output / "history" / p.name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(p, dest)
    replicas = []
    for folder in sorted(source.glob("replica-*")):
        if not (folder / "analysis.json").exists():
            continue
        destination = output / folder.name
        destination.mkdir(exist_ok=True)
        for name in [
            "analysis.json",
            "residue-rmsf.csv",
            "site-dccm.csv",
            "sidechain-torsions.csv",
            "state.csv",
            "result.json",
        ]:
            p = folder / name
            if p.exists():
                shutil.copyfile(p, destination / name)
        replicas.append(folder.name)
        evidence[folder.name] = {
            p.name: digest_file(p) for p in folder.iterdir() if p.is_file()
        }
    write_json(output / "artifact-checksums.json", evidence)
    sites = json.loads((source / "H1-static-sites.json").read_text())
    rows = []
    for structure, value in sites.items():
        for record in value["values"]["residues"]:
            rows.append(
                {
                    "structure": structure,
                    **{
                        k: v
                        for k, v in record.items()
                        if k != "heavy_atom_coordinates_A"
                    },
                }
            )
    pd.DataFrame(rows).to_csv(output / "site-residues.csv", index=False)
    refined = (source / "conditions-factorial-focused.csv").exists()
    factorial = pd.read_csv(
        source
        / (
            "conditions-factorial-focused.csv"
            if refined
            else "conditions-factorial.csv"
        )
    )
    grid = factorial[factorial.ionic_strength_M == 0.15].pivot(
        index="ph", columns="temperature_K", values="field_norm_V_per_A"
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    im = axes[0].imshow(grid.values, origin="lower", aspect="auto")
    axes[0].set_xticks(range(3), ["15", "25", "37"])
    axes[0].set_yticks(range(3), grid.index)
    axes[0].set(
        xlabel="Temperature, C",
        ylabel="pH",
        title="APBS129 focused field norm, V/A"
        if refined
        else "APBS65 field norm, V/A (grid-sensitive)",
    )
    fig.colorbar(im, ax=axes[0])
    angles = pd.read_csv(source / "torsions.csv")
    for resid, frame in angles.groupby("resid"):
        axes[1].plot(frame.rotation_deg, frame.delta_potential_kj_mol, label=str(resid))
    axes[1].set(
        xlabel="chi1 rotation, degrees",
        ylabel="Conditional delta potential, kJ/mol",
        title="Frozen environment: not PMF or a reaction barrier",
    )
    axes[1].set_yscale("symlog", linthresh=10)
    axes[1].legend(title="Residue", fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "static-results.png", dpi=160)
    plt.close(fig)
    if replicas:
        fig, axes = plt.subplots(2, 1, figsize=(10, 7))
        for name in replicas:
            table = pd.read_csv(source / name / "residue-rmsf.csv")
            a = table[table.chain == 0]
            axes[0].plot(a.resid, a.rmsf_nm, label=name)
            table = pd.read_csv(source / name / "state.csv")
            axes[1].plot(table.iloc[:, 0] * 0.002, table["Temperature (K)"], label=name)
        axes[0].set(
            xlabel="Author residue, chain A",
            ylabel="CA RMSF, nm",
            title="Short independent apo-dimer trajectories",
        )
        axes[0].legend()
        axes[1].axhline(298.15, color="black", linestyle="--")
        axes[1].set(xlabel="Simulation time, ps", ylabel="Temperature, K")
        fig.tight_layout()
        fig.savefig(output / "replicas.png", dpi=160)
        plt.close(fig)
    coupling = pd.read_csv(source / "torsion-coupling.csv")
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(projection="3d")
    matrix = coupling.pivot(
        index="Tyr14_rotation_deg",
        columns="Asp99_rotation_deg",
        values="coupling_kj_mol",
    )
    xx, yy = np.meshgrid(matrix.columns, matrix.index)
    ax.plot_surface(xx, yy, matrix.values, cmap="coolwarm")
    ax.set(
        xlabel="Asp99 rotation, deg",
        ylabel="Tyr14 rotation, deg",
        zlabel="Coupling, kJ/mol",
        title="Conditional torsion coupling: not chemical PES",
    )
    fig.tight_layout()
    fig.savefig(output / "coupling-surface.png", dpi=160)
    plt.close(fig)
    diary = [
        "# Дневник вычислительных экспериментов",
        "",
        "Время ниже: Europe/Moscow. Полные параметры и UTC — events.jsonl.",
        "",
        "| Время | Эксперимент | Событие | Подробности |",
        "| --- | --- | --- | --- |",
    ]
    events = []
    for journal_path in source.rglob("events.jsonl"):
        for line in journal_path.read_text().splitlines():
            events.append(
                {**json.loads(line), "journal": str(journal_path.relative_to(source))}
            )
    events.sort(key=lambda event: event["utc"])
    (output / "all-events.jsonl").write_text(
        "".join(json.dumps(event, ensure_ascii=False) + "\n" for event in events)
    )
    for event in events:
        time = (
            datetime.fromisoformat(event["utc"])
            .astimezone(ZoneInfo("Europe/Moscow"))
            .strftime("%Y-%m-%d %H:%M:%S")
        )
        detail = event.get(
            "reason",
            event.get(
                "error", f"{event['seconds']:.1f} s" if "seconds" in event else ""
            ),
        )
        diary.append(
            f"| {time} | {event['experiment']} | {event['event']} | {str(detail).replace('|', '/')} |"
        )
    (output / "DIARY.md").write_text("\n".join(diary) + "\n")
    summary = {
        "replicas_exported": replicas,
        "primary_grid": json.loads(
            (source / "H2-conditions-factorial.json").read_text()
        ),
        "torsion_coupling": json.loads(
            (source / "H4-conditional-torsions.json").read_text()
        ),
        "mechanism": json.loads((source / "H5-mechanism-status.json").read_text()),
    }
    if (source / "H3-replicated-md.json").exists():
        summary["replica_tests"] = json.loads(
            (source / "H3-replicated-md.json").read_text()
        )
    if (source / "H2c-grid-refinement.json").exists():
        summary["grid_refinement"] = json.loads(
            (source / "H2c-grid-refinement.json").read_text()
        )
    write_json(output / "summary.json", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.source, args.output)


if __name__ == "__main__":
    main()
