"""Batch physical feature extraction independent of response labels."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re

import pandas as pd

from enzymekinetics_prediction.learning.features import canonical_smiles
from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR
from .artifacts import digest_file, write_json
from .quantum import QuantumSettings, molecular_descriptors
from .structure import pocket_features, protein_id
from .tables import write_table


def collect_pockets(prediction_dir, pdb_dir, output, workers=4):
    if not 1 <= workers <= 8:
        raise ValueError("Pocket extraction supports 1–8 worker threads")
    candidates = {}
    duplicates = 0
    for path in sorted(Path(prediction_dir).rglob("*_predictions.csv")):
        try:
            uid = protein_id(path)
        except ValueError:
            continue
        match = re.search(r"model_v(\d+)", path.name)
        priority = int(match.group(1)) if match else 0
        if uid in candidates:
            duplicates += 1
        if uid not in candidates or priority > candidates[uid][0]:
            candidates[uid] = (priority, path)
    if not candidates:
        raise ValueError("No P2Rank prediction files found")
    failures = []

    def extract(item):
        uid, (_, path) = item
        pdb = Path(pdb_dir) / f"{uid}.pdb"
        try:
            result = pocket_features(path, pdb if pdb.is_file() else None)
            values = {
                f"pocket_{k}": v
                for k, v in result["values"].items()
                if not isinstance(v, list)
            }
            units = {
                f"pocket_{k}": result["units"][k]
                for k in result["values"]
                if not isinstance(result["values"][k], list)
            }
            return {"uniprot": uid, **values}, units, result["inputs"]
        except (
            ValueError,
            OSError,
            pd.errors.ParserError,
            pd.errors.EmptyDataError,
        ) as exc:
            return uid, str(exc), None

    rows, all_units, sources = [], {}, []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for row, units, source in pool.map(extract, sorted(candidates.items())):
            if isinstance(row, dict):
                rows.append(row)
                all_units.update(units)
                sources.append(source)
            else:
                failures.append({"uniprot": row, "error": units})
    if not rows:
        raise ValueError("No usable pockets extracted")
    frame = pd.DataFrame(rows)
    write_table(
        output,
        frame,
        "pocket",
        ["uniprot"],
        all_units,
        {
            "sources": sources,
            "failures": failures,
            "duplicates_resolved_by_highest_version": duplicates,
        },
    )
    return {
        "rows": len(frame),
        "failures": len(failures),
        "missing_composition": int(frame["pocket_charge_proxy"].isna().sum()),
    }


def collect_quantum(
    dataset, output, max_molecules=32, settings=QuantumSettings(), cache_dir=None
):
    frame = pd.read_csv(dataset)
    col = "reactant_smiles" if "reactant_smiles" in frame else "substrate_smiles"
    # Selection uses molecular frequency, never target values.
    counts = Counter(canonical_smiles(str(s)) for s in frame[col].dropna())
    selected = sorted(counts, key=lambda s: (-counts[s], s))
    if max_molecules > 0:
        selected = selected[:max_molecules]
    rows, failures, sources = [], [], []
    units = {
        "qm_energy_hartree": "Hartree",
        "qm_homo_eV": "eV",
        "qm_lumo_eV": "eV",
        "qm_gap_eV": "eV",
        "qm_dipole_norm_D": "Debye",
    }
    if settings.fukui:
        units.update(
            qm_fukui_plus_max="dimensionless", qm_fukui_minus_max="dimensionless"
        )
    for i, smiles in enumerate(selected, 1):
        try:
            result = molecular_descriptors(smiles, settings, cache_dir)
            v = result["values"]
            row = {
                "canonical_smiles": smiles,
                **{
                    f"qm_{k}": v[k]
                    for k in [
                        "energy_hartree",
                        "homo_eV",
                        "lumo_eV",
                        "gap_eV",
                        "dipole_norm_D",
                    ]
                },
            }
            if settings.fukui:
                row.update(
                    qm_fukui_plus_max=max(v["fukui_plus"]),
                    qm_fukui_minus_max=max(v["fukui_minus"]),
                )
            rows.append(row)
            sources.append({"smiles": smiles, "artifact_key": result["cache_key"]})
        except (ValueError, RuntimeError) as exc:
            failures.append({"smiles": smiles, "error": str(exc)})
        if i % 10 == 0:
            print(f"QM {i}/{len(selected)}: successful={len(rows)}", flush=True)
    if not rows:
        raise ValueError("No QM calculations succeeded")
    write_table(
        output,
        pd.DataFrame(rows),
        "qm",
        ["canonical_smiles"],
        units,
        {
            "dataset_sha256": digest_file(dataset),
            "sources": sources,
            "failures": failures,
            "selection": "most frequent unlabeled canonical SMILES",
            "requested": len(selected),
        },
    )
    return {
        "rows": len(rows),
        "failures": len(failures),
        "covered_dataset_rows": sum(counts[r["canonical_smiles"]] for r in rows),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="kind", required=True)
    p = sub.add_parser("pockets")
    p.add_argument(
        "--predictions", type=Path, default=DATA_DIR / "processed/p2rank_out"
    )
    p.add_argument("--pdb-dir", type=Path, default=DATA_DIR / "alphafold_pdb")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--workers", type=int, default=4)
    q = sub.add_parser("quantum")
    q.add_argument("--dataset", type=Path, required=True)
    q.add_argument("--output", type=Path, required=True)
    q.add_argument("--max-molecules", type=int, default=32)
    q.add_argument("--cache-dir", type=Path, default=OUTPUT_DIR / "cache/quantum")
    q.add_argument("--no-fukui", action="store_true")
    args = parser.parse_args()
    if args.kind == "pockets":
        result = collect_pockets(
            args.predictions, args.pdb_dir, args.output, args.workers
        )
    else:
        result = collect_quantum(
            args.dataset,
            args.output,
            args.max_molecules,
            QuantumSettings(fukui=not args.no_fukui),
            args.cache_dir,
        )
    write_json(args.output.with_suffix(".summary.json"), result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
