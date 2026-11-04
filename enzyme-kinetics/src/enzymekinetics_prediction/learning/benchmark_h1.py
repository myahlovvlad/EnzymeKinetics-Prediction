"""Bounded H1 development experiment; delegates all fitting to learning.training."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from . import training
from .embeddings import embed_sequences, embedding_identity, cached_sources
from .features import enzyme_substrate_pairs, sha256
from .predict import predict


def scores(y, prediction):
    result = training.metrics(y, prediction)
    result["spearman"] = (
        float(spearmanr(y, prediction).statistic)
        if len(y) > 2 and np.ptp(y) > 0 and np.ptp(prediction) > 0
        else None
    )
    return result


def cluster_intervals(y, prediction, groups, seed=42, replicates=500):
    """Percentile cluster bootstrap, conditional on the fitted model and split."""
    y, prediction, groups = np.asarray(y), np.asarray(prediction), np.asarray(groups)
    keys = np.unique(groups)
    if len(keys) < 2 or replicates < 20:
        raise ValueError("Need >=2 evaluation clusters and >=20 bootstrap draws")
    rng = np.random.default_rng(seed)
    draws = {k: [] for k in ("mae", "rmse", "r2", "spearman")}
    for _ in range(replicates):
        ids = np.concatenate(
            [np.flatnonzero(groups == g) for g in rng.choice(keys, len(keys))]
        )
        for key, value in scores(y[ids], prediction[ids]).items():
            if key in draws and value is not None and np.isfinite(value):
                draws[key].append(value)
    return {
        k: {
            "ci95": np.quantile(v, [0.025, 0.975]).tolist() if v else None,
            "valid_draws": len(v),
        }
        for k, v in draws.items()
    }


def interval_radius(y, prediction, groups, alpha=0.1):
    """Calibration cluster maxima; infinite radius if finite-sample rank unavailable."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0,1)")
    residual = pd.DataFrame(
        {"group": groups, "residual": np.abs(np.asarray(y) - prediction)}
    )
    values = np.sort(residual.groupby("group").residual.max().to_numpy())
    rank = math.ceil((len(values) + 1) * (1 - alpha))
    return float(values[rank - 1]) if 0 < rank <= len(values) else None


def paired_ablation(y, b0, b1, groups, replicates=500):
    """Paired B1-minus-B0 MAE, resampling the same independent cluster units."""
    difference = np.abs(np.asarray(y) - b1) - np.abs(np.asarray(y) - b0)
    groups = np.asarray(groups)
    keys = np.unique(groups)
    rng = np.random.default_rng(42)
    draws = []
    for _ in range(replicates):
        ids = np.concatenate(
            [np.flatnonzero(groups == g) for g in rng.choice(keys, len(keys))]
        )
        draws.append(float(difference[ids].mean()))
    return {
        "delta_mae": float(difference.mean()),
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "direction": "B1 minus B0; negative favors B1",
    }


def prepare(frame, config):
    """Feature-only eligibility/order; four disjoint provided sequence-cluster groups."""
    group = config["group_by"]
    training.split_groups(frame, group)
    frame = frame.copy()
    frame["original_source_record"] = frame.source_record
    frame = frame[frame.sequence.str.len().between(1, config["max_length"])].copy()
    frame["order"] = frame.sequence.map(
        lambda s: hashlib.sha256(s.encode()).hexdigest()
    )
    frame = (
        frame.sort_values(["order", "source_record"])
        .drop_duplicates("sequence")
        .head(config["max_rows"])
    )
    if len(frame) < 40:
        raise ValueError("Insufficient eligible development records")
    remaining, evaluation = training.grouped_split(
        frame, config["split_seed"], 0.2, group
    )
    pool = frame.iloc[remaining]
    fitting, calibration = training.grouped_split(
        pool, config["split_seed"] + 1, 0.25, group
    )
    parts = {
        "trainval": pool.iloc[fitting].copy(),
        "calibration": pool.iloc[calibration].copy(),
        "development_evaluation": frame.iloc[evaluation].copy(),
    }
    return parts


def leakage_audit(parts, group):
    def keys(frame, key):
        if key == "pair":
            return set(enzyme_substrate_pairs(frame))
        if key == "scaffold":
            column = (
                "reactant_smiles" if "reactant_smiles" in frame else "substrate_smiles"
            )
            values = set()
            for smiles in frame[column].fillna(""):
                molecule = Chem.MolFromSmiles(str(smiles))
                if molecule is not None:
                    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
                    if scaffold:
                        values.add(scaffold)
            return values
        return set(frame[key].dropna().astype(str))

    result = {
        "pretraining_overlap": "UNKNOWN: complete training membership not available",
        "sequence_identity": "provided cluster IDs only; not independently recomputed",
        "scaffold_empty": "acyclic/invalid molecules excluded from scaffold-set counts",
    }
    names = list(parts)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            overlaps = {
                k: len(keys(parts[a], k) & keys(parts[b], k))
                for k in ["sequence", group, "pair", "ec", "scaffold"]
            }
            if any(overlaps[k] for k in ("sequence", group, "pair")):
                raise ValueError("Development partitions leak sequence/cluster/pair")
            result[a + ":" + b] = overlaps
    for name, frame in parts.items():
        result[name] = {
            "rows": len(frame),
            "clusters": int(frame[group].nunique()),
            "missing_ph": int(frame.ph.isna().sum()),
            "missing_temperature": int(frame.temperature.isna().sum()),
        }
    return result


def run(config, data_dir, output):
    output, data_dir = Path(output), Path(data_dir)
    if output.exists():
        raise FileExistsError(output)
    if (
        config["encoder"] != "facebook/esm2_t6_8M_UR50D"
        or len(config["seeds"]) < 2
        or len(set(config["seeds"])) != len(config["seeds"])
        or not set(config["targets"]) <= {"kcat", "km"}
        or not 1 <= config["epochs"] <= 50
        or not 40 <= config["max_rows"]
        or not 1 <= config["threads"]
        or not 1 <= config["max_length"]
        or config["variants"] != ["B0", "B1"]
        or config["max_rows"] > 256
        or config["max_length"] > 256
        or config["threads"] > 2
    ):
        raise ValueError("This pilot permits B0/B1, <=256 rows/residues, <=2 threads")
    output.mkdir(parents=True)
    training.write_json(output / "config.json", config)
    started = time.monotonic()
    results = {}
    for target in config["targets"]:
        source = data_dir / "Baseline" / target / f"{target}-random_trainval.csv"
        frame, dropped = training.load_table(source)
        parts = prepare(frame, config)
        target_dir = output / target
        prepared = target_dir / "prepared"
        destination = prepared / "Baseline" / target
        destination.mkdir(parents=True)
        feature_columns = [
            c
            for c in [
                "sequence",
                "uniprot",
                "ec",
                "ph",
                "temperature",
                "reactant_smiles",
                "substrate_smiles",
            ]
            if c in frame
        ]
        parts["trainval"][[*feature_columns, "log10_value", config["group_by"]]].to_csv(
            destination / source.name, index=False
        )
        # Save the complete explicit partition, including source row lineage.
        fitting, validation = training.grouped_split(
            parts["trainval"], config["split_seed"], 0.25, config["group_by"]
        )
        audit_parts = {
            "train": parts["trainval"].iloc[fitting],
            "validation": parts["trainval"].iloc[validation],
            **{k: v for k, v in parts.items() if k != "trainval"},
        }
        partition = pd.concat([v.assign(partition=k) for k, v in audit_parts.items()])
        partition[
            ["original_source_record", "sequence", config["group_by"], "partition"]
        ].to_csv(target_dir / "partition.csv", index=False)
        audit = leakage_audit(audit_parts, config["group_by"])
        training.write_json(target_dir / "audit.json", audit)
        embedding = target_dir / "esm2.pt"
        encoder_seconds = None
        case_results = {}
        ensemble_predictions = {}
        for variant in config["variants"]:
            cache_identity = None
            if variant == "B1":
                encoding_started = time.monotonic()
                embed_sequences(
                    partition.sequence,
                    embedding,
                    config["encoder"],
                    config["revision"],
                    batch_size=2,
                    threads=config["threads"],
                    local_files_only=True,
                )
                encoder_seconds = time.monotonic() - encoding_started
                cache_identity = embedding_identity(
                    partition.sequence,
                    config["encoder"],
                    config["revision"],
                    cached_sources(config["encoder"], config["revision"]),
                )
            ensemble_eval, ensemble_cal = [], []
            for seed in config["seeds"]:
                directory = target_dir / variant / str(seed)
                args = SimpleNamespace(
                    target=target,
                    output=directory,
                    data_dir=prepared,
                    seed=seed,
                    split_seed=config["split_seed"],
                    validation_fraction=0.25,
                    group_by=config["group_by"],
                    ec_level=2,
                    split_scheme="provided-test",
                    development_only=True,
                    models=config["heads"],
                    embedding_path=embedding
                    if variant == "B1"
                    else prepared / "absent-embedding.pt",
                    physical_table=[],
                    embedding_cache_identity=cache_identity,
                    embedding_dimension=320 if cache_identity is not None else None,
                    epochs=config["epochs"],
                    patience=5,
                    batch_size=64,
                    threads=config["threads"],
                )
                t0 = time.monotonic()
                training.run(args)
                # No labels cross the inference API boundary.
                cal = predict(
                    directory,
                    parts["calibration"][feature_columns],
                    data_dir=prepared,
                )
                evaluation, flags = predict(
                    directory,
                    parts["development_evaluation"][feature_columns],
                    data_dir=prepared,
                    return_flags=True,
                )
                if variant == "B0" and seed == config["seeds"][0]:
                    repeat_args = SimpleNamespace(**vars(args))
                    repeat_args.output = target_dir / "B0-reproducibility"
                    training.run(repeat_args)
                    repeated = predict(
                        repeat_args.output,
                        parts["development_evaluation"][feature_columns],
                        data_dir=prepared,
                    )
                    np.testing.assert_array_equal(evaluation, repeated)
                    training.write_json(
                        target_dir / "reproducibility.json",
                        {
                            "seed": seed,
                            "prediction_max_absolute_difference": float(
                                np.max(np.abs(evaluation - repeated))
                            ),
                            "same_seed_retrained": True,
                        },
                    )
                ensemble_eval.append(evaluation)
                ensemble_cal.append(cal)
                radius = interval_radius(
                    parts["calibration"].log10_value.to_numpy(),
                    cal,
                    parts["calibration"][config["group_by"]],
                )
                y = parts["development_evaluation"].log10_value.to_numpy()
                metrics = scores(y, evaluation)
                metrics.update(
                    ci=cluster_intervals(
                        y,
                        evaluation,
                        parts["development_evaluation"][config["group_by"]],
                        replicates=config["bootstrap_draws"],
                    ),
                    interval_radius_log10=radius,
                    empirical_interval_coverage=float(
                        np.mean(np.abs(y - evaluation) <= radius)
                    )
                    if radius is not None
                    else None,
                    runtime_seconds=time.monotonic() - t0,
                    process_peak_rss_MiB=resource.getrusage(
                        resource.RUSAGE_SELF
                    ).ru_maxrss
                    / 1024,
                    missing_embedding_fraction=float(flags[:, 1].mean()),
                    selected_model=json.loads(
                        (directory / "metadata.json").read_text()
                    )["selected_model"],
                )
                case_results[variant + ":" + str(seed)] = metrics
            pred = np.mean(ensemble_eval, axis=0)
            ensemble_predictions[variant] = pred
            cal = np.mean(ensemble_cal, axis=0)
            y = parts["development_evaluation"].log10_value.to_numpy()
            radius = interval_radius(
                parts["calibration"].log10_value.to_numpy(),
                cal,
                parts["calibration"][config["group_by"]],
            )
            case_results[variant + ":ensemble"] = {
                **scores(y, pred),
                "ci": cluster_intervals(
                    y,
                    pred,
                    parts["development_evaluation"][config["group_by"]],
                    replicates=config["bootstrap_draws"],
                ),
                "interval_radius_log10": radius,
                "empirical_interval_coverage": float(
                    np.mean(np.abs(y - pred) <= radius)
                )
                if radius is not None
                else None,
                "interval_width_log10": 2 * radius if radius is not None else None,
                "cluster_interval_coverage": float(
                    pd.DataFrame(
                        {
                            "cluster": parts["development_evaluation"][
                                config["group_by"]
                            ],
                            "covered": np.abs(y - pred) <= radius,
                        }
                    )
                    .groupby("cluster")
                    .covered.min()
                    .mean()
                )
                if radius is not None
                else None,
            }
            pd.DataFrame(
                {
                    "source_record": parts[
                        "development_evaluation"
                    ].original_source_record,
                    "cluster": parts["development_evaluation"][config["group_by"]],
                    "observed_log10": y,
                    "prediction_log10": pred,
                }
            ).to_csv(
                target_dir / (variant + "-development-predictions.csv"), index=False
            )
        results[target] = {
            "source_sha256": sha256(source),
            "source_rows": len(frame),
            "dropped_rows": dropped,
            "encoder_seconds": encoder_seconds,
            "cases": case_results,
            "paired_embedding_ablation": paired_ablation(
                y,
                ensemble_predictions["B0"],
                ensemble_predictions["B1"],
                parts["development_evaluation"][config["group_by"]],
                config["bootstrap_draws"],
            ),
            "unit": "s^-1 (upstream convention)"
            if target == "kcat"
            else "dataset-native concentration; assay units unverified",
        }
        training.write_json(output / "results.json", results)
    training.write_json(
        output / "manifest.json",
        {
            "scope": "development_pilot_only",
            "official_test_accessed": False,
            "code_sha256": {
                p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")
            },
            "output_sha256": {
                str(p.relative_to(output)): sha256(p)
                for p in sorted(output.rglob("*"))
                if p.is_file()
            },
            "elapsed_seconds": time.monotonic() - started,
            "peak_rss_scope": "process lifetime high-water mark, not per-case allocation",
            "scientific_status": "exploratory_not_validated",
        },
    )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    with threadpool_limits(limits=config["threads"]):
        run(config, args.data_dir, args.output)


if __name__ == "__main__":
    main()
