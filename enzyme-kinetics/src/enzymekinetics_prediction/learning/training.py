"""Train, select on grouped validation, then evaluate frozen models once on test."""

import argparse
import copy
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import random
import re
import time

import joblib
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
import torch
from torch import nn

from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR
from .features import FEATURE_VERSION, FeatureBuilder, enzyme_substrate_pairs, sha256
from .models import KineticsMLP, PhysicsAdd3DGraphNet
from .chemometrics import ConsensusMBPLS, OPLSRegressor, ECLocalRegressor


def metrics(y, predicted):
    if not len(y):
        return {"n": 0, "mae": None, "rmse": None, "r2": None, "within_factor_10": None}
    return {
        "n": len(y),
        "mae": float(mean_absolute_error(y, predicted)),
        "rmse": float(np.sqrt(mean_squared_error(y, predicted))),
        "r2": float(r2_score(y, predicted)) if len(y) > 1 else None,
        "within_factor_10": float(np.mean(np.abs(y - predicted) <= 1)),
    }


def split_groups(frame, group_by="sequence", ec_level=2):
    if group_by == "ec":
        return (
            frame["ec"]
            .fillna("unknown")
            .astype(str)
            .map(lambda v: ".".join(v.split(".")[:ec_level]))
        )
    if group_by not in frame:
        raise ValueError(f"Split-group column is missing: {group_by}")
    values = frame[group_by]
    if values.isna().any() or values.astype(str).str.strip().eq("").any():
        raise ValueError("Split groups cannot be missing")
    return values.astype(str).str.strip()


def grouped_split(
    frame, seed=42, validation_fraction=0.2, group_by="sequence", ec_level=2
):
    groups = split_groups(frame, group_by, ec_level)
    if groups.nunique() < 2:
        raise ValueError("At least two independent groups are needed")
    splitter = GroupShuffleSplit(
        n_splits=1, test_size=validation_fraction, random_state=seed
    )
    train, validation = next(splitter.split(frame, groups=groups))
    if set(groups.iloc[train]) & set(groups.iloc[validation]):
        raise RuntimeError("Training and validation groups overlap")
    return train, validation


def model_predict(model, x, frame):
    return (
        model.predict(x, ec=frame.ec)
        if isinstance(model, ECLocalRegressor)
        else model.predict(x).reshape(-1)
    )


def parameter_record(value):
    if isinstance(value, dict):
        return {str(k): parameter_record(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [parameter_record(v) for v in value]
    if hasattr(value, "get_params"):
        return {
            "class": type(value).__name__,
            "parameters": parameter_record(value.get_params(deep=False)),
        }
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def load_table(path):
    frame = pd.read_csv(path)
    required = {"sequence", "uniprot", "ec", "log10_value"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    if not {"reactant_smiles", "substrate_smiles"} & set(frame.columns):
        raise ValueError(f"No substrate SMILES column in {path}")
    frame["source_dataset"] = Path(path).name
    frame["source_record"] = np.arange(1, len(frame) + 1)
    y = pd.to_numeric(frame["log10_value"], errors="coerce")
    valid = (
        np.isfinite(y)
        & frame["sequence"].notna()
        & frame["sequence"].astype(str).str.strip().ne("")
    )
    dropped = int((~valid).sum())
    frame = frame.loc[valid].copy().reset_index(drop=True)
    frame["sequence"] = frame["sequence"].astype(str).str.strip()
    frame["log10_value"] = y.loc[valid].to_numpy()
    if len(frame) < 2:
        raise ValueError(f"Not enough valid records: {path}")
    return frame, dropped


def read_embeddings(path, expected_identity=None, expected_dimension=None):
    strict_cache_key = None
    if expected_identity is not None:
        from .embeddings import validate_embedding_cache

        strict_cache_key = validate_embedding_cache(
            path, expected_identity, expected_dimension
        )
    if not path.is_file():
        return {}, {"status": "absent", "path": str(path)}
    # Never silently retry with unrestricted pickle loading.
    raw = torch.load(path, map_location="cpu", weights_only=True)
    if not isinstance(raw, dict) or not raw:
        raise ValueError("Embeddings must be a nonempty dictionary of tensor vectors")
    vectors = {}
    collisions = 0
    # Prefer AlphaFold over ESMFold when both representations exist for one ID.
    for key in sorted(raw, key=lambda k: (not str(k).startswith("AF-"), str(k))):
        value = raw[key]
        if not isinstance(value, torch.Tensor):
            raise ValueError("Unsupported embeddings format: expected tensor values")
        vector = value.detach().cpu().numpy().astype(np.float32).reshape(-1)
        if not np.isfinite(vector).all():
            raise ValueError("Nonfinite embedding")
        vectors[str(key)] = vector
        match = re.match(r"^(?:ESMFold-|AF-)([A-Z0-9]+)(?:-|\.)", str(key))
        if match:
            uid = match.group(1)
            if uid in vectors:
                collisions += 1
            else:
                vectors[uid] = vector
    if len({len(v) for v in vectors.values()}) != 1:
        raise ValueError("Inconsistent embedding lengths")
    return vectors, {
        "status": "loaded",
        "path": str(path),
        "sha256": sha256(path),
        "entries": len(raw),
        "dimension": len(next(iter(vectors.values()))),
        **(
            {"strict_cache_key": strict_cache_key}
            if strict_cache_key is not None
            else {}
        ),
        "duplicate_uniprot_representations": collisions,
        "id_policy": "parse AF-/ESMFold- filename; prefer AF when both exist",
    }


def neural_predict(model, x, geometry):
    model.eval()
    with torch.no_grad():
        return np.concatenate(
            [
                model(
                    torch.from_numpy(x[i : i + 512]),
                    torch.from_numpy(geometry[i : i + 512]),
                ).numpy()
                for i in range(0, len(x), 512)
            ]
        )


def train_neural(
    name, x, geometry, y, xv, gv, yv, *, seed, epochs, patience, batch_size
):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    model = (
        KineticsMLP(x.shape[1]) if name == "mlp" else PhysicsAdd3DGraphNet(x.shape[1])
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
    mean, scale = float(y.mean()), max(float(y.std()), 1e-6)
    tx, tg, ty = (
        torch.from_numpy(x),
        torch.from_numpy(geometry),
        torch.from_numpy((y - mean) / scale),
    )
    generator = torch.Generator().manual_seed(seed)
    history, best, best_state, stalled = [], float("inf"), None, 0
    best_epoch = 0
    for epoch in range(1, epochs + 1):
        model.train()
        permutation = torch.randperm(len(x), generator=generator)
        train_loss = 0.0
        for ids in permutation.split(batch_size):
            optimizer.zero_grad()
            loss = nn.functional.mse_loss(model(tx[ids], tg[ids]), ty[ids])
            if not torch.isfinite(loss):
                raise ValueError("Nonfinite neural loss")
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            train_loss += float(loss.detach()) * len(ids)
        predicted = neural_predict(model, xv, gv) * scale + mean
        score = float(mean_absolute_error(yv, predicted))
        history.append(
            {
                "epoch": epoch,
                "train_standardized_mse": train_loss / len(x),
                "validation_mae_log10": score,
            }
        )
        if score < best - 1e-5:
            best, best_state, best_epoch, stalled = (
                score,
                copy.deepcopy(model.state_dict()),
                epoch,
                0,
            )
        else:
            stalled += 1
        if epoch == 1 or epoch % 10 == 0:
            print(f"{name}: epoch={epoch} validation_MAE={score:.4f}", flush=True)
        if stalled >= patience:
            break
    model.load_state_dict(best_state)
    return model, {
        "name": name,
        "input_dim": x.shape[1],
        "target_mean": mean,
        "target_scale": scale,
        "best_epoch": best_epoch,
        "epochs_run": len(history),
        "history": history,
    }


def write_json(path, value):
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    )


def run(args):
    package_root = Path(__file__).parent.parent
    code_hashes = {
        str(p.relative_to(package_root)): sha256(p) for p in package_root.rglob("*.py")
    }
    if args.epochs < 1 or args.threads < 1 or args.batch_size < 1 or args.patience < 1:
        raise ValueError("Epochs, threads, patience and batch size must be positive")
    torch.set_num_threads(args.threads)
    torch.use_deterministic_algorithms(True)
    started = time.monotonic()
    output = args.output or OUTPUT_DIR / "training" / (
        args.target + "-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    )
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to overwrite nonempty run directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    train_path = (
        args.data_dir / "Baseline" / args.target / f"{args.target}-random_trainval.csv"
    )
    test_path = (
        args.data_dir / "Baseline" / args.target / f"{args.target}-random_test.csv"
    )
    frame, dropped_train = load_table(train_path)
    development_only = getattr(args, "development_only", False)
    if (
        development_only
        and getattr(args, "split_scheme", "provided-test") != "provided-test"
    ):
        raise ValueError("Development-only mode cannot repartition the official test")
    test, dropped_test = (None, 0) if development_only else load_table(test_path)
    split_seed = getattr(args, "split_seed", None)
    split_seed = args.seed if split_seed is None else split_seed
    group_by = getattr(args, "group_by", "sequence")
    ec_level = getattr(args, "ec_level", 2)
    split_scheme = getattr(args, "split_scheme", "provided-test")
    if split_scheme == "grouped-three-way":
        combined = pd.concat([frame, test], ignore_index=True)
        remain_ids, test_ids = grouped_split(
            combined, split_seed, 0.2, group_by, ec_level
        )
        frame, test = (
            combined.iloc[remain_ids].reset_index(drop=True),
            combined.iloc[test_ids].reset_index(drop=True),
        )
    train_ids, validation_ids = grouped_split(
        frame,
        split_seed + 1 if split_scheme == "grouped-three-way" else split_seed,
        args.validation_fraction,
        group_by,
        ec_level,
    )
    fit_frame, val_frame = frame.iloc[train_ids], frame.iloc[validation_ids]
    partitions = [("train", fit_frame), ("validation", val_frame)]
    if not development_only:
        partitions.append(("test", test))
    split_frame = pd.concat(
        [records.assign(partition=name) for name, records in partitions]
    )
    split_frame["split_group"] = split_groups(split_frame, group_by, ec_level)
    split_frame[
        [
            "source_dataset",
            "source_record",
            "partition",
            "uniprot",
            "sequence",
            "split_group",
        ]
    ].to_csv(output / "split.csv", index=False)
    embeddings, embedding_info = read_embeddings(
        getattr(args, "embedding_path", None)
        or args.data_dir / "CatPred-DB/data/catpred_progres_embeds_dict.pt",
        getattr(args, "embedding_cache_identity", None),
        getattr(args, "embedding_dimension", None),
    )
    builder = FeatureBuilder(
        args.data_dir, embeddings, getattr(args, "physical_table", None)
    ).fit(fit_frame)
    print("Building train, validation and test features", flush=True)
    matrices, geometries, coverage = {}, {}, {}
    physical_coverage = {}
    for name, records in partitions:
        matrices[name] = builder.transform(records)
        geometries[name] = builder.last_geometry_.copy()
        coverage[name] = dict(
            zip(
                [
                    "invalid_smiles_fraction",
                    "missing_embedding_fraction",
                    "missing_structure_fraction",
                ],
                map(float, builder.last_flags_.mean(axis=0)),
                strict=True,
            )
        )
        physical_coverage[name] = builder.last_physical_coverage_.copy()
        print(f"{name}: {matrices[name].shape}, coverage={coverage[name]}", flush=True)
    preprocessing = make_pipeline(
        SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()
    )
    preprocessing.fit(matrices["train"])
    matrices = {
        name: np.clip(preprocessing.transform(value), -10, 10).astype(np.float32)
        for name, value in matrices.items()
    }
    x, xv = matrices["train"], matrices["validation"]
    g, gv = geometries["train"], geometries["validation"]
    y = fit_frame.log10_value.to_numpy(dtype=np.float32)
    yv = val_frame.log10_value.to_numpy(dtype=np.float32)
    results, fitted, specs, validation_predictions = {}, {}, {}, {}
    for name in args.models:
        t0 = time.monotonic()
        print(f"Training {name}", flush=True)
        if name in {"mlp", "radial"}:
            model, spec = train_neural(
                name,
                x,
                g,
                y,
                xv,
                gv,
                yv,
                seed=args.seed,
                epochs=args.epochs,
                patience=args.patience,
                batch_size=args.batch_size,
            )
            pred = (
                neural_predict(model, xv, gv) * spec["target_scale"]
                + spec["target_mean"]
            )
            torch.save(
                {
                    "state_dict": model.state_dict(),
                    "spec": {k: v for k, v in spec.items() if k != "history"},
                },
                output / f"{name}.pt",
            )
            write_json(output / f"{name}_history.json", spec["history"])
        else:
            if name == "mean":
                candidates = [DummyRegressor(strategy="mean")]
            elif name == "ridge":
                candidates = [Ridge(alpha=a) for a in (1.0, 10.0, 100.0)]
            elif name == "pls":
                candidates = [
                    PLSRegression(
                        n_components=min(n, x.shape[1], len(x) - 1), scale=False
                    )
                    for n in (8, 16)
                ]
            elif name == "mbpls":
                candidates = [
                    ConsensusMBPLS(builder.feature_blocks_, n_components=n)
                    for n in (8, 16)
                ]
            elif name == "opls":
                candidates = [OPLSRegressor(n_orthogonal=n) for n in (0, 1, 2)]
            elif name == "ec_mbpls":
                candidates = [
                    ECLocalRegressor(
                        ConsensusMBPLS(builder.feature_blocks_, n_components=8),
                        min_group_size=100,
                        ec_level=ec_level,
                    )
                ]
            elif name == "hist":
                candidates = [
                    HistGradientBoostingRegressor(
                        max_iter=150,
                        max_leaf_nodes=15,
                        l2_regularization=10.0,
                        early_stopping=False,
                        random_state=args.seed,
                    )
                ]
            else:
                raise ValueError(name)
            candidate_scores = []
            for candidate in candidates:
                if isinstance(candidate, ECLocalRegressor):
                    candidate.fit(x, y, ec=fit_frame.ec)
                else:
                    candidate.fit(x, y)
                candidate_scores.append(
                    float(
                        mean_absolute_error(yv, model_predict(candidate, xv, val_frame))
                    )
                )
            model = candidates[int(np.argmin(candidate_scores))]
            pred = model_predict(model, xv, val_frame)
            spec = {
                "name": name,
                "selected_parameters": parameter_record(model.get_params(deep=False)),
                "candidate_validation_mae": candidate_scores,
            }
            joblib.dump(model, output / f"{name}.joblib")
        fitted[name], specs[name] = model, spec
        validation_predictions[name] = pred
        results[name] = {
            "validation": metrics(yv, pred),
            "training_seconds": time.monotonic() - t0,
        }
        print(
            f"{name}: validation_MAE={results[name]['validation']['mae']:.4f}",
            flush=True,
        )
    selected = min(results, key=lambda name: results[name]["validation"]["mae"])
    if development_only:
        # A separate opt-in research boundary: never open/hash the official test.
        builder.last_geometry_ = None
        builder.last_flags_ = None
        joblib.dump(
            {"builder": builder, "preprocessing": preprocessing},
            output / "preprocessing.joblib",
        )
        predictions = val_frame[
            ["source_dataset", "source_record", "log10_value"]
        ].copy()
        for name, values in validation_predictions.items():
            predictions[name] = values
        predictions.to_csv(output / "validation-predictions.csv", index=False)
        write_json(
            output / "metadata.json",
            {
                "target": args.target,
                "selected_model": selected,
                "models": specs,
                "evaluation_scope": "development_only",
                "official_test_accessed": False,
                "seed": args.seed,
                "split_seed": split_seed,
                "group_by": group_by,
                "trainval_sha256": sha256(train_path),
                "code_sha256_at_start": code_hashes,
                "coverage": coverage,
                "embeddings": embedding_info,
                "feature_blocks": builder.feature_blocks_,
                "rows": {"train": len(y), "validation": len(yv)},
                "versions": {
                    p: importlib.metadata.version(p)
                    for p in ["torch", "numpy", "pandas", "rdkit", "scikit-learn"]
                },
                "elapsed_seconds": time.monotonic() - started,
                "limitations": [
                    "Development validation; not an official test estimate"
                ],
            },
        )
        write_json(output / "metrics.json", results)
        return output
    # Selection is frozen before any test targets are evaluated.
    print(
        f"Selected by validation: {selected}; evaluating frozen models on test",
        flush=True,
    )
    y_test = test.log10_value.to_numpy(dtype=np.float32)
    novel = ~test.sequence.isin(set(frame.sequence)).to_numpy()
    trainval_pairs = set(enzyme_substrate_pairs(frame))
    pair_disjoint = np.asarray(
        [pair not in trainval_pairs for pair in enzyme_substrate_pairs(test)]
    )
    predictions = test[
        ["source_dataset", "source_record", "uniprot", "log10_value"]
    ].copy()
    predictions["unseen_exact_sequence_vs_trainval"] = novel
    predictions["canonical_pair_disjoint_vs_trainval"] = pair_disjoint
    for name, model in fitted.items():
        if name in {"mlp", "radial"}:
            pred = neural_predict(model, matrices["test"], geometries["test"])
            pred = pred * specs[name]["target_scale"] + specs[name]["target_mean"]
        else:
            pred = model_predict(model, matrices["test"], test)
        results[name]["test"] = metrics(y_test, pred)
        results[name]["test_unseen_exact_sequence"] = metrics(
            y_test[novel], pred[novel]
        )
        results[name]["test_pair_disjoint"] = metrics(
            y_test[pair_disjoint], pred[pair_disjoint]
        )
        predictions[name] = pred
    predictions.to_csv(output / "predictions.csv", index=False)
    builder.last_geometry_ = None
    builder.last_flags_ = None
    joblib.dump(
        {"builder": builder, "preprocessing": preprocessing},
        output / "preprocessing.joblib",
    )
    write_json(output / "geometry_provenance.json", builder.geometry_provenance_)
    metadata = {
        "target": args.target,
        "target_column": "log10_value",
        "target_scale": "log10",
        "target_unit": "s^-1"
        if args.target == "kcat"
        else "dataset-native concentration unit; verify assay provenance",
        "selected_model": selected,
        "selection_metric": f"validation MAE, group={group_by}",
        "split_scheme": split_scheme,
        "group_by": group_by,
        "ec_level": ec_level,
        "test_group_overlap_rows": int(
            split_groups(test, group_by, ec_level)
            .isin(set(split_groups(frame, group_by, ec_level)))
            .sum()
        ),
        "feature_version": FEATURE_VERSION,
        "code_sha256_at_start": code_hashes,
        "canonical_pair_overlap_test_rows": int((~pair_disjoint).sum()),
        "neural_training": {
            "max_epochs": args.epochs,
            "patience": args.patience,
            "batch_size": args.batch_size,
            "optimizer": "AdamW",
            "learning_rate": 1e-3,
            "weight_decay": 1e-3,
            "target_standardization": "mean/std fitted on train only",
        },
        "seed": args.seed,
        "threads": args.threads,
        "models": specs,
        "coverage": coverage,
        "physical_coverage": physical_coverage,
        "feature_blocks": builder.feature_blocks_,
        "physical_tables": [t["metadata"] for t in builder.physical_lookup_],
        "embeddings": embedding_info,
        "trainval_sha256": sha256(train_path),
        "test_sha256": sha256(test_path),
        "trainval_path": str(train_path),
        "test_path": str(test_path),
        "rows": {"train": len(y), "validation": len(yv), "test": len(test)},
        "dropped_rows": {"trainval": dropped_train, "test": dropped_test},
        "limitations": [
            "Single-seed exploratory benchmark",
            "Exact-sequence novelty is not sequence-cluster OOD",
            "Radial geometry is not a molecular interaction graph or a QM barrier",
            "No empirical claim of reaction-mechanism simulation",
            "Products and target-derived barriers excluded",
        ],
        "versions": {
            p: importlib.metadata.version(p)
            for p in ["torch", "numpy", "pandas", "rdkit", "scikit-learn"]
        },
        "elapsed_seconds": time.monotonic() - started,
    }
    write_json(output / "metadata.json", metadata)
    write_json(output / "metrics.json", results)
    lines = [
        f"# {args.target}: completed training benchmark",
        "",
        f"Selected by validation: **{selected}**.",
        "",
        "| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, r in results.items():
        novel_mae = r["test_unseen_exact_sequence"]["mae"]
        lines.append(
            f"| {name} | {r['validation']['mae']:.4f} | {r['test']['mae']:.4f} | {r['test']['rmse']:.4f} | {r['test']['r2']:.4f} | {novel_mae if novel_mae is not None else 'n/a'} |"
        )
    lines += [
        "",
        "All errors are in log10 target units. Test was not used for model selection.",
        "The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.",
        "Single seed; no claim of calibrated uncertainty or superiority to published models.",
        "",
        f"Rows: {metadata['rows']}. Wall time: {metadata['elapsed_seconds']:.1f} s.",
    ]
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Saved run: {output}", flush=True)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=["kcat", "km", "ki"], default="kcat")
    parser.add_argument(
        "--models",
        nargs="+",
        choices=[
            "mean",
            "ridge",
            "pls",
            "hist",
            "mlp",
            "radial",
            "mbpls",
            "opls",
            "ec_mbpls",
        ],
        default=["mean", "ridge", "pls", "hist", "mlp", "radial"],
    )
    parser.add_argument(
        "--group-by",
        choices=["sequence", "ec", "sequence_40cluster", "sequence_60cluster"],
        default="sequence",
    )
    parser.add_argument("--ec-level", type=int, choices=[1, 2, 3, 4], default=2)
    parser.add_argument(
        "--split-scheme",
        choices=["provided-test", "grouped-three-way"],
        default="provided-test",
    )
    parser.add_argument("--physical-table", action="append", type=Path)
    parser.add_argument("--embedding-path", type=Path)
    parser.add_argument(
        "--development-only",
        action="store_true",
        help="Train/select using trainval only; do not open the official test",
    )
    parser.add_argument(
        "--split-seed", type=int, help="Fix partitions independently of the model seed"
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--threads", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    args = parser.parse_args()
    with threadpool_limits(limits=args.threads):
        run(args)


if __name__ == "__main__":
    main()
