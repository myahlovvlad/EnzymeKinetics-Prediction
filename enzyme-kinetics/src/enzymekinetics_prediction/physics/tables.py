"""Physical feature tables with explicit join keys, units and checksums."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .artifacts import digest_file, write_json

ALLOWED_JOINS = {"uniprot", "canonical_smiles", "ph", "temperature", "ionic_strength_M"}
BLOCKS = {"pocket", "qm", "electro", "md"}


def write_table(path, frame, block, join_on, units, provenance):
    path = Path(path)
    if (
        block not in BLOCKS
        or not set(join_on) <= ALLOWED_JOINS
        or not join_on
        or len(set(join_on)) != len(join_on)
        or not set(join_on) <= set(frame.columns)
    ):
        raise ValueError("Invalid physical table block or keys")
    columns = [c for c in frame if c not in join_on]
    if set(columns) != set(units) or any(
        not c.startswith(block + "_") for c in columns
    ):
        raise ValueError(
            "Only explicitly named physical features with units are allowed"
        )
    if frame.duplicated(join_on).any():
        raise ValueError("Physical table must have unique join keys")
    for c in columns:
        numeric = pd.to_numeric(frame[c], errors="raise")
        if np.isinf(numeric).any():
            raise ValueError("Infinite physical feature")
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    write_json(
        path.with_suffix(path.suffix + ".meta.json"),
        {
            "schema_version": 1,
            "block": block,
            "join_on": list(join_on),
            "feature_columns": columns,
            "units": units,
            "csv_sha256": digest_file(path),
            "generated_without_labels": True,
            "provenance": provenance,
        },
    )
    return path


def key_for(row, join_on):
    values = []
    for name in join_on:
        value = row.get(name)
        if pd.isna(value):
            values.append(None)
        elif name in {"ph", "temperature", "ionic_strength_M"}:
            values.append(float(value))
        else:
            values.append(str(value).strip())
    return tuple(values)


def load_tables(paths):
    tables = []
    for path in paths or []:
        path = Path(path)
        meta = json.loads(path.with_suffix(path.suffix + ".meta.json").read_text())
        if (
            meta.get("schema_version") != 1
            or meta.get("block") not in BLOCKS
            or not meta.get("generated_without_labels")
        ):
            raise ValueError("Unsupported or label-dependent physical table")
        if meta["csv_sha256"] != digest_file(path):
            raise ValueError("Physical table checksum mismatch")
        joins, features = meta["join_on"], meta["feature_columns"]
        if (
            not joins
            or len(set(joins)) != len(joins)
            or not set(joins) <= ALLOWED_JOINS
            or not features
            or len(set(features)) != len(features)
            or set(features) != set(meta["units"])
        ):
            raise ValueError("Invalid table keys or units")
        if any(not c.startswith(meta["block"] + "_") for c in features):
            raise ValueError("Invalid physical feature name")
        frame = pd.read_csv(path)
        if set(frame.columns) != set(joins + features) or frame.duplicated(joins).any():
            raise ValueError("Unexpected columns or duplicate physical join keys")
        matrix = frame[features].to_numpy(dtype=float)
        if np.isinf(matrix).any():
            raise ValueError("Infinite physical table values")
        lookup = {
            key_for(row, joins): matrix[i].astype(np.float32)
            for i, (_, row) in enumerate(frame.iterrows())
        }
        tables.append(
            {
                "metadata": meta,
                "lookup": lookup,
                "feature_count": len(features),
                "path": str(path),
            }
        )
    return tables
