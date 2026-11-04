"""Verify frozen-model reload against recorded predictions, using original source rows."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from enzymekinetics_prediction.learning.predict import predict
from enzymekinetics_prediction.paths import DATA_DIR


def verify(run, data_dir, records=12):
    meta = json.loads((run / "metadata.json").read_text())
    saved = pd.read_csv(run / "predictions.csv").head(records)
    if saved.empty:
        raise ValueError("No predictions to verify")
    rows, cache = [], {}
    for _, row in saved.iterrows():
        filename = row.get("source_dataset", meta["target"] + "-random_test.csv")
        path = data_dir / "Baseline" / meta["target"] / filename
        if path not in cache:
            cache[path] = pd.read_csv(path)
        # Training source_record is a one-based DATA row, excluding the header.
        rows.append(cache[path].iloc[int(row.source_record) - 1])
    frame = pd.DataFrame(rows).drop(columns="log10_value")
    errors = {}
    for model in meta["models"]:
        actual = predict(run, frame, model, data_dir)
        expected = saved[model].to_numpy()
        np.testing.assert_allclose(
            actual, expected, atol=1e-5, rtol=1e-5, err_msg=f"{run.name}/{model}"
        )
        errors[model] = float(np.max(np.abs(actual - expected)))
    return {
        "run": run.name,
        "records": len(saved),
        "max_absolute_error": errors,
        "status": "passed",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--records", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.records < 1:
        raise ValueError("Positive record count required")
    results = []
    for run in args.runs:
        results.append(verify(run, args.data_dir, args.records))
        print(json.dumps(results[-1]), flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
