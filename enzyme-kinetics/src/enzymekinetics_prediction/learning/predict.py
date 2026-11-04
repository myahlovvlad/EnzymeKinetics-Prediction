"""Load a trusted local training run and predict log10 kinetics from CSV."""

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

from enzymekinetics_prediction.paths import DATA_DIR
from .models import KineticsMLP, PhysicsAdd3DGraphNet
from .training import neural_predict, model_predict


def predict(run_dir, frame, model_name=None, data_dir=DATA_DIR, return_flags=False):
    run_dir = Path(run_dir)
    metadata = json.loads((run_dir / "metadata.json").read_text())
    name = model_name or metadata["selected_model"]
    if name not in metadata["models"]:
        raise ValueError(f"Unknown model: {name}")
    if frame.empty:
        raise ValueError("Input CSV is empty")
    # joblib files must only be loaded from trusted local runs.
    bundle = joblib.load(run_dir / "preprocessing.joblib")
    builder = bundle["builder"]
    builder.data_dir = data_dir
    x = builder.transform(frame)
    x = np.clip(bundle["preprocessing"].transform(x), -10, 10).astype(np.float32)
    if name in {"mlp", "radial"}:
        checkpoint = torch.load(
            run_dir / f"{name}.pt", weights_only=True, map_location="cpu"
        )
        spec = checkpoint["spec"]
        model = (
            KineticsMLP(spec["input_dim"])
            if name == "mlp"
            else PhysicsAdd3DGraphNet(spec["input_dim"])
        )
        model.load_state_dict(checkpoint["state_dict"])
        values = neural_predict(model, x, builder.last_geometry_)
        values = values * spec["target_scale"] + spec["target_mean"]
    else:
        values = model_predict(joblib.load(run_dir / f"{name}.joblib"), x, frame)
    return (values, builder.last_flags_.copy()) if return_flags else values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    args = parser.parse_args()
    frame = pd.read_csv(args.input)
    values, flags = predict(
        args.run, frame, args.model, args.data_dir, return_flags=True
    )
    result = frame.drop(columns=["log10_value"], errors="ignore").copy()
    result["prediction_log10"] = values
    result["prediction_dataset_units"] = np.power(10.0, values)
    result["invalid_smiles"] = flags[:, 0].astype(bool)
    result["missing_embedding"] = flags[:, 1].astype(bool)
    result["missing_structure"] = flags[:, 2].astype(bool)
    result["prediction_target"] = json.loads((args.run / "metadata.json").read_text())[
        "target"
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Saved {len(result)} predictions: {args.output}")


if __name__ == "__main__":
    main()
