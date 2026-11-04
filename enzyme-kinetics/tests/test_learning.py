"""Tests for leakage boundaries, geometry invariance, fitting and artifact reload."""

from argparse import Namespace
import json

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from enzymekinetics_prediction.learning.features import FeatureBuilder, radial_quantiles  # noqa: E402
from enzymekinetics_prediction.learning.predict import predict  # noqa: E402
from enzymekinetics_prediction.learning.training import grouped_split, run, train_neural  # noqa: E402


def sample_frame(n=60):
    return pd.DataFrame(
        {
            "sequence": ["ACDEFGHIK" + "A" * (i % 10) for i in range(n)],
            "uniprot": [f"P{i % 10}" for i in range(n)],
            "ec": ["3.1.1.1"] * n,
            "reactant_smiles": ["CCO" if i % 2 else "CC(=O)O" for i in range(n)],
            "temperature": [25.0] * n,
            "ph": [7.0] * n,
            "log10_value": [float(i % 10) / 5 for i in range(n)],
        }
    )


def test_split_excludes_shared_protein_sequences():
    frame = sample_frame()
    train, val = grouped_split(frame)
    assert set(train).isdisjoint(val)
    assert len(train) + len(val) == len(frame)
    assert set(frame.iloc[train].sequence).isdisjoint(frame.iloc[val].sequence)
    assert np.array_equal(train, grouped_split(frame)[0])


def test_features_do_not_use_labels_or_product_and_report_missingness(tmp_path):
    frame = sample_frame(10)
    builder = FeatureBuilder(tmp_path).fit(frame)
    x = builder.transform(frame)
    changed = frame.assign(log10_value=123456.0, value=987.0, product_smiles="N")
    np.testing.assert_array_equal(x, builder.transform(changed))
    assert builder.last_flags_[:, 2].all()  # missing structure explicitly recorded
    invalid = frame.copy()
    invalid.loc[0, "reactant_smiles"] = "invalid"
    builder.transform(invalid)
    assert builder.last_flags_[0, 0] == 1
    unknown = frame.assign(ec="9.9.9.9")
    assert builder.transform(unknown).shape == x.shape


def test_radial_geometry_invariant_to_rigid_motion_and_node_order():
    coordinates = np.random.default_rng(42).normal(size=(50, 3))
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    np.testing.assert_allclose(
        radial_quantiles(coordinates),
        radial_quantiles(coordinates[::-1] @ rotation + 23.0),
        atol=1e-6,
    )


def test_neural_regressor_can_learn_signal():
    torch.set_num_threads(1)
    rng = np.random.default_rng(42)
    x = rng.normal(size=(100, 4)).astype(np.float32)
    g = np.zeros((100, 32), dtype=np.float32)
    y = (2 * x[:, 0] - x[:, 1]).astype(np.float32)
    _, spec = train_neural(
        "mlp",
        x[:80],
        g[:80],
        y[:80],
        x[80:],
        g[80:],
        y[80:],
        seed=42,
        epochs=60,
        patience=20,
        batch_size=20,
    )
    assert min(h["validation_mae_log10"] for h in spec["history"]) < 0.4


def test_training_checkpoint_reloads_identical_predictions(tmp_path):
    folder = tmp_path / "data" / "Baseline" / "kcat"
    folder.mkdir(parents=True)
    train = sample_frame()
    test = sample_frame(12)
    test["sequence"] += "WW"
    train.to_csv(folder / "kcat-random_trainval.csv", index=False)
    test.to_csv(folder / "kcat-random_test.csv", index=False)
    output = tmp_path / "run"
    args = Namespace(
        data_dir=tmp_path / "data",
        output=output,
        target="kcat",
        models=["mean", "ridge", "mlp", "radial", "mbpls", "opls", "ec_mbpls"],
        seed=42,
        validation_fraction=0.2,
        threads=1,
        epochs=3,
        patience=2,
        batch_size=20,
    )
    run(args)
    saved = pd.read_csv(output / "predictions.csv")
    for name in args.models:
        reloaded = predict(
            output, test.drop(columns="log10_value"), name, tmp_path / "data"
        )
        np.testing.assert_allclose(reloaded, saved[name], atol=1e-6)
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["rows"]["test"] == len(test)
    assert metadata["selected_model"] in args.models


def test_embedding_filename_normalization_prefers_alphafold(tmp_path):
    from enzymekinetics_prediction.learning.training import read_embeddings

    path = tmp_path / "embeddings.pt"
    torch.save(
        {
            "ESMFold-P12345-v1.pdb": torch.zeros(3),
            "AF-P12345-F1-model_v4.pdb": torch.ones(3),
        },
        path,
    )
    embeddings, info = read_embeddings(path)
    np.testing.assert_array_equal(embeddings["P12345"], np.ones(3))
    assert info["duplicate_uniprot_representations"] == 1


def test_pair_overlap_detects_equivalent_smiles():
    from enzymekinetics_prediction.learning.features import enzyme_substrate_pairs

    a = sample_frame(1).assign(reactant_smiles="CCO")
    b = sample_frame(1).assign(reactant_smiles="OCC")
    assert enzyme_substrate_pairs(a) == enzyme_substrate_pairs(b)
