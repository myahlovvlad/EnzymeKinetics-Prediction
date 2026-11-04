"""Scientific boundaries of the small H1 experiment, with no network or weights."""

import json
import sys
import copy
from types import SimpleNamespace

import joblib
import numpy as np
import pandas as pd
import pytest
import torch

from enzymekinetics_prediction.learning.benchmark_h1 import (
    cluster_intervals,
    interval_radius,
    prepare,
    scores,
    paired_ablation,
)
from enzymekinetics_prediction.learning.embeddings import (
    embed_sequences,
    embedding_identity,
    cached_sources,
    validate_embedding_cache,
    SOURCE_FILES,
)
from enzymekinetics_prediction.learning.predict import predict
from enzymekinetics_prediction.learning import training


def frame():
    return pd.DataFrame(
        {
            "sequence": ["A" * i + "G" for i in range(1, 61)],
            "uniprot": [f"P{i}" for i in range(60)],
            "ec": [f"1.{i % 3}.1.1" for i in range(60)],
            "reactant_smiles": ["CCO"] * 60,
            "sequence_40cluster": [str(i) for i in range(60)],
            "ph": [7.0] * 60,
            "temperature": [25.0] * 60,
            "log10_value": np.linspace(-1, 2, 60),
            "source_record": np.arange(1, 61),
        }
    )


def test_partition_is_label_independent_and_cluster_disjoint():
    config = {
        "group_by": "sequence_40cluster",
        "max_length": 256,
        "max_rows": 60,
        "split_seed": 42,
    }
    original = frame()
    changed = original.copy()
    changed.log10_value = changed.log10_value.iloc[::-1].to_numpy()
    a, b = prepare(original, config), prepare(changed, config)
    for name in a:
        assert a[name].source_record.tolist() == b[name].source_record.tolist()
    names = list(a)
    for i, name in enumerate(names):
        for other in names[i + 1 :]:
            assert not set(a[name].sequence_40cluster) & set(
                a[other].sequence_40cluster
            )


def test_independent_metric_reference_and_bootstrap():
    y = np.arange(12, dtype=float)
    p = y + 1
    result = scores(y, p)
    assert result["mae"] == result["rmse"] == 1
    assert result["spearman"] == pytest.approx(1)
    ci = cluster_intervals(y, p, np.repeat(np.arange(6), 2), replicates=40)
    assert ci == cluster_intervals(y, p, np.repeat(np.arange(6), 2), replicates=40)
    assert ci["mae"]["ci95"] == [1.0, 1.0]
    assert scores(y, np.ones(12))["spearman"] is None


def test_conformal_rank_and_cluster_maxima():
    assert interval_radius(np.arange(1, 10), np.zeros(9), np.arange(9)) == 9
    assert interval_radius(np.arange(8), np.zeros(8), np.arange(8)) is None
    # Two rows in one cluster contribute one maximum, not two independent scores.
    assert (
        interval_radius(np.array([1, 5, 2, 4]), np.zeros(4), [0, 0, 1, 1], alpha=0.5)
        == 5
    )


def test_paired_ablation_uses_same_cluster_draws():
    y = np.arange(12, dtype=float)
    errors = np.arange(1, 13, dtype=float)
    groups = np.repeat(np.arange(3), [1, 3, 8])
    result = paired_ablation(y, y + errors, y + errors + 1, groups, 40)
    assert result["delta_mae"] == 1.0
    assert result["ci95"] == [1.0, 1.0]
    assert paired_ablation(y, y + errors, y + errors, groups, 40)["ci95"] == [0.0, 0.0]


def test_development_training_never_opens_official_test(tmp_path, monkeypatch):
    root = tmp_path / "Baseline/kcat"
    root.mkdir(parents=True)
    frame().to_csv(root / "kcat-random_trainval.csv", index=False)
    (root / "kcat-random_test.csv").write_text("SEALED: must never be parsed")
    original = training.load_table

    def guarded(path):
        assert path.name != "kcat-random_test.csv"
        return original(path)

    monkeypatch.setattr(training, "load_table", guarded)
    splits = []
    for seed in (11, 33):
        args = SimpleNamespace(
            target="kcat",
            data_dir=tmp_path,
            output=tmp_path / str(seed),
            seed=seed,
            split_seed=42,
            development_only=True,
            split_scheme="provided-test",
            group_by="sequence_40cluster",
            ec_level=2,
            validation_fraction=0.25,
            epochs=1,
            patience=1,
            batch_size=8,
            threads=1,
            models=["ridge"],
            embedding_path=tmp_path / "none.pt",
            physical_table=[],
        )
        training.run(args)
        meta = json.loads((args.output / "metadata.json").read_text())
        assert meta["official_test_accessed"] is False
        assert set(meta["rows"]) == {"train", "validation"}
        split = pd.read_csv(args.output / "split.csv")
        splits.append(split)
        bundle = joblib.load(args.output / "preprocessing.joblib")
        assert len(bundle["builder"].encoder_.categories_[0]) <= 3
        values = predict(
            args.output, frame().drop(columns="log10_value"), data_dir=tmp_path
        )
        assert values.shape == (60,) and np.isfinite(values).all()
    pd.testing.assert_frame_equal(*splits)


def test_offline_encoder_requires_pin_before_loading(tmp_path):
    with pytest.raises(ValueError, match="immutable revision"):
        embed_sequences(["ACDE"], tmp_path / "x.pt", local_files_only=True)


def test_legacy_embeddings_remain_readable_but_not_a_strict_cache(tmp_path):
    path = tmp_path / "legacy.pt"
    torch.save({"AC": torch.ones(2)}, path)
    path.with_suffix(".pt.meta.json").write_text('{"provider":"legacy"}')
    vectors, info = training.read_embeddings(path)
    np.testing.assert_array_equal(vectors["AC"], np.ones(2))
    assert info["status"] == "loaded"
    expected = embedding_identity(
        ["AC"],
        "facebook/esm2_t6_8M_UR50D",
        "a" * 40,
        {name: "b" * 64 for name in SOURCE_FILES},
    )
    with pytest.raises(ValueError, match="legacy"):
        training.read_embeddings(path, expected, 2)


def test_strict_cache_rejects_changed_inputs_and_payload(tmp_path, monkeypatch):
    """A toy frozen encoder checks the same boundary, without external weights."""
    sources = tmp_path / "sources"
    sources.mkdir()
    for name in SOURCE_FILES:
        (sources / name).write_text(name)
    revision = "a" * 40
    model_id = "facebook/esm2_t6_8M_UR50D"

    def network_forbidden(*args, **kwargs):
        pytest.fail("Offline provider attempted a network model lookup")

    monkeypatch.setitem(
        sys.modules,
        "huggingface_hub",
        SimpleNamespace(
            model_info=network_forbidden,
            try_to_load_from_cache=lambda model, name, revision: str(sources / name),
        ),
    )

    class Tokenizer:
        @classmethod
        def from_pretrained(cls, name, **kwargs):
            assert name == model_id and kwargs["revision"] == revision
            assert kwargs["local_files_only"] and not kwargs["trust_remote_code"]
            return cls()

        def __call__(self, sequences, **kwargs):
            assert not kwargs["truncation"]
            # [special, residue, residue, special]; special values must not be pooled.
            return {
                "input_ids": torch.tensor([[100, 2, 4, 100]] * len(sequences)),
                "attention_mask": torch.ones((len(sequences), 4), dtype=torch.long),
                "special_tokens_mask": torch.tensor([[1, 0, 0, 1]] * len(sequences)),
            }

    class Encoder(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.scale = torch.nn.Parameter(torch.ones(1))
            self.config = SimpleNamespace(max_position_embeddings=1026)

        @classmethod
        def from_pretrained(cls, name, **kwargs):
            assert kwargs["local_files_only"] and kwargs["use_safetensors"]
            return cls()

        def forward(self, input_ids, attention_mask):
            assert not self.training and not torch.is_grad_enabled()
            assert not self.scale.requires_grad
            return SimpleNamespace(
                last_hidden_state=input_ids[:, :, None].float().repeat(1, 1, 2)
                * self.scale
            )

    monkeypatch.setitem(
        sys.modules,
        "transformers",
        SimpleNamespace(AutoTokenizer=Tokenizer, AutoModel=Encoder),
    )
    path = tmp_path / "vectors.pt"
    meta = embed_sequences(
        ["AC", "DE"], path, model_id, revision, local_files_only=True
    )
    identity = embedding_identity(
        ["DE", "AC", "AC"], model_id, revision, cached_sources(model_id, revision)
    )
    assert validate_embedding_cache(path, identity, 2) == meta["cache_key"]
    assert torch.equal(
        torch.load(path, weights_only=True)["AC"], torch.tensor([3.0, 3.0])
    )
    for key, value in [
        ("revision", "b" * 40),
        ("algorithm", "changed"),
        ("sequence_sha256", "f" * 64),
        (
            "source_sha256",
            {**identity["source_sha256"], "tokenizer_config.json": "f" * 64},
        ),
    ]:
        changed = copy.deepcopy(identity)
        changed[key] = value
        with pytest.raises(ValueError, match="identity"):
            validate_embedding_cache(path, changed, 2)
    (sources / "special_tokens_map.json").write_text("changed tokenizer")
    with pytest.raises(ValueError, match="identity"):
        validate_embedding_cache(
            path,
            embedding_identity(
                ["AC", "DE"], model_id, revision, cached_sources(model_id, revision)
            ),
            2,
        )
    path.write_bytes(path.read_bytes() + b"altered")
    with pytest.raises(ValueError, match="checksum"):
        validate_embedding_cache(path, identity, 2)
    with pytest.raises(ValueError, match="already exists"):
        embed_sequences(["AC"], path, model_id, revision, local_files_only=True)
