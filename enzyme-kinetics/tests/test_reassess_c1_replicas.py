"""Registered replica-level statistics and replay guards, without backend runs."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from enzymekinetics_prediction.research.statistics import autocorrelation_diagnostics

spec = importlib.util.spec_from_file_location(
    "reassess_c1_replicas",
    Path(__file__).resolve().parents[1] / "scripts/diagnostics/reassess_c1_replicas.py",
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_exact_replica_comparison_and_two_family_bh():
    cat, pos = module.COMPARISONS
    results = [{cat: -1.0, pos: value} for value in (-1.0, 1.0, 1.0)]
    tests = module.compare_replicas(results)
    assert tests[cat]["n_replicates"] == 3
    assert tests[cat]["p_two_sided_exact"] == 0.25
    assert tests[pos]["p_two_sided_exact"] == 1.0
    assert tests[cat]["q_bh"] == 0.5
    assert tests[pos]["q_bh"] == 1.0
    assert tests[cat]["mean_difference"] == -1.0


def test_reported_acf_matches_independent_alternating_sequence_reference():
    values = [1, -1, 1, -1]
    original = autocorrelation_diagnostics(values)
    extended = autocorrelation_diagnostics(values, include_acf=True)
    np.testing.assert_allclose(extended.pop("acf"), [1, -0.75, 0.5, -0.25], atol=1e-15)
    assert extended == original
    assert extended["effective_samples"] == 4
    constant = autocorrelation_diagnostics([2] * 4, include_acf=True)
    assert constant["acf"] is None and constant["effective_samples"] is None


def test_replay_requires_inputs_and_never_overwrites_a_run(tmp_path):
    sources = [
        tmp_path / name for name in ("replica-101", "replica-202", "replica-303")
    ]
    output = tmp_path / "new-results"
    with pytest.raises(FileNotFoundError, match="BLOCKED"):
        module.reassess(sources, output)
    assert not output.exists()
    output.mkdir()
    sentinel = output / "historical.txt"
    sentinel.write_text("preserve")
    with pytest.raises(FileExistsError):
        module.reassess(sources, output)
    assert sentinel.read_text() == "preserve"
