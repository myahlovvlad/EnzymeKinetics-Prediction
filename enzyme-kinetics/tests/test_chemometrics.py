"""Projection and fallback properties independent of enzyme labels."""

import numpy as np
import pytest
from sklearn.linear_model import Ridge

from enzymekinetics_prediction.learning.chemometrics import (
    ConsensusMBPLS,
    OPLSRegressor,
    ECLocalRegressor,
)
from enzymekinetics_prediction.learning.training import grouped_split, split_groups


def signal():
    x = np.random.default_rng(42).normal(size=(100, 8))
    y = x[:, 0] * 2 - x[:, 3]
    return x, y


def test_mbpls_invariant_to_uniform_block_rescaling():
    x, y = signal()
    blocks = {"sequence": (0, 3), "molecule": (3, 8)}
    a = ConsensusMBPLS(blocks, 4).fit(x, y)
    scaled = x.copy()
    scaled[:, :3] *= 100
    b = ConsensusMBPLS(blocks, 4).fit(scaled, y)
    np.testing.assert_allclose(a.predict(x), b.predict(scaled), atol=1e-9)
    assert set(a.block_coefficients()) == set(blocks)
    with pytest.raises(ValueError, match="partition"):
        ConsensusMBPLS({"a": (0, 5), "b": (4, 8)}).fit(x, y)


def test_opls_removes_only_response_orthogonal_variation():
    x, y = signal()
    model = OPLSRegressor(2).fit(x, y)
    standardized = (x - model.mean_) / model.scale_
    filtered = model.transform(x)
    np.testing.assert_allclose(
        filtered.T @ (y - y.mean()), standardized.T @ (y - y.mean()), atol=1e-9
    )
    assert np.mean(np.abs(model.predict(x) - y)) < 0.1


def test_ec_global_fallback_for_unseen_and_small_groups():
    x, y = signal()
    ec = ["1.1.1.1"] * 90 + ["2.2.2.2"] * 10
    model = ECLocalRegressor(Ridge(), min_group_size=20).fit(x, y, ec)
    assert set(model.local_) == {"1.1"}
    np.testing.assert_allclose(
        model.predict(x, ["9.9"] * 100), model.global_.predict(x)
    )
    np.testing.assert_allclose(
        model.predict(x[90:], ec[90:]), model.global_.predict(x[90:])
    )


def test_ec_grouping_respects_requested_level():
    import pandas as pd

    frame = pd.DataFrame(
        {"ec": [f"{i}.1.{j}.1" for i in range(1, 9) for j in range(2)]}
    )
    train, val = grouped_split(frame, group_by="ec", ec_level=2)
    groups = split_groups(frame, "ec", 2)
    assert set(groups.iloc[train]).isdisjoint(groups.iloc[val])
