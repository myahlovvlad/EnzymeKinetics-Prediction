"""Explicit single-y block-scaled PLS, OPLS, and EC-local/global fallback models."""

import copy

import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.cross_decomposition import PLSRegression


def arrays(x, y=None):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all():
        raise ValueError("Expected a finite feature matrix")
    if y is None:
        return x
    y = np.asarray(y, dtype=float).reshape(-1)
    if len(y) != len(x) or len(y) < 2 or not np.isfinite(y).all():
        raise ValueError("Invalid scalar response")
    return x, y


class ConsensusMBPLS(RegressorMixin, BaseEstimator):
    """PLS on trace-normalized semantic blocks, an explicit consensus formulation.

    Each block is centered and scaled to unit training Frobenius norm per row;
    standard PLS is then fitted to the concatenation. Block loadings/contributions
    are exposed. This does not assert superiority over equally scaled PLS.
    """

    def __init__(self, blocks=None, n_components=8):
        self.blocks = blocks
        self.n_components = n_components

    def fit(self, x, y):
        x, y = arrays(x, y)
        self.n_features_in_ = x.shape[1]
        blocks = self.blocks or {"all": (0, x.shape[1])}
        covered = []
        self.block_slices_ = {}
        self.mean_ = x.mean(axis=0)
        self.scales_ = np.ones(x.shape[1])
        for name, (start, stop) in blocks.items():
            if not 0 <= start < stop <= x.shape[1]:
                raise ValueError("Invalid block slice")
            covered.extend(range(start, stop))
            self.block_slices_[name] = (start, stop)
            trace = np.mean(
                np.sum((x[:, start:stop] - self.mean_[start:stop]) ** 2, axis=1)
            )
            self.scales_[start:stop] = np.sqrt(trace) if trace > 1e-15 else 1.0
        if sorted(covered) != list(range(x.shape[1])):
            raise ValueError("Blocks must partition features exactly without overlap")
        transformed = (x - self.mean_) / self.scales_
        if self.n_components < 1:
            raise ValueError("n_components must be positive")
        n = min(self.n_components, x.shape[1], len(x) - 1)
        self.model_ = PLSRegression(n_components=n, scale=False).fit(transformed, y)
        self.block_loadings_ = {
            name: self.model_.x_loadings_[a:b].copy()
            for name, (a, b) in self.block_slices_.items()
        }
        return self

    def predict(self, x):
        x = arrays(x)
        return self.model_.predict((x - self.mean_) / self.scales_).reshape(-1)

    def block_coefficients(self):
        coefficients = np.asarray(self.model_.coef_).reshape(-1) / self.scales_
        return {
            name: coefficients[a:b].tolist()
            for name, (a, b) in self.block_slices_.items()
        }


class OPLSRegressor(RegressorMixin, BaseEstimator):
    """Single-y OPLS: remove y-orthogonal X variation, then fit one predictive LV."""

    def __init__(self, n_orthogonal=1):
        self.n_orthogonal = n_orthogonal

    def fit(self, x, y):
        x, y = arrays(x, y)
        if self.n_orthogonal < 0:
            raise ValueError("n_orthogonal cannot be negative")
        self.n_features_in_ = x.shape[1]
        self.mean_ = x.mean(axis=0)
        self.scale_ = x.std(axis=0)
        self.scale_[self.scale_ < 1e-12] = 1
        filtered = (x - self.mean_) / self.scale_
        response = y - y.mean()
        self.orthogonal_ = []
        for _ in range(min(self.n_orthogonal, len(x) - 2, x.shape[1] - 1)):
            w = filtered.T @ response
            norm = np.linalg.norm(w)
            if norm < 1e-12:
                break
            w /= norm
            t = filtered @ w
            p = filtered.T @ t / (t @ t)
            wo = p - w * (w @ p)
            if np.linalg.norm(wo) < 1e-12:
                break
            wo /= np.linalg.norm(wo)
            to = filtered @ wo
            if to @ to < 1e-12:
                break
            po = filtered.T @ to / (to @ to)
            self.orthogonal_.append((wo.copy(), po.copy()))
            filtered -= np.outer(to, po)
        self.model_ = PLSRegression(n_components=1, scale=False).fit(filtered, y)
        return self

    def transform(self, x):
        filtered = (arrays(x) - self.mean_) / self.scale_
        for w, p in self.orthogonal_:
            filtered -= np.outer(filtered @ w, p)
        return filtered

    def predict(self, x):
        return self.model_.predict(self.transform(x)).reshape(-1)


class ECLocalRegressor(RegressorMixin, BaseEstimator):
    """Global model plus separately fitted projections for sufficiently large EC groups."""

    def __init__(self, estimator=None, min_group_size=100, ec_level=2):
        self.estimator = estimator
        self.min_group_size = min_group_size
        self.ec_level = ec_level

    def groups(self, ec):
        return np.asarray([".".join(str(v).split(".")[: self.ec_level]) for v in ec])

    def fit(self, x, y, ec):
        x, y = arrays(x, y)
        if len(ec) != len(x) or self.min_group_size < 3 or self.ec_level < 1:
            raise ValueError("Invalid EC grouping configuration")
        prototype = self.estimator if self.estimator is not None else ConsensusMBPLS()
        self.global_ = clone(prototype).fit(x, y)
        self.local_ = {}
        groups = self.groups(ec)
        self.group_sizes_ = {}
        for group in np.unique(groups):
            mask = groups == group
            self.group_sizes_[group] = int(mask.sum())
            if mask.sum() >= self.min_group_size and np.std(y[mask]) > 1e-12:
                self.local_[group] = copy.deepcopy(
                    clone(prototype).fit(x[mask], y[mask])
                )
        return self

    def predict(self, x, ec):
        x = arrays(x)
        if len(ec) != len(x):
            raise ValueError("EC label count mismatch")
        predicted = self.global_.predict(x)
        groups = self.groups(ec)
        for group, model in self.local_.items():
            mask = groups == group
            if mask.any():
                predicted[mask] = model.predict(x[mask])
        return predicted
