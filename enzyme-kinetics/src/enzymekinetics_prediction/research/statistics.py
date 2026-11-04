"""Replicate-level tests and correlation diagnostics; frames are not replicates."""

import itertools
import numpy as np


def paired_sign_flip(differences):
    d = np.asarray(differences, float)
    if d.ndim != 1 or not 2 <= len(d) <= 20 or not np.isfinite(d).all():
        raise ValueError("Require 2..20 finite independent paired differences")
    observed = abs(d.mean())
    samples = np.array(
        [
            abs(np.mean(d * np.array(sign)))
            for sign in itertools.product([-1, 1], repeat=len(d))
        ]
    )
    return {
        "n_replicates": len(d),
        "mean_difference": float(d.mean()),
        "p_two_sided_exact": float(np.mean(samples >= observed - 1e-14)),
        "assumption": "Independent replicate differences exchangeable under sign reversal; not frames/atoms",
    }


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, float)
    if (
        p.ndim != 1
        or not len(p)
        or not np.isfinite(p).all()
        or np.any((p < 0) | (p > 1))
    ):
        raise ValueError("Invalid p-values")
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def autocorrelation_diagnostics(values):
    x = np.asarray(values, float)
    if x.ndim != 1 or len(x) < 4 or not np.isfinite(x).all():
        raise ValueError("Need >=4 finite sequential samples")
    x = x - x.mean()
    variance = x @ x
    if variance < 1e-25:
        return {
            "tau_integrated_frames": None,
            "effective_samples": None,
            "status": "constant_or_numerically_unresolved",
        }
    # Initial-positive paired sequence, conservative truncation for short paths.
    fft = np.fft.rfft(x, n=2 * len(x))
    acf = np.fft.irfft(fft * fft.conjugate())[: len(x)] / variance
    total = 0.0
    for lag in range(1, len(x) - 1, 2):
        pair = acf[lag] + acf[lag + 1]
        if pair <= 0:
            break
        total += pair
    tau = max(1.0, 1 + 2 * total)
    return {
        "tau_integrated_frames": float(tau),
        "effective_samples": float(len(x) / tau),
        "status": "exploratory_short_series_estimate",
    }
