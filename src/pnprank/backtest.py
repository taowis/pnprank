"""Temporal split, interval simulation and evaluation metrics."""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Area under the ROC curve via the Mann-Whitney U statistic (ties averaged).

    NaN scores are dropped. Returns NaN if either class is empty.
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    ok = ~np.isnan(scores)
    s, y = scores[ok], labels[ok]
    n_pos, n_neg = y.sum(), (~y).sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def bootstrap_ci(
    scores: np.ndarray,
    labels: np.ndarray,
    stat=auroc,
    n_resamples: int = 1000,
    ci: float = 0.95,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Point estimate and percentile bootstrap CI, resampling positives and negatives
    separately (stratified) so every resample has both classes."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(labels), np.flatnonzero(~labels)
    est = stat(scores, labels)
    boots = np.empty(n_resamples)
    for i in range(n_resamples):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        boots[i] = stat(scores[idx], labels[idx])
    lo, hi = np.nanpercentile(boots, [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100])
    return est, float(lo), float(hi)


def paired_bootstrap_diff(
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    labels: np.ndarray,
    stat=auroc,
    n_resamples: int = 1000,
    ci: float = 0.95,
    seed: int = 0,
) -> tuple[float, float, float, float]:
    """stat(a) - stat(b) on the same stratified resamples.

    Returns (difference, ci_low, ci_high, fraction of resamples with difference <= 0).
    """
    a = np.asarray(scores_a, dtype=float)
    b = np.asarray(scores_b, dtype=float)
    y = np.asarray(labels, dtype=bool)
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y), np.flatnonzero(~y)
    d0 = stat(a, y) - stat(b, y)
    diffs = np.empty(n_resamples)
    for i in range(n_resamples):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        diffs[i] = stat(a[idx], y[idx]) - stat(b[idx], y[idx])
    lo, hi = np.nanpercentile(diffs, [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100])
    return d0, float(lo), float(hi), float(np.mean(diffs <= 0))
