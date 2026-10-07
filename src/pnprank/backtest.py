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


# --- temporal back-test -----------------------------------------------------------


def training_labels(universe: list[str], truth, cutoff: int, tiers: set[str]):
    """Training labels at a cutoff, using only what was known before it.

    1 = truth genes in `tiers` discovered before `cutoff`. Truth genes of other tiers
    discovered before `cutoff` are dropped (they were known). Every other gene is 0,
    including genes discovered at or after the cutoff and candidates whose status is
    only known today.
    """
    import pandas as pd

    from pnprank.features import known_before

    pos = known_before(truth, cutoff, tiers)
    drop = known_before(truth, cutoff, {"1", "2"}) - pos
    genes = [g for g in universe if g not in drop]
    return pd.Series([int(g in pos) for g in genes], index=genes)


def oof_scores_at_cutoff(X_static, y, kind: str, seed: int, folds: int, inner_folds: int,
                         label_fn=None, n_jobs: int = 1):
    """Score every gene labelled 0 out of fold.

    Label-derived features are built once from the positives in `y` (pre-cutoff genes).
    Genes labelled 0 are split into `folds`; each fold is scored by a model trained on
    all positives and the remaining folds, with hyperparameters tuned by inner CV on that
    training set. Positives are not scored.

    Returns (scores, best params per fold, fold index per gene). Scores from different
    fold models are on different scales (hyperparameters differ), so they must be
    compared within fold (`within_fold_percentile`) before pooling.
    """
    import pandas as pd
    from joblib import Parallel, delayed
    from sklearn.model_selection import KFold

    from pnprank.model import fit_tuned, with_label_features

    known = set(y.index[y == 1])
    X = with_label_features(X_static, label_fn, known).loc[y.index]
    zeros = y.index[y == 0].to_numpy()
    splits = list(KFold(n_splits=folds, shuffle=True, random_state=seed).split(zeros))

    def run(k, tr, te):
        # positives in y's order (never set order: it varies with hash randomisation)
        train = list(y.index[y == 1]) + list(zeros[tr])
        est, best = fit_tuned(X.loc[train], y.loc[train].to_numpy(), kind, seed + k, inner_folds)
        return pd.Series(est.predict_proba(X.loc[zeros[te]])[:, 1], index=zeros[te]), best

    out = Parallel(n_jobs=n_jobs)(delayed(run)(k, tr, te) for k, (tr, te) in enumerate(splits))
    fold = pd.Series(0, index=zeros)
    for k, (_, te) in enumerate(splits):
        fold.iloc[te] = k
    return pd.concat([s for s, _ in out]).reindex(zeros), [b for _, b in out], fold


def within_fold_percentile(scores, fold, covariate=None, n_bins: int = 10):
    """Percentile (0-1] of each score among genes scored by the same fold model.

    With `covariate`, genes are additionally grouped into quantile bins of it (computed
    over all genes), so the percentile is also matched on that covariate (e.g. expression
    level). Missing scores or covariates give NaN.
    """
    import pandas as pd

    df = pd.DataFrame({"s": scores, "g": fold.reindex(scores.index)})
    if covariate is not None:
        c = covariate.reindex(scores.index)
        df = df[c.notna()]
        b = pd.qcut(c[c.notna()].rank(method="first"), n_bins, labels=False)
        df["g"] = df["g"] * n_bins + b
    df = df.dropna()
    return df.groupby("g")["s"].rank(pct=True).reindex(scores.index)


def evaluation_labels(truth, negatives: set[str], cutoff: int, tiers: set[str],
                      window: int | None = None):
    """Evaluation set at a cutoff: truth genes in `tiers` discovered at or after `cutoff`
    (and before cutoff + window, if given) are 1; `negatives` are 0."""
    import pandas as pd

    t = truth[truth.tier.isin(tiers)]
    yr = t.discovery_year.astype(int)
    sel = yr >= cutoff
    if window is not None:
        sel &= yr < cutoff + window
    pos = set(t.ensembl_gene_id[sel])
    genes = sorted(pos | (negatives - pos))
    return pd.Series([int(g in pos) for g in genes], index=genes)


def percentile_rank(scores, genes) -> dict:
    """Percentile (0-1, 1 = top) of each gene in `genes` among all scored genes."""
    r = scores.rank(pct=True)
    return {g: float(r[g]) for g in genes if g in r.index and not np.isnan(r[g])}
