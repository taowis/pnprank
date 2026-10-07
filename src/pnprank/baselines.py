"""Simple baselines the model must beat, and controls for expression-level confounding."""

from __future__ import annotations

import pandas as pd


def within_bin_percentile(score: pd.Series, covariate: pd.Series, n_bins: int = 10) -> pd.Series:
    """Rank `score` within quantile bins of `covariate` (e.g. expression level).

    Returns each gene's percentile (0-1] among genes in the same covariate bin, so a
    downstream AUROC reflects the score after matching on the covariate. Genes with a
    missing score or covariate get NaN.
    """
    df = pd.DataFrame({"s": score, "c": covariate}).dropna()
    df["bin"] = pd.qcut(df["c"].rank(method="first"), n_bins, labels=False)
    pct = df.groupby("bin")["s"].rank(pct=True)
    return pct.reindex(score.index)
