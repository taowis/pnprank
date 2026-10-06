"""Co-expression with known genes; must use only pre-cutoff genes.

This is the only label-derived feature family. Every function here takes the set of
known genes explicitly, so the caller decides what "known" means at a cutoff, and
`tests/test_leakage.py` checks that post-cutoff genes cannot change the result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def standardise_profiles(logexpr: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene's profile across cell types (rows = genes).

    Genes with zero variance across cell types get NaN (no shape to correlate).
    """
    x = logexpr.to_numpy(dtype=float)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (x - mu) / sd
    z[sd[:, 0] == 0] = np.nan
    return pd.DataFrame(z, index=logexpr.index, columns=logexpr.columns)


def known_gene_similarity(logexpr: pd.DataFrame, known: set[str]) -> pd.Series:
    """Mean Pearson correlation of each gene's cell-type profile with known genes'.

    A known gene's own profile is excluded from its mean (leave-self-out), so a
    gene's label never feeds its own feature. Requires >= 3 cell types.
    """
    if logexpr.shape[1] < 3:
        raise ValueError("need at least 3 cell types for a correlation profile")
    z = standardise_profiles(logexpr)
    k = [g for g in z.index if g in known and not z.loc[g].isna().any()]
    if not k:
        return pd.Series(np.nan, index=z.index, name="coexp_known")
    zz = np.nan_to_num(z.to_numpy())
    zk = z.loc[k].to_numpy()
    n = z.shape[1]
    corr = zz @ zk.T / n  # genes x known; Pearson r because rows are z-scored
    total = corr.sum(axis=1)
    count = np.full(len(z), len(k), dtype=float)
    pos = {g: i for i, g in enumerate(k)}
    for gi, g in enumerate(z.index):
        if g in pos:  # leave self out
            total[gi] -= corr[gi, pos[g]]
            count[gi] -= 1
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / count
    out[z.isna().any(axis=1).to_numpy()] = np.nan
    return pd.Series(out, index=z.index, name="coexp_known")
