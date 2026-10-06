"""Cell-type specificity scores per gene."""

from __future__ import annotations

import numpy as np
import pandas as pd


def fraction_specificity(expr: pd.DataFrame) -> pd.DataFrame:
    """Share of a gene's expression that falls in each cell type.

    expr is genes x cell types on a linear scale (e.g. mean CPM). Each row sums to 1;
    genes with zero expression everywhere get NaN. This is the specificity used by
    EWCE / MAGMA cell typing (Skene & Grant 2016).
    """
    total = expr.sum(axis=1)
    return expr.div(total.where(total > 0), axis=0)


def tau(expr: pd.DataFrame, log: bool = True) -> pd.Series:
    """Tau tissue-specificity index (0 = uniform, 1 = one cell type) per gene.

    Computed on log1p(expr) by default, as in Yanai et al. 2005 / Kryuchkova-Mostacci &
    Robinson-Rechavi 2017. Genes with zero expression everywhere get NaN.
    """
    x = np.log1p(expr.to_numpy()) if log else expr.to_numpy()
    mx = x.max(axis=1, keepdims=True)
    n = x.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (1 - x / mx).sum(axis=1) / (n - 1)
    t[mx[:, 0] == 0] = np.nan
    return pd.Series(t, index=expr.index, name="tau")
