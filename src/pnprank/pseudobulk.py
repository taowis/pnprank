"""Aggregate single-cell counts to cell-type pseudobulk profiles."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def pseudobulk(
    counts: sp.spmatrix,
    genes: list[str],
    cell_type: pd.Series,
    donor: pd.Series,
    min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts per (cell type, donor).

    `counts` is genes x cells (raw integer counts), aligned with `cell_type` and `donor`.
    Groups with fewer than `min_cells` cells are dropped. Returns (sums, n_cells): sums is
    genes x groups with a MultiIndex (cell_type, donor) on columns.
    """
    counts = sp.csc_matrix(counts)
    if counts.shape[1] != len(cell_type) or len(cell_type) != len(donor):
        raise ValueError("counts columns, cell_type and donor must align")
    groups = pd.MultiIndex.from_arrays(
        [cell_type.to_numpy(), donor.to_numpy()], names=["cell_type", "donor"]
    )
    codes, uniques = pd.factorize(groups)
    sizes = np.bincount(codes, minlength=len(uniques))
    keep = np.flatnonzero(sizes >= min_cells)
    # indicator matrix cells x groups
    ind = sp.csr_matrix(
        (np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniques))
    )
    sums = (counts @ ind).toarray()[:, keep]
    cols = pd.MultiIndex.from_tuples([uniques[k] for k in keep], names=["cell_type", "donor"])
    n_cells = pd.DataFrame({"n_cells": sizes[keep]}, index=cols)
    return pd.DataFrame(sums, index=genes, columns=cols), n_cells


def cpm(sums: pd.DataFrame) -> pd.DataFrame:
    """Counts per million within each pseudobulk sample."""
    lib = sums.sum(axis=0)
    if (lib == 0).any():
        raise ValueError("empty pseudobulk sample")
    return sums / lib * 1e6


def mean_over_donors(expr: pd.DataFrame, min_donors: int = 1) -> pd.DataFrame:
    """Average per-donor profiles within each cell type (genes x cell types).

    Cell types seen in fewer than `min_donors` donors are dropped.
    """
    n = expr.columns.get_level_values("cell_type").value_counts()
    keep = n[n >= min_donors].index
    out = expr.T.groupby(level="cell_type").mean().T
    return out[[c for c in out.columns if c in keep]]
