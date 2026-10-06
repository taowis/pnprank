import numpy as np
import pandas as pd
import scipy.sparse as sp

from pnprank.pseudobulk import cpm, mean_over_donors, pseudobulk
from pnprank.specificity import fraction_specificity, tau


def test_pseudobulk_sums_and_min_cells():
    # 2 genes x 5 cells
    counts = sp.csc_matrix(np.array([[1, 2, 0, 5, 1], [0, 1, 3, 0, 0]]))
    ct = pd.Series(["A", "A", "B", "B", "C"])
    donor = pd.Series(["d1", "d1", "d1", "d1", "d1"])
    sums, n = pseudobulk(counts, ["g1", "g2"], ct, donor, min_cells=2)
    assert list(sums.columns.get_level_values("cell_type")) == ["A", "B"]  # C dropped
    assert sums.loc["g1", ("A", "d1")] == 3 and sums.loc["g2", ("B", "d1")] == 3
    assert n["n_cells"].tolist() == [2, 2]


def test_cpm_and_mean_over_donors():
    cols = pd.MultiIndex.from_tuples(
        [("A", "d1"), ("A", "d2"), ("B", "d1")], names=["cell_type", "donor"]
    )
    sums = pd.DataFrame([[1, 3, 0], [1, 1, 2]], index=["g1", "g2"], columns=cols)
    c = cpm(sums)
    assert np.allclose(c.sum(axis=0), 1e6)
    m = mean_over_donors(c, min_donors=2)
    assert list(m.columns) == ["A"]
    assert np.isclose(m.loc["g1", "A"], (0.5e6 + 0.75e6) / 2)


def test_specificity_scores():
    expr = pd.DataFrame(
        [[10, 0, 0], [5, 5, 5], [0, 0, 0]], index=["spec", "flat", "off"], columns=["A", "B", "C"]
    )
    f = fraction_specificity(expr)
    assert f.loc["spec", "A"] == 1 and np.isclose(f.loc["flat", "B"], 1 / 3)
    assert f.loc["off"].isna().all()
    t = tau(expr)
    assert np.isclose(t["spec"], 1) and np.isclose(t["flat"], 0) and np.isnan(t["off"])
