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


def test_auroc_and_bootstrap():
    from pnprank.backtest import auroc, bootstrap_ci

    s = np.array([0.9, 0.8, 0.1, 0.2, np.nan])
    y = np.array([1, 1, 0, 0, 1], dtype=bool)
    assert auroc(s, y) == 1.0
    assert auroc(np.array([1, 1, 1, 1.0]), np.array([1, 0, 1, 0], dtype=bool)) == 0.5
    est, lo, hi = bootstrap_ci(s, y, n_resamples=50, seed=1)
    assert est == 1.0 and lo <= est <= hi


def test_within_bin_percentile_removes_covariate():
    from pnprank.backtest import auroc
    from pnprank.baselines import within_bin_percentile

    rng = np.random.default_rng(0)
    cov = pd.Series(rng.normal(size=2000))
    labels = (cov + rng.normal(scale=0.5, size=2000)) > 1  # labels driven by covariate only
    score = cov + rng.normal(scale=0.1, size=2000)  # score is just the covariate
    assert auroc(score.to_numpy(), labels.to_numpy()) > 0.85
    matched = within_bin_percentile(score, cov, n_bins=20)
    assert abs(auroc(matched.to_numpy(), labels.to_numpy()) - 0.5) < 0.1


def test_paired_bootstrap_diff():
    from pnprank.backtest import paired_bootstrap_diff

    rng = np.random.default_rng(0)
    y = np.r_[np.ones(50), np.zeros(500)].astype(bool)
    good = y + rng.normal(scale=0.5, size=len(y))
    noise = rng.normal(size=len(y))
    d, lo, hi, p = paired_bootstrap_diff(good, noise, y, n_resamples=200, seed=1)
    assert d > 0.3 and lo > 0 and p == 0
    d, lo, hi, p = paired_bootstrap_diff(good, good, y, n_resamples=50, seed=1)
    assert d == 0 and lo == 0 and hi == 0
