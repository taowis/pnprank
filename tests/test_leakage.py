"""No gene discovered at or after the cutoff may influence any feature or label.

Must pass before any back-test result is reported.
"""

import numpy as np
import pandas as pd

from pnprank.coexpression import known_gene_similarity
from pnprank.features import known_before, label_features
from pnprank.backtest import evaluation_labels, oof_scores_at_cutoff, training_labels
from pnprank.model import cross_validated_scores
from pnprank.truth import known_before as known_before_rows


def _truth():
    return pd.DataFrame(
        {
            "ensembl_gene_id": ["G0", "G1", "G2", "G3", "G4"],
            "discovery_year": ["1995", "2005", "2009", "2010", "2018"],
            "tier": ["1", "1", "2", "1", "1"],
        }
    )


def _profiles(seed=0, n_genes=60, n_types=6):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    return {
        "atlasA": pd.DataFrame(rng.gamma(2, 50, (n_genes, n_types)), index=genes,
                               columns=[f"t{j}" for j in range(n_types)]),
        "atlasB": pd.DataFrame(rng.gamma(2, 50, (n_genes, n_types)), index=genes,
                               columns=[f"u{j}" for j in range(n_types)]),
    }


def test_known_before_is_strict():
    rows = [
        {"ensembl_gene_id": "ENSG_A", "discovery_year": "2009"},
        {"ensembl_gene_id": "ENSG_B", "discovery_year": "2010"},
        {"ensembl_gene_id": "ENSG_C", "discovery_year": "2011"},
    ]
    assert known_before_rows(rows, 2010) == {"ENSG_A"}
    assert known_before(_truth(), 2010, {"1"}) == {"G0", "G1"}
    assert known_before(_truth(), 2010, {"1", "2"}) == {"G0", "G1", "G2"}


def test_label_features_ignore_post_cutoff_genes():
    """Changing, removing or adding post-cutoff truth must not change any feature."""
    profiles = _profiles()
    universe = list(profiles["atlasA"].index)
    cutoff = 2010
    t = _truth()
    base = label_features(profiles, universe, known_before(t, cutoff, {"1"}))

    t_removed = t[t.discovery_year.astype(int) < cutoff]
    t_moved = t.assign(discovery_year=lambda d: d.discovery_year.where(
        d.discovery_year.astype(int) < cutoff, "2030"))
    t_added = pd.concat([t, pd.DataFrame({"ensembl_gene_id": ["G30", "G31"],
                                          "discovery_year": ["2012", "2020"], "tier": ["1", "1"]})])
    for variant in (t_removed, t_moved, t_added):
        other = label_features(profiles, universe, known_before(variant, cutoff, {"1"}))
        pd.testing.assert_frame_equal(base, other)


def test_label_features_do_change_with_pre_cutoff_genes():
    """Sanity: the test above would catch leakage, because known genes do matter."""
    profiles = _profiles()
    universe = list(profiles["atlasA"].index)
    a = label_features(profiles, universe, {"G0", "G1"})
    b = label_features(profiles, universe, {"G0", "G1", "G3"})
    assert not a.equals(b)


def test_coexpression_leaves_self_out():
    prof = np.log1p(_profiles()["atlasA"])
    s = known_gene_similarity(prof, {"G0"})
    # G0's only known gene is itself, so after leaving self out it has no comparison
    assert np.isnan(s["G0"])
    assert not np.isnan(s["G1"])


def test_cv_label_features_use_training_fold_only():
    """Every call to the label-feature function inside CV gets only training positives."""
    profiles = _profiles(n_genes=80)
    universe = list(profiles["atlasA"].index)
    rng = np.random.default_rng(1)
    X_static = pd.DataFrame(rng.normal(size=(80, 3)), index=universe, columns=["a", "b", "c"])
    y = pd.Series(0, index=universe)
    y.iloc[:20] = 1
    seen = []

    def label_fn(known):
        seen.append(set(known))
        return label_features(profiles, universe, known)

    oof, _ = cross_validated_scores(X_static, y, "logistic", seed=0, outer_folds=4,
                                    inner_folds=3, label_fn=label_fn)
    positives = set(y.index[y == 1])
    assert len(seen) == 4
    for known in seen:
        assert known < positives and len(known) == 15  # 3/4 of 20 positives per fold
    assert not oof.isna().any()


def _backtest_scores(truth, cutoff, tiers=frozenset({"1"})):
    profiles = _profiles(n_genes=80)
    universe = list(profiles["atlasA"].index)
    rng = np.random.default_rng(2)
    X_static = pd.DataFrame(rng.normal(size=(80, 3)), index=universe, columns=["a", "b", "c"])
    y = training_labels(universe, truth, cutoff, set(tiers))
    scores, _ = oof_scores_at_cutoff(
        X_static, y, "logistic", seed=0, folds=3, inner_folds=2,
        label_fn=lambda known: label_features(profiles, universe, known))
    return y, scores


def _truth_backtest():
    years = [1995, 1998, 2001, 2003, 2004, 2006, 2007, 2009, 2012, 2015, 2018, 2021]
    return pd.DataFrame({"ensembl_gene_id": [f"G{i}" for i in range(12)],
                         "discovery_year": [str(y) for y in years],
                         "tier": ["1"] * 10 + ["2", "2"]})


def test_backtest_pipeline_end_to_end():
    """Back-test scores at a cutoff must not change when post-cutoff truth is removed,
    moved later, or extended; and positives (pre-cutoff genes) are never scored."""
    cutoff = 2008
    t = _truth_backtest()
    y, base = _backtest_scores(t, cutoff)
    assert set(y.index[y == 1]) == {"G0", "G1", "G2", "G3", "G4", "G5", "G6"}
    assert not set(base.index) & set(y.index[y == 1])
    assert not base.isna().any()
    pre = t.discovery_year.astype(int) < cutoff
    variants = [
        t[pre],
        t.assign(discovery_year=t.discovery_year.where(pre, "2040")),
        pd.concat([t, pd.DataFrame({"ensembl_gene_id": ["G40", "G41"],
                                    "discovery_year": ["2010", "2019"], "tier": ["1", "2"]})]),
        t.assign(tier=t.tier.where(pre, "2")),
    ]
    for v in variants:
        yv, sv = _backtest_scores(v, cutoff)
        pd.testing.assert_series_equal(y, yv)
        pd.testing.assert_series_equal(base, sv)


def test_backtest_drops_pre_cutoff_other_tier_only():
    t = _truth_backtest().assign(tier=["1"] * 9 + ["2", "2", "2"])
    t.loc[9, "discovery_year"] = "2000"  # tier 2, known before the cutoff
    y = training_labels([f"G{i}" for i in range(20)], t, 2008, {"1"})
    assert "G9" not in y.index           # known (tier 2) before cutoff: dropped
    assert y["G10"] == 0 and y["G11"] == 0  # tier 2 after cutoff: unknown then, label 0
    assert y["G8"] == 0                  # tier 1 discovered 2012: label 0 at 2008


def test_evaluation_labels_window():
    t = _truth_backtest()
    neg = {f"N{i}" for i in range(5)}
    e = evaluation_labels(t, neg, 2008, {"1"})
    assert set(e.index[e == 1]) == {"G7", "G8", "G9"}
    w = evaluation_labels(t, neg, 2008, {"1"}, window=5)
    assert set(w.index[w == 1]) == {"G7", "G8"}
