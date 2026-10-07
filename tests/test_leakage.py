"""No gene discovered at or after the cutoff may influence any feature or label.

Must pass before any back-test result is reported.
"""

import numpy as np
import pandas as pd
import pytest

from pnprank.coexpression import known_gene_similarity
from pnprank.features import known_before, label_features
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


@pytest.mark.skip(reason="activated in Stage 5 once the back-test pipeline exists")
def test_backtest_pipeline_end_to_end():
    pass
