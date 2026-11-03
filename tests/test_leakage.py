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
    # 4 outer folds, each followed by its 3 inner splits
    assert len(seen) == 4 * (1 + 3)
    outer = seen[::4]
    for known in outer:
        assert known < positives and len(known) == 15  # 3/4 of 20 positives per fold
    for k in range(4):
        for inner in seen[4 * k + 1: 4 * k + 4]:
            assert inner < outer[k] and len(inner) == 10  # 2/3 of 15
    assert not oof.isna().any()


def test_inner_cv_label_features_exclude_inner_validation_positives():
    """While hyperparameters are chosen, label features are rebuilt from the inner-training
    positives only, so no inner-validation positive is ever a seed."""
    from sklearn.model_selection import StratifiedKFold

    from pnprank.model import fit_tuned, with_label_features

    profiles = _profiles(n_genes=90)
    universe = list(profiles["atlasA"].index)
    rng = np.random.default_rng(3)
    X_static = pd.DataFrame(rng.normal(size=(90, 3)), index=universe, columns=["a", "b", "c"])
    y = pd.Series(0, index=universe)
    y.iloc[::4] = 1
    seen = []

    def label_fn(known):
        seen.append(set(known))
        return label_features(profiles, universe, known)

    X = with_label_features(X_static, lambda k: label_features(profiles, universe, k),
                            set(y.index[y == 1]))
    fit_tuned(X, y.to_numpy(), "logistic", seed=5, folds=3, label_fn=label_fn)
    splits = StratifiedKFold(n_splits=3, shuffle=True, random_state=5).split(X, y)
    assert len(seen) == 3
    for known, (tr, va) in zip(seen, splits):
        val_pos = set(y.index[va][y.iloc[va] == 1])
        assert val_pos and not (known & val_pos)
        assert known == set(y.index[tr][y.iloc[tr] == 1])


def test_fit_tuned_without_label_features_unchanged():
    """The non-atlas path (no label features) is still a plain grid search."""
    from sklearn.model_selection import GridSearchCV, StratifiedKFold

    from pnprank.model import GRIDS, fit_tuned, make_model

    rng = np.random.default_rng(4)
    X = pd.DataFrame(rng.normal(size=(120, 4)), columns=list("abcd"))
    y = (X.a + rng.normal(scale=2, size=120) > 1).astype(int).to_numpy()
    est, best = fit_tuned(X, y, "logistic", seed=1, folds=3)
    gs = GridSearchCV(make_model("logistic", 1), GRIDS["logistic"], scoring="roc_auc",
                      cv=StratifiedKFold(3, shuffle=True, random_state=1)).fit(X, y)
    assert best == gs.best_params_
    assert np.allclose(est.predict_proba(X), gs.best_estimator_.predict_proba(X))


def _backtest_scores(truth, cutoff, tiers=frozenset({"1"})):
    profiles = _profiles(n_genes=80)
    universe = list(profiles["atlasA"].index)
    rng = np.random.default_rng(2)
    X_static = pd.DataFrame(rng.normal(size=(80, 3)), index=universe, columns=["a", "b", "c"])
    y = training_labels(universe, truth, cutoff, set(tiers))
    scores, _, _ = oof_scores_at_cutoff(
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


def test_within_fold_percentile_ignores_scale_between_folds():
    """Scores on different scales per fold give the same percentiles as unscaled ones."""
    from pnprank.backtest import within_fold_percentile

    idx = [f"G{i}" for i in range(40)]
    rng = np.random.default_rng(3)
    s = pd.Series(rng.normal(size=40), index=idx)
    fold = pd.Series(np.repeat([0, 1], 20), index=idx)
    shifted = s.where(fold == 0, s * 10 - 5)
    pd.testing.assert_series_equal(within_fold_percentile(s, fold), within_fold_percentile(shifted, fold))
    cov = pd.Series(rng.normal(size=40), index=idx)
    m = within_fold_percentile(s, fold, cov, n_bins=2)
    assert m.notna().all() and m.max() == 1.0


def test_backtest_scores_reproducible_across_hash_seeds(tmp_path):
    """Back-test scores must not depend on Python's string-hash randomisation."""
    import subprocess
    import sys

    script = tmp_path / "run.py"
    script.write_text(
        "import numpy as np, pandas as pd\n"
        "from pnprank.backtest import oof_scores_at_cutoff\n"
        "idx=[f'G{i}' for i in range(120)]\n"
        "rng=np.random.default_rng(0)\n"
        "X=pd.DataFrame(rng.normal(size=(120,4)),index=idx)\n"
        "y=pd.Series(0,index=idx); y.iloc[::6]=1\n"
        "s,_,_=oof_scores_at_cutoff(X,y,'logistic',0,3,3)\n"
        "print(','.join(f'{v:.12f}' for v in s))\n")
    outs = []
    for h in ("1", "2", "3"):
        r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                           env={**__import__("os").environ, "PYTHONHASHSEED": h}, check=True)
        outs.append(r.stdout)
    assert outs[0] == outs[1] == outs[2]


# --- comparators (docs/stage_reports/comparators.md) ---------------------------------


def _comparator_inputs():
    import scipy.sparse as sp

    from pnprank.comparators import adjacency

    universe = [f"G{i}" for i in range(40)]
    rng = np.random.default_rng(5)
    go = sp.csr_matrix((rng.random((40, 15)) < 0.25).astype(float))
    edges = pd.DataFrame(rng.choice(universe, size=(80, 2)), columns=["a", "b"])
    edges = edges[edges.a != edges.b]
    return universe, go, adjacency(edges, universe)


def test_comparators_ignore_post_cutoff_genes():
    """GO similarity and network propagation scores at a cutoff must not change when
    post-cutoff truth is removed, moved later, extended or re-tiered; the seed set is
    exactly the pre-cutoff tier-1 genes; seeds are never scored."""
    from pnprank.comparators import go_similarity_scores, rwr_scores, seed_scores_at_cutoff

    universe, go, adj = _comparator_inputs()
    cutoff = 2008
    t = _truth_backtest()
    pre = t.discovery_year.astype(int) < cutoff
    variants = [
        t[pre],
        t.assign(discovery_year=t.discovery_year.where(pre, "2040")),
        pd.concat([t, pd.DataFrame({"ensembl_gene_id": ["G30", "G31"],
                                    "discovery_year": ["2010", "2019"], "tier": ["1", "2"]})]),
        t.assign(tier=t.tier.where(pre, "2")),
    ]
    seen = []

    def go_fn(known):
        seen.append(frozenset(known))
        return go_similarity_scores(go, universe, known)

    def ppi_fn(known):
        seen.append(frozenset(known))
        return rwr_scores(adj, universe, known, restart=0.5)

    for fn in (go_fn, ppi_fn):
        y = training_labels(universe, t, cutoff, {"1"})
        base = seed_scores_at_cutoff(fn, y)
        assert not set(base.index) & set(y.index[y == 1])
        assert base.notna().all()
        for v in variants:
            yv = training_labels(universe, v, cutoff, {"1"})
            pd.testing.assert_series_equal(base, seed_scores_at_cutoff(fn, yv))
    assert set(seen) == {frozenset({"G0", "G1", "G2", "G3", "G4", "G5", "G6"})}


def test_comparators_do_change_with_pre_cutoff_genes():
    """Sanity: the test above would catch leakage, because the seeds do matter."""
    from pnprank.comparators import go_similarity_scores, rwr_scores

    universe, go, adj = _comparator_inputs()
    a, b = {"G0", "G1"}, {"G0", "G1", "G7"}
    assert not go_similarity_scores(go, universe, a).equals(go_similarity_scores(go, universe, b))
    assert not rwr_scores(adj, universe, a).equals(rwr_scores(adj, universe, b))


def test_comparator_annotations_filtered_to_pre_cutoff(tmp_path):
    """The GO annotation reader drops annotations dated on or after the cutoff."""
    import gzip

    from pnprank.comparators import read_gaf

    row = ["UniProt", "P1", "A", "", "GO:1", "r", "IDA", "", "P", "", "", "protein",
           "taxon:9606"]
    p = tmp_path / "a.gaf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(row + ["20071231", "x"]) + "\n")
        fh.write("\t".join(row[:1] + ["P2"] + row[2:] + ["20080101", "x"]) + "\n")
    assert read_gaf(p, before="20080101").db_id.tolist() == ["P1"]


def test_profiles_from_static_round_trip():
    """Profiles rebuilt from the static logcpm columns give the same label features."""
    from pnprank.features import atlas_features, profiles_from_static

    profiles = _profiles(n_genes=50)
    universe = list(profiles["atlasA"].index) + ["absent"]
    X = atlas_features(profiles, universe)
    rebuilt = profiles_from_static(X, ["atlasA", "atlasB"])
    known = {"G1", "G5", "G9", "G20"}
    a = label_features(profiles, universe, known)
    b = label_features(rebuilt, universe, known)
    assert "absent" not in rebuilt["atlasA"].index
    pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-12, atol=1e-12)
