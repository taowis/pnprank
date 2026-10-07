"""Stage 5: temporal back-test and interval simulation (design in docs/decisions.md).

At each cutoff, models are trained on genes known before it only and every other gene
is scored out of fold (`pnprank.backtest.oof_scores_at_cutoff`). Scores are cached in
results/backtest/scores/ so the interval simulation reuses cutoffs already run.

Usage: python workflow/scripts/stage5_backtest.py [n_jobs]
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score

from joblib import Parallel, delayed

from pnprank.backtest import (
    auroc,
    bootstrap_ci,
    evaluation_labels,
    oof_scores_at_cutoff,
    paired_bootstrap_diff,
    training_labels,
    within_fold_percentile,
)
from pnprank.baselines import within_bin_percentile
from pnprank.features import atlas_profiles, label_features, load_celltype_map

n_jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 1
cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
fc, mc, bc = cfg["features"], cfg["model"], cfg["backtest"]
seed, nboot, ci = cfg["seeds"]["global"], cfg["bootstrap"]["n_resamples"], cfg["bootstrap"]["ci"]
out = Path("results/backtest")
(out / "scores").mkdir(parents=True, exist_ok=True)
root = Path(paths["data_root"])

X_static = pd.read_csv(root / "processed" / "features_static.tsv.gz", sep="\t", index_col=0)
missing = [a for a in fc["atlases"] if not any(c.startswith(a + "|") for c in X_static.columns)]
if missing:
    raise SystemExit(f"features_static.tsv.gz lacks atlases {missing}; rerun stage4_model.py")
universe = list(X_static.index)
profiles = atlas_profiles(fc["atlases"], root / "interim", load_celltype_map(fc["celltype_map"]),
                          cfg["pseudobulk"]["min_cells_per_group"], cfg["pseudobulk"]["min_donors"])
non_atlas = [c for c in X_static.columns if c.startswith(("gnomad|", "gtex|"))]

truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
negatives = set(universe) - set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id)
expr = X_static["gtex|mean_all|logtpm"]
baselines = {"gtex_nerve": X_static["gtex|Nerve_Tibial|logtpm"],
             "gtex_mean": X_static["gtex|mean_all|logtpm"],
             "neg_loeuf": -X_static["gnomad|loeuf"]}
TIERS = {"tier1": {"1"}, "tier1+2": {"1", "2"}}


def label_fn(known):
    return label_features(profiles, universe, known)


def scores_at(set_name: str, cutoff: int) -> pd.DataFrame:
    """Out-of-fold model scores plus baselines for all genes scored at this cutoff."""
    f = out / "scores" / f"{set_name}_{cutoff}.tsv.gz"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0)
    y = training_labels(universe, truth, cutoff, TIERS[set_name])
    t0 = time.time()
    full, p_full, fold = oof_scores_at_cutoff(X_static, y, "logistic", seed, bc["scoring_folds"],
                                              mc["inner_folds"], label_fn, n_jobs)
    non, p_non, fold2 = oof_scores_at_cutoff(X_static[non_atlas], y, "logistic", seed,
                                             bc["scoring_folds"], mc["inner_folds"], None, n_jobs)
    assert fold.equals(fold2)
    df = pd.DataFrame({"fold": fold, "full": full, "non_atlas": non})
    for k, v in baselines.items():
        df[k] = v.reindex(df.index)
    df.to_csv(f, sep="\t", compression="gzip")
    print(f"{set_name} {cutoff}: n_pos_train={int(y.sum())} scored={len(df)} "
          f"C_full={[p['logisticregression__C'] for p in p_full]} "
          f"C_non={[p['logisticregression__C'] for p in p_non]} {time.time() - t0:.0f}s", flush=True)
    return df


def ap(s, y):
    ok = ~np.isnan(s)
    return average_precision_score(y[ok], s[ok])


SCORES = ["full", "non_atlas", "gtex_nerve", "gtex_mean", "neg_loeuf"]

# Score transforms. Primary: percentile within the fold model that scored the gene
# (fold models are tuned separately, so raw scores differ in scale between folds; see
# docs/decisions.md). "pooled" is the pre-specified analysis (raw scores pooled across
# folds) and is kept for transparency.
TRANSFORMS = ("fold_rank", "fold_rank_expr_matched", "pooled", "pooled_expr_matched")


def transformed(df: pd.DataFrame) -> dict:
    e = expr.reindex(df.index)
    out = {}
    for s in SCORES:
        out[("fold_rank", s)] = within_fold_percentile(df[s], df["fold"])
        out[("fold_rank_expr_matched", s)] = within_fold_percentile(df[s], df["fold"], e,
                                                                    bc["expression_bins"])
        out[("pooled", s)] = df[s]
        out[("pooled_expr_matched", s)] = within_bin_percentile(df[s], e, bc["expression_bins"])
    return out


def evaluate(set_name: str, cutoff: int, df: pd.DataFrame, k0: int):
    tr = transformed(df)
    perf, pairs = [], []
    k = k0
    for window in (None, bc["interval_window_years"]):
        ev = evaluation_labels(truth, negatives, cutoff, TIERS[set_name], window)
        ev = ev[ev.index.isin(df.index)]
        lab = ev.to_numpy().astype(bool)
        wname = "all_after" if window is None else f"within_{window}y"
        for t in TRANSFORMS:
            for s in SCORES:
                sc = tr[(t, s)].reindex(ev.index).to_numpy(float)
                a = bootstrap_ci(sc, lab, auroc, nboot, ci, seed + k)
                k += 1
                perf.append({"positives": set_name, "cutoff": cutoff, "eval": wname,
                             "transform": t, "score": s, "n_pos": int(lab.sum()),
                             "n_neg": int((~lab).sum()), "auroc": round(a[0], 3),
                             "auroc_lo": round(a[1], 3), "auroc_hi": round(a[2], 3),
                             "avg_precision": round(ap(sc, lab), 4)})
                if s == "full":
                    continue
                d, lo, hi, p = paired_bootstrap_diff(
                    tr[(t, "full")].reindex(ev.index).to_numpy(float), sc, lab, auroc, nboot, ci,
                    seed)
                pairs.append({"positives": set_name, "cutoff": cutoff, "eval": wname,
                              "transform": t, "a": "full", "b": s, "n_pos": int(lab.sum()),
                              "auroc_diff": round(d, 3), "ci_lo": round(lo, 3),
                              "ci_hi": round(hi, 3), "frac_le0": round(p, 3)})
    return perf, pairs


# --- back-test at fixed cutoffs ---------------------------------------------------
jobs = [(set_name, cutoff) for set_name in TIERS for cutoff in bc["cutoff_years"]]
cached = {j: scores_at(*j) for j in jobs}
res = Parallel(n_jobs=n_jobs)(delayed(evaluate)(sn, c, cached[(sn, c)], 1000 * i)
                              for i, (sn, c) in enumerate(jobs))
pd.DataFrame([r for p, _ in res for r in p]).to_csv(out / "backtest_performance.tsv", sep="\t",
                                                    index=False)
pd.DataFrame([r for _, q in res for r in q]).to_csv(out / "backtest_paired_differences.tsv",
                                                    sep="\t", index=False)
print("back-test evaluation written", flush=True)

# --- interval simulation (tier 1): rank each gene in its discovery year -------------
t1 = truth[truth.tier == "1"].assign(year=lambda d: d.discovery_year.astype(int))
lo_y, hi_y = bc["interval_simulation_years"]
rows = []
for year in sorted(set(t1.year[(t1.year >= lo_y) & (t1.year <= hi_y)])):
    df = scores_at("tier1", year)
    tr = transformed(df)
    for _, r in t1[t1.year == year].iterrows():
        g = r.ensembl_gene_id
        row = {"ensembl_gene_id": g, "hgnc_symbol": r.hgnc_symbol, "year": year,
               "category": r.disease_category, "n_scored": len(df),
               "fold_size": int((df["fold"] == df.at[g, "fold"]).sum()) if g in df.index else np.nan}
        for s in SCORES:
            row[f"pct_{s}"] = round(float(tr[("fold_rank", s)].get(g, np.nan)), 4)
        row["pct_full_expr_matched"] = round(float(tr[("fold_rank_expr_matched", "full")].get(g, np.nan)), 4)
        row["pct_full_pooled"] = round(float(df["full"].rank(pct=True).get(g, np.nan)), 4)
        rows.append(row)
pd.DataFrame(rows).to_csv(out / "interval_simulation.tsv", sep="\t", index=False)
print("done")
