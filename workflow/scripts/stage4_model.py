"""Stage 4: build the full feature matrix and evaluate models by nested CV on current labels.

This is NOT the temporal back-test (Stage 5). It answers: with all currently known genes,
do atlas features add signal beyond simple non-atlas baselines?

Usage: python workflow/scripts/stage4_model.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score

from pnprank.backtest import auroc, bootstrap_ci, paired_bootstrap_diff
from pnprank.features import (
    atlas_features,
    atlas_profiles,
    gnomad_features,
    gtex_features,
    label_features,
    load_celltype_map,
)
from pnprank.model import cross_validated_scores, fit_tuned, with_label_features

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
fc, mc = cfg["features"], cfg["model"]
seed, nboot = cfg["seeds"]["global"], cfg["bootstrap"]["n_resamples"]
out = Path("results/tables/stage4")
out.mkdir(parents=True, exist_ok=True)
interim = Path(paths["data_root"]) / "interim"

# --- universe and labels ----------------------------------------------------------
hgnc = pd.read_csv(Path(cfg["truth"]["sources_dir"]) / "hgnc_complete_set.txt", sep="\t", dtype=str)
pc = hgnc[(hgnc.status == "Approved") & (hgnc.locus_group == "protein-coding gene")]
pc = pc.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id")
symbol = dict(zip(pc.ensembl_gene_id, pc.symbol))
truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)

# --- features ---------------------------------------------------------------------
cmap = load_celltype_map(fc["celltype_map"])
profiles = atlas_profiles(fc["atlases"], interim, cmap, cfg["pseudobulk"]["min_cells_per_group"],
                          cfg["pseudobulk"]["min_donors"])
universe_all = sorted(pc.ensembl_gene_id)
X_atlas = atlas_features(profiles, universe_all)
X_gnomad = gnomad_features(cfg["external"]["gnomad_constraint"]["dest"], universe_all)
X_gtex = gtex_features(cfg["external"]["gtex_median_tpm"]["dest"], universe_all, fc["gtex_tissues"])
# keep genes seen in at least one atlas or in GTEx (drops genes with no expression data at all)
seen = X_atlas.notna().any(axis=1) | X_gtex.notna().any(axis=1)
universe_all = [g for g in universe_all if seen[g]]
X_static = pd.concat([X_atlas, X_gnomad, X_gtex], axis=1).loc[universe_all]
X_static.to_csv(Path(paths["data_root"]) / "processed" / "features_static.tsv.gz", sep="\t",
                compression="gzip")


def label_fn(known: set[str]) -> pd.DataFrame:
    return label_features(profiles, universe_all, known)


def labelled(tiers: set[str]) -> pd.Series:
    """Positives: truth genes in `tiers`. Negatives: protein-coding genes that are not
    candidates at all. Other candidates (other tiers, unverified, excluded) are dropped."""
    pos = set(truth.ensembl_gene_id[truth.tier.isin(tiers)])
    drop = set(cands.ensembl_gene_id) - pos
    genes = [g for g in universe_all if g not in drop]
    return pd.Series([int(g in pos) for g in genes], index=genes)


def evaluate(scores: pd.Series, y: pd.Series, name: str, k: int) -> dict:
    s, lab = scores.reindex(y.index).to_numpy(float), y.to_numpy().astype(bool)
    a = bootstrap_ci(s, lab, auroc, nboot, cfg["bootstrap"]["ci"], seed + k)

    def ap(sc, yy):
        ok = ~np.isnan(sc)
        return average_precision_score(yy[ok], sc[ok])

    p = bootstrap_ci(s, lab, ap, nboot, cfg["bootstrap"]["ci"], seed + 1000 + k)
    return {"score": name, "n_pos": int(lab.sum()), "n_neg": int((~lab).sum()),
            "prevalence": round(lab.mean(), 4),
            "auroc": round(a[0], 3), "auroc_lo": round(a[1], 3), "auroc_hi": round(a[2], 3),
            "avg_precision": round(p[0], 4), "ap_lo": round(p[1], 4), "ap_hi": round(p[2], 4)}


non_atlas = [c for c in X_static.columns if c.startswith(("gnomad|", "gtex|"))]
results, params_log = [], {}
k = 0
for set_name, tiers in (("tier1", {"1"}), ("tier1+2", {"1", "2"})):
    y = labelled(tiers)
    # single-feature baselines (no fitting)
    single = {
        "baseline:gtex_mean_all_logtpm": X_static["gtex|mean_all|logtpm"],
        "baseline:gtex_nerve_tibial_logtpm": X_static["gtex|Nerve_Tibial|logtpm"],
        "baseline:gnomad_neg_loeuf": -X_static["gnomad|loeuf"],
    }
    for name, sc in single.items():
        results.append({"positives": set_name, **evaluate(sc, y, name, k)})
        k += 1
    # fitted models by nested CV
    specs = [
        ("logistic:non_atlas", X_static[non_atlas], None, "logistic"),
        ("logistic:full", X_static, label_fn, "logistic"),
        ("hgb:non_atlas", X_static[non_atlas], None, "hgb"),
        ("hgb:full", X_static, label_fn, "hgb"),
    ]
    for name, X, lf, kind in specs:
        oof, params = cross_validated_scores(X.loc[y.index], y, kind, seed, mc["cv_folds"],
                                             mc["inner_folds"], lf)
        oof.to_csv(out / f"oof_{set_name}_{name.replace(':', '_')}.tsv.gz", sep="\t", header=["score"])
        params_log[f"{set_name}|{name}"] = params
        results.append({"positives": set_name, **evaluate(oof, y, name, k)})
        k += 1
        print(results[-1], flush=True)

# paired comparisons on identical resamples: does adding atlas features help?
pairs = []
for set_name, tiers in (("tier1", {"1"}), ("tier1+2", {"1", "2"})):
    y = labelled(tiers)
    oofs = {n: pd.read_csv(out / f"oof_{set_name}_{n.replace(':', '_')}.tsv.gz", sep="\t",
                           index_col=0)["score"].reindex(y.index)
            for n in ("logistic:full", "logistic:non_atlas", "hgb:full", "hgb:non_atlas")}
    for a, b in (("logistic:full", "logistic:non_atlas"), ("hgb:full", "hgb:non_atlas"),
                 ("logistic:full", "hgb:full")):
        d, lo, hi, p = paired_bootstrap_diff(oofs[a].to_numpy(), oofs[b].to_numpy(),
                                             y.to_numpy().astype(bool), auroc, nboot,
                                             cfg["bootstrap"]["ci"], seed)
        pairs.append({"positives": set_name, "a": a, "b": b, "auroc_diff": round(d, 3),
                      "ci_lo": round(lo, 3), "ci_hi": round(hi, 3), "frac_resamples_le0": round(p, 3)})
pd.DataFrame(pairs).to_csv(out / "paired_auroc_differences.tsv", sep="\t", index=False)
print(pd.DataFrame(pairs).to_string())

res = pd.DataFrame(results)
res.to_csv(out / "cv_performance.tsv", sep="\t", index=False)
json.dump(params_log, open(out / "cv_best_params.json", "w"), indent=1)

# --- final model on all current tier-1 labels: ranking and coefficients ---------------
y = labelled({"1"})
known = set(y.index[y == 1])
X_full = with_label_features(X_static, label_fn, known)
est, best = fit_tuned(X_full.loc[y.index], y.to_numpy(), "logistic", seed, mc["inner_folds"])
lr = est.named_steps["logisticregression"]
imp = est.named_steps["simpleimputer"]
names = list(X_full.columns) + [f"missing:{X_full.columns[i]}" for i in imp.indicator_.features_]
coef = pd.Series(lr.coef_[0], index=names).sort_values(key=abs, ascending=False)
coef.to_csv(out / "logistic_full_coefficients.tsv", sep="\t", header=["coef_standardised"])
score = pd.Series(est.predict_proba(X_full.loc[universe_all])[:, 1], index=universe_all)
status = pd.Series("negative", index=universe_all)
for t in ("1", "2"):
    status[status.index.isin(truth.ensembl_gene_id[truth.tier == t])] = f"truth_tier{t}"
status[status.index.isin(set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id))] = "candidate_not_verified"
rank = pd.DataFrame({"hgnc_symbol": [symbol[g] for g in universe_all], "score": score,
                     "status": status}).sort_values("score", ascending=False)
rank.insert(0, "rank", range(1, len(rank) + 1))
rank.index.name = "ensembl_gene_id"
rank.to_csv(out / "ranking_current_logistic_full.tsv.gz", sep="\t", compression="gzip")
json.dump({"atlases": fc["atlases"], "n_universe": len(universe_all),
           "n_features_static": X_static.shape[1], "n_features_full": X_full.shape[1],
           "final_logistic_params": best,
           "cell_types": {a: list(p.columns) for a, p in profiles.items()}},
          open(out / "summary.json", "w"), indent=1)
print(res.to_string())
