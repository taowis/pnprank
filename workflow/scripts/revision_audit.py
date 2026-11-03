"""Summary tables for manuscript revision 3, derived from committed results only.

No model is refitted. Inputs are the cached back-test scores in results/backtest/scores/,
the committed result tables and the truth set. Outputs go to results/tables/revision_v3/.

1. eval_sets_tier1.tsv: per cutoff, evaluation positives, background genes with a score
   for each method, prevalence, average precision and AP / prevalence (primary analysis:
   tier 1, all genes discovered at or after the cutoff, within-fold percentiles).
2. sensitivity_eval_only.tsv: evaluation-only sensitivity analyses on the same scores
   (training is unchanged): (a) background widened to every scored gene that is not an
   evaluation positive, i.e. candidate-list genes and later tier-2 genes kept as
   negatives; (b) the two uncertain tier-1 genes (KIF1A, SCN9A) removed from the
   evaluation positives.
3. uncertain_truth_genes.tsv: the low-confidence and manual re-check genes, their tier,
   reason and the role they play at each cutoff.
4. locus_summary.tsv: linkage case-study relative ranks with and without DHMN1.
5. top50_fisher.tsv: the post hoc Fisher test recomputed from the committed table.

Usage: python workflow/scripts/revision_audit.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact
from sklearn.metrics import average_precision_score

from pnprank.backtest import auroc, bootstrap_ci, evaluation_labels, within_fold_percentile

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
bc = cfg["backtest"]
seed, nboot, ci = cfg["seeds"]["global"], cfg["bootstrap"]["n_resamples"], cfg["bootstrap"]["ci"]
out = Path("results/tables/revision_v3")
out.mkdir(parents=True, exist_ok=True)

truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
SCORES = ["full", "non_atlas", "gtex_nerve", "gtex_mean", "neg_loeuf"]


def load(cutoff: int) -> pd.DataFrame:
    return pd.read_csv(f"results/backtest/scores/tier1_{cutoff}.tsv.gz", sep="\t", index_col=0)


def ap(s, y):
    ok = ~np.isnan(s)
    return float(average_precision_score(y[ok], s[ok]))


# --- 1. evaluation sets ------------------------------------------------------------
# The universe is every gene with features; scored genes at a cutoff are the universe
# minus training positives and dropped tier-2 genes. Background genes ("negatives") are
# universe genes not on any candidate list and not in the truth set.
universe = set(pd.read_csv("results/tables/pnprank_ranking.tsv.gz", sep="\t").ensembl_gene_id)
negatives = universe - set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id)
assert len(universe) == 19229 and len(negatives) == 18982, (len(universe), len(negatives))

rows, sens = [], []
k = 0
for c in bc["cutoff_years"]:
    df = load(c)
    expr = df["gtex_mean"]
    ev = evaluation_labels(truth, negatives, c, {"1"})
    ev = ev[ev.index.isin(df.index)]
    lab = ev.to_numpy().astype(bool)
    row = {"cutoff": c, "n_pos": int(lab.sum()), "n_background": int((~lab).sum())}
    for s in SCORES:
        pct = within_fold_percentile(df[s], df["fold"]).reindex(ev.index).to_numpy(float)
        ok = ~np.isnan(pct)
        row[f"n_pos_scored_{s}"] = int((lab & ok).sum())
        row[f"n_background_scored_{s}"] = int((~lab & ok).sum())
        if s in ("full", "non_atlas"):
            row[f"ap_{s}"] = round(ap(pct, lab), 4)
    pm = within_fold_percentile(df["full"], df["fold"], expr, bc["expression_bins"])
    pm = pm.reindex(ev.index).to_numpy(float)
    row["n_background_expr_matched"] = int((~lab & ~np.isnan(pm)).sum())
    row["n_pos_expr_matched"] = int((lab & ~np.isnan(pm)).sum())
    row["prevalence"] = round(lab.mean(), 5)
    row["ap_random_baseline"] = row["prevalence"]
    row["ap_full_over_prevalence"] = round(row["ap_full"] / lab.mean(), 2)
    row["ap_non_atlas_over_prevalence"] = round(row["ap_non_atlas"] / lab.mean(), 2)
    rows.append(row)

    # --- 2. evaluation-only sensitivity analyses ----------------------------------
    pct_full = within_fold_percentile(df["full"], df["fold"])
    pct_non = within_fold_percentile(df["non_atlas"], df["fold"])
    pos = set(ev.index[ev == 1])
    variants = {
        "primary": ev,
        "background_all_scored": pd.Series([int(g in pos) for g in df.index], index=df.index),
        "drop_uncertain_tier1": ev[~ev.index.isin(
            truth.ensembl_gene_id[truth.hgnc_symbol.isin(["KIF1A", "SCN9A"])])],
    }
    for name, e in variants.items():
        y = e.to_numpy().astype(bool)
        for s, p in (("full", pct_full), ("non_atlas", pct_non)):
            est, lo, hi = bootstrap_ci(p.reindex(e.index).to_numpy(float), y, auroc, nboot, ci,
                                       seed + k)
            k += 1
            sens.append({"cutoff": c, "analysis": name, "score": s, "n_pos": int(y.sum()),
                         "n_background": int((~y).sum()), "auroc": round(est, 3),
                         "auroc_lo": round(lo, 3), "auroc_hi": round(hi, 3)})

pd.DataFrame(rows).to_csv(out / "eval_sets_tier1.tsv", sep="\t", index=False)
pd.DataFrame(sens).to_csv(out / "sensitivity_eval_only.tsv", sep="\t", index=False)

# --- 3. uncertain truth genes ------------------------------------------------------
recheck = {"ATM", "FXN", "KIF1A", "MAG", "PEX7", "PNPLA6", "SCARB2", "SLC52A2", "UCHL1", "XPA"}
low = {"ADA2", "AGXT", "ERCC6", "GLA", "HMBS", "MT-TL1", "MTTP", "PPOX", "SAMD9L", "SCN9A",
       "SPG21", "TTC19", "VAPB", "VPS13A", "XPA"}
why_recheck = {g: "discovery paper read from PMC PDFs fetched by an automated downloader"
               for g in ("ATM", "FXN", "KIF1A", "MAG", "XPA")}
why_recheck.update({g: "full-text quote came through a summarising web fetch; not verified "
                       "verbatim" for g in ("PEX7", "PNPLA6", "SCARB2", "SLC52A2", "UCHL1")})
urows = []
for _, r in truth[truth.hgnc_symbol.isin(recheck | low)].iterrows():
    g, yr, tier = r.hgnc_symbol, int(r.discovery_year), r.tier
    reasons = []
    if g in recheck:
        reasons.append("manual re-check: " + why_recheck[g])
    if g in low:
        reasons.append("low-confidence label")
    roles = []
    for c in bc["cutoff_years"]:
        if tier == "1":
            roles.append(f"{c}:{'train+' if yr < c else 'test+'}")
        else:
            roles.append(f"{c}:{'dropped' if yr < c else 'label0/tier1+2 test+'}")
    urows.append({"hgnc_symbol": g, "tier": tier, "category": r.disease_category,
                  "discovery_year": yr, "pmid": r.pmid, "reason": "; ".join(reasons),
                  "role_by_cutoff_tier1_analysis": " ".join(roles)})
pd.DataFrame(urows).sort_values(["tier", "hgnc_symbol"]).to_csv(
    out / "uncertain_truth_genes.tsv", sep="\t", index=False)

# --- 4. linkage case study with and without DHMN1 -------------------------------------
loc = pd.read_csv("results/case_studies/locus_ranking.tsv", sep="\t")
lrows = []
for label, d in (("all_7", loc), ("without_DHMN1", loc[loc.locus_id != "DHMN1"])):
    for s in ("full", "non_atlas", "gtex_nerve"):
        rel = d[f"rank_{s}"] / d["n_scored"]
        lrows.append({"set": label, "score": s, "n_loci": len(d),
                      "median_relative_rank": round(float(rel.median()), 3),
                      "n_top_half": int((rel <= 0.5).sum())})
pd.DataFrame(lrows).to_csv(out / "locus_summary.tsv", sep="\t", index=False)

# --- 5. top-50 Fisher test -------------------------------------------------------------
t50 = pd.read_csv("results/tables/top50_literature.tsv", sep="\t")
e = t50["class"] == "E"
low_non = t50["percentile_non_atlas"] < 0.9
tab = [[int((e & low_non).sum()), int((e & ~low_non).sum())],
       [int((~e & low_non).sum()), int((~e & ~low_non).sum())]]
odds, p = fisher_exact(tab)
pd.DataFrame([{"E_nonatlas_lt0.9": tab[0][0], "E_nonatlas_ge0.9": tab[0][1],
               "other_nonatlas_lt0.9": tab[1][0], "other_nonatlas_ge0.9": tab[1][1],
               "odds_ratio": round(odds, 2), "p_two_sided": round(p, 3)}]).to_csv(
    out / "top50_fisher.tsv", sep="\t", index=False)
print("written to", out)
