"""Stage 5: summary tables, failure analysis and figures from stage5_backtest.py outputs.

Usage: python workflow/scripts/stage5_summarise.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from pnprank.plotting import auroc_by_cutoff, percentile_strip

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
out, fig = Path("results/backtest"), Path("results/figures")
X = pd.read_csv(Path(paths["data_root"]) / "processed" / "features_static.tsv.gz", sep="\t",
                index_col=0)

perf = pd.read_csv(out / "backtest_performance.tsv", sep="\t")
for set_name in ("tier1", "tier1+2"):
    d = perf[(perf.positives == set_name) & (perf["eval"] == "all_after")]
    auroc_by_cutoff(d, [("raw", f"{set_name}: raw AUROC"),
                        ("expr_matched", f"{set_name}: expression-matched")],
                    str(fig / f"backtest_auroc_{set_name.replace('+', '_')}.png"))

# --- interval simulation ------------------------------------------------------------
sim = pd.read_csv(out / "interval_simulation.tsv", sep="\t")
cols = ["pct_full", "pct_non_atlas", "pct_gtex_nerve", "pct_gtex_mean", "pct_neg_loeuf",
        "pct_full_expr_matched"]
rng = np.random.default_rng(cfg["seeds"]["global"])


def summary(d: pd.DataFrame, c: str) -> dict:
    v = d[c].dropna().to_numpy()
    boots = [np.median(rng.choice(v, len(v))) for _ in range(cfg["bootstrap"]["n_resamples"])]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"n": len(v), "median_pct": round(float(np.median(v)), 3),
            "median_lo": round(float(lo), 3), "median_hi": round(float(hi), 3),
            "top1pct": int((v >= 0.99).sum()), "top5pct": int((v >= 0.95).sum()),
            "top10pct": int((v >= 0.90).sum()), "frac_top10pct": round(float((v >= 0.90).mean()), 3)}


rows = []
for group, d in [("all", sim)] + [(c, sim[sim.category == c]) for c in sorted(set(sim.category))]:
    for c in cols:
        rows.append({"group": group, "score": c.removeprefix("pct_"), **summary(d, c)})
pd.DataFrame(rows).to_csv(out / "interval_summary.tsv", sep="\t", index=False)
percentile_strip(sim, cols[:5], str(fig / "interval_simulation_percentiles.png"))

# --- failure analysis: genes ranked in the bottom half at discovery -------------------
spec_cols = [c for c in X.columns if c.endswith("|spec")]
top_spec = X[spec_cols].idxmax(axis=1).str.replace("|spec", "", regex=False)
fail = sim[sim.pct_full < 0.5].sort_values("pct_full").copy()
fail["gtex_nerve_logtpm"] = X["gtex|Nerve_Tibial|logtpm"].reindex(fail.ensembl_gene_id).round(2).to_numpy()
fail["gtex_mean_logtpm"] = X["gtex|mean_all|logtpm"].reindex(fail.ensembl_gene_id).round(2).to_numpy()
fail["loeuf"] = X["gnomad|loeuf"].reindex(fail.ensembl_gene_id).round(2).to_numpy()
fail["n_atlases_detected"] = X[[c for c in X.columns if c.endswith("|tau")]].notna().sum(
    axis=1).reindex(fail.ensembl_gene_id).to_numpy()
fail["most_specific_type"] = top_spec.reindex(fail.ensembl_gene_id).to_numpy()
fail.to_csv(out / "failure_analysis.tsv", sep="\t", index=False)

# reference: the same feature summaries for genes ranked in the top 10% at discovery
ok = sim[sim.pct_full >= 0.9]
ref = pd.DataFrame({
    "group": ["ranked_top10pct", "ranked_bottom_half", "all_genes"],
    "n": [len(ok), len(fail), np.nan],
    "median_gtex_nerve_logtpm": [X["gtex|Nerve_Tibial|logtpm"].reindex(ok.ensembl_gene_id).median(),
                                 fail.gtex_nerve_logtpm.median(), X["gtex|Nerve_Tibial|logtpm"].median()],
    "median_loeuf": [X["gnomad|loeuf"].reindex(ok.ensembl_gene_id).median(), fail.loeuf.median(),
                     X["gnomad|loeuf"].median()],
}).round(3)
ref.to_csv(out / "failure_feature_contrast.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
print(fail[["hgnc_symbol", "year", "category", "pct_full", "gtex_nerve_logtpm", "loeuf",
            "n_atlases_detected", "most_specific_type"]].to_string())
print(ref.to_string())
