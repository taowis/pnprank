"""Figure and headline table for the comparator back-test (comparators_backtest.py).

Usage: python workflow/scripts/comparators_summarise.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from pnprank.plotting import auroc_by_cutoff

out = Path("results/backtest/comparators")
perf = pd.read_csv(out / "comparator_performance.tsv", sep="\t")
pairs = pd.read_csv(out / "comparator_paired_differences.tsv", sep="\t")

# --- figure: AUROC by cutoff, tier 1 and tier 1+2, raw and expression-matched ---------
SHOWN = ["full", "non_atlas", "gtex_nerve", "go_similarity", "ppi_rwr", "mis_z"]
d = perf[(perf["eval"] == "all_after") & perf.score.isin(SHOWN)].copy()
d["adjust"] = d["positives"] + "|" + d["transform"]
panels = [("tier1|fold_rank", "Tier 1: AUROC"),
          ("tier1|fold_rank_expr_matched", "Tier 1: expression-matched"),
          ("tier1+2|fold_rank", "Tier 1+2: AUROC"),
          ("tier1+2|fold_rank_expr_matched", "Tier 1+2: expression-matched")]
auroc_by_cutoff(d, panels, str(out / "comparators_auroc.png"))

# --- headline table: PNPrank minus each comparator ---------------------------------
p = pairs[(pairs["eval"] == "all_after")].copy()
p["cell"] = (p.auroc_diff.map("{:+.3f}".format) + " [" + p.ci_lo.map("{:+.3f}".format) + ", "
             + p.ci_hi.map("{:+.3f}".format) + "]")
p["excludes_zero"] = (p.ci_lo > 0) | (p.ci_hi < 0)
p[["positives", "transform", "cutoff", "b", "n_pos", "cell", "excludes_zero"]].rename(
    columns={"b": "comparator", "cell": "auroc_diff_ci"}).to_csv(
    out / "comparator_headline_differences.tsv", sep="\t", index=False)
print(p.groupby(["positives", "transform", "b"]).excludes_zero.sum().unstack().to_string())
