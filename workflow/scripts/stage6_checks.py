"""Stage 6: descriptive checks of the released ranking (post hoc, not pre-specified).

- Held-out tier-2 genes (never used in training) vs unlabelled genes: AUROC of each
  percentile column, and the paired difference full vs non-atlas.
- Composition of the top of the list: status and GTEx mean-expression decile.

Usage: python workflow/scripts/stage6_checks.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from pnprank.backtest import auroc, bootstrap_ci, paired_bootstrap_diff
from pnprank.ranking import load_ranking

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
seed, nboot, ci = cfg["seeds"]["global"], cfg["bootstrap"]["n_resamples"], cfg["bootstrap"]["ci"]
out = Path("results/tables/stage6")
out.mkdir(parents=True, exist_ok=True)
t = load_ranking("results/tables/pnprank_ranking.tsv.gz")

d = t[t.status.isin(["known_tier2", "unlabelled"])]
rows = []
for c in ["percentile", "percentile_expr_matched", "percentile_non_atlas"]:
    ok = d[c].notna()
    lab = (d.status[ok] == "known_tier2").to_numpy()
    a, lo, hi = bootstrap_ci(d[c][ok].to_numpy(float), lab, auroc, nboot, ci, seed)
    rows.append({"check": "tier2_held_out_auroc", "score": c, "n_pos": int(lab.sum()),
                 "value": round(a, 3), "ci_lo": round(lo, 3), "ci_hi": round(hi, 3)})
dd = d.dropna(subset=["percentile", "percentile_non_atlas"])
lab = (dd.status == "known_tier2").to_numpy()
diff, lo, hi, _ = paired_bootstrap_diff(dd.percentile.to_numpy(float),
                                        dd.percentile_non_atlas.to_numpy(float), lab, auroc,
                                        nboot, ci, seed)
rows.append({"check": "tier2_held_out_auroc_diff", "score": "percentile - percentile_non_atlas",
             "n_pos": int(lab.sum()), "value": round(diff, 3), "ci_lo": round(lo, 3),
             "ci_hi": round(hi, 3)})
for n in (100, 500):
    top = t[t["rank"] <= n]
    for s, k in top.status.value_counts().items():
        rows.append({"check": f"top{n}_status", "score": s, "value": int(k)})
    rows.append({"check": f"top{n}_frac_expr_decile10", "score": "gtex_mean_decile",
                 "value": round(float((top.gtex_mean_decile == 10).mean()), 3)})
    rows.append({"check": f"top{n}_frac_expr_decile_le5", "score": "gtex_mean_decile",
                 "value": round(float((top.gtex_mean_decile <= 5).mean()), 3)})
res = pd.DataFrame(rows)
res.to_csv(out / "release_checks.tsv", sep="\t", index=False)
print(res.to_string(index=False))
