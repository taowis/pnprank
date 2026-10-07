"""Stage 5: repeat the Stage 4 current-label comparison (logistic full vs non-atlas,
tier 1) over several CV seeds to measure split-to-split variance.

Usage: python workflow/scripts/stage5_repeated_cv.py [n_jobs]
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml
from joblib import Parallel, delayed

from pnprank.backtest import auroc
from pnprank.features import atlas_profiles, label_features, load_celltype_map
from pnprank.model import cross_validated_scores

n_jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 1
cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
fc, mc = cfg["features"], cfg["model"]
root = Path(paths["data_root"])
out = Path("results/backtest")
out.mkdir(parents=True, exist_ok=True)

X_static = pd.read_csv(root / "processed" / "features_static.tsv.gz", sep="\t", index_col=0)
universe = list(X_static.index)
profiles = atlas_profiles(fc["atlases"], root / "interim", load_celltype_map(fc["celltype_map"]),
                          cfg["pseudobulk"]["min_cells_per_group"], cfg["pseudobulk"]["min_donors"])
non_atlas = [c for c in X_static.columns if c.startswith(("gnomad|", "gtex|"))]
truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
pos = set(truth.ensembl_gene_id[truth.tier == "1"])
drop = set(cands.ensembl_gene_id) - pos
genes = [g for g in universe if g not in drop]
y = pd.Series([int(g in pos) for g in genes], index=genes)


def label_fn(known):
    return label_features(profiles, universe, known)


def one(seed: int) -> dict:
    full, _ = cross_validated_scores(X_static.loc[genes], y, "logistic", seed, mc["cv_folds"],
                                     mc["inner_folds"], label_fn)
    non, _ = cross_validated_scores(X_static.loc[genes, non_atlas], y, "logistic", seed,
                                    mc["cv_folds"], mc["inner_folds"], None)
    a, b = auroc(full.to_numpy(), y.to_numpy()), auroc(non.to_numpy(), y.to_numpy())
    return {"cv_seed": seed, "auroc_full": round(a, 4), "auroc_non_atlas": round(b, 4),
            "diff": round(a - b, 4)}


base = cfg["seeds"]["global"]
seeds = [base + i for i in range(cfg["backtest"]["repeated_cv_seeds"])]
rows = Parallel(n_jobs=n_jobs)(delayed(one)(s) for s in seeds)
df = pd.DataFrame(rows)
df.to_csv(out / "repeated_cv_tier1.tsv", sep="\t", index=False)
print(df.to_string())
print(df[["auroc_full", "auroc_non_atlas", "diff"]].describe().round(4).to_string())
