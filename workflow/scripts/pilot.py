"""Stage 3 feasibility pilot on one atlas.

Builds cell-type pseudobulk and specificity for the pilot atlas and asks whether known
neuropathy genes (all discovery years; this is NOT a temporal back-test) are enriched in
specific cell types, compared with simple expression baselines.

Usage: python workflow/scripts/pilot.py [config/config.yaml]
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from pnprank.atlases import collapse_to_ensembl, load_mtx, symbol_to_ensembl
from pnprank.backtest import auroc, bootstrap_ci
from pnprank.baselines import within_bin_percentile
from pnprank.pseudobulk import cpm, mean_over_donors, pseudobulk
from pnprank.specificity import fraction_specificity, tau

cfg = yaml.safe_load(Path(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
atlas = cfg["pilot"]["atlas"]
seed = cfg["seeds"]["global"]
nboot = cfg["bootstrap"]["n_resamples"]
out = Path("results/pilot")
out.mkdir(parents=True, exist_ok=True)
interim = Path(paths["data_root"]) / "interim" / atlas

# --- load and map genes -------------------------------------------------------
m, symbols, cells = load_mtx(interim / "counts")
hgnc_path = Path(cfg["truth"]["sources_dir"]) / "hgnc_complete_set.txt"
mapping, how = symbol_to_ensembl(hgnc_path)
m, ens, unmapped = collapse_to_ensembl(m, symbols, mapping)
meta_path = next((Path(paths["data_root"]) / "raw" / atlas).glob("*snRNA-seq_human_barcode_meta.csv.gz"))
meta = pd.read_csv(meta_path).set_index("V1").loc[cells]
how_counts = pd.Series([how.get(s, "unmapped") for s in symbols]).value_counts().to_dict()

# --- pseudobulk and specificity --------------------------------------------------
pb = cfg["pseudobulk"]
sums, ncells = pseudobulk(
    m, ens, meta[cfg["pilot"]["cell_type_column"]].reset_index(drop=True),
    meta[cfg["pilot"]["donor_column"]].astype(str).reset_index(drop=True),
    pb["min_cells_per_group"],
)
expr_d = cpm(sums)
expr = mean_over_donors(expr_d, pb["min_donors"])
spec = fraction_specificity(expr)
tau_s = tau(expr)
ncells.to_csv(out / "pseudobulk_groups.tsv", sep="\t")

# donor consistency: Spearman of specificity between donor pairs, per cell type
donors = sorted(expr_d.columns.get_level_values("donor").unique())
cons = []
for i, a in enumerate(donors):
    for b in donors[i + 1:]:
        ca = expr_d.xs(a, axis=1, level="donor")
        cb = expr_d.xs(b, axis=1, level="donor")
        common = [c for c in ca.columns if c in cb.columns and c in expr.columns]
        ta, tb = fraction_specificity(ca[common]), fraction_specificity(cb[common])
        for c in common:
            x = pd.concat([ta[c], tb[c]], axis=1).dropna()
            cons.append({"donor_a": a, "donor_b": b, "cell_type": c, "n_genes": len(x),
                         "spearman": round(x.rank().corr().iloc[0, 1], 3)})
pd.DataFrame(cons).to_csv(out / "donor_consistency.tsv", sep="\t", index=False)

# --- labels -------------------------------------------------------------------
hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str)
pc = set(hgnc.loc[(hgnc.status == "Approved") & (hgnc.locus_group == "protein-coding gene"),
                  "ensembl_gene_id"].dropna())
truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
universe = [g for g in expr.index if g in pc]
pos_all = set(truth.ensembl_gene_id)
neg = set(universe) - set(cands.ensembl_gene_id)  # drop every candidate, incl. unverified
missing_truth = sorted(set(truth.hgnc_symbol[~truth.ensembl_gene_id.isin(universe)]))


def evaluate(score: pd.Series, positives: set[str], name: str, seed_offset: int) -> dict:
    genes = [g for g in universe if g in positives or g in neg]
    s = score.reindex(genes).to_numpy(dtype=float)
    y = np.array([g in positives for g in genes])
    est, lo, hi = bootstrap_ci(s, y, auroc, nboot, cfg["bootstrap"]["ci"], seed + seed_offset)
    return {"score": name, "n_pos": int(y.sum()), "n_neg": int((~y).sum()),
            "auroc": round(est, 3), "ci_low": round(lo, 3), "ci_high": round(hi, 3)}


baselines = {
    "mean_cpm_all_types": expr.mean(axis=1),
    "max_cpm_any_type": expr.max(axis=1),
    "tau": tau_s,
    "max_specificity": spec.max(axis=1),
}
sets = {"all": pos_all, "tier1": set(truth.ensembl_gene_id[truth.tier == "1"]),
        "tier2": set(truth.ensembl_gene_id[truth.tier == "2"])}
for cat in sorted(truth.disease_category.unique()):
    sets[f"cat_{cat}"] = set(truth.ensembl_gene_id[truth.disease_category == cat])

rows = []
k = 0
for set_name, positives in sets.items():
    for name, score in baselines.items():
        rows.append({"positives": set_name, **evaluate(score, positives, name, k)})
        k += 1
    for c in spec.columns:
        rows.append({"positives": set_name, **evaluate(spec[c], positives, f"specificity:{c}", k)})
        k += 1
    # expression-matched: specificity ranked within deciles of mean expression
    level = expr.mean(axis=1).reindex(universe)
    for c in ["max_specificity", *spec.columns]:
        sc = baselines[c] if c == "max_specificity" else spec[c]
        matched = within_bin_percentile(sc.reindex(universe), level, n_bins=10)
        rows.append({"positives": set_name, **evaluate(matched, positives, f"expr_matched:{c}", k)})
        k += 1
res = pd.DataFrame(rows)
res.to_csv(out / "auroc_by_celltype.tsv", sep="\t", index=False)

# top truth genes per cell type (sanity listing)
sym = dict(zip(truth.ensembl_gene_id, truth.hgnc_symbol))
top = []
for c in spec.columns:
    s = spec.loc[[g for g in spec.index if g in pos_all], c].dropna().sort_values(ascending=False)
    top.append({"cell_type": c, "top_truth_genes": ", ".join(f"{sym[g]}({v:.2f})" for g, v in s.head(10).items())})
pd.DataFrame(top).to_csv(out / "top_truth_genes_by_celltype.tsv", sep="\t", index=False)

json.dump({
    "atlas": atlas, "n_cells": len(cells), "n_symbols": len(symbols),
    "symbol_mapping": how_counts, "n_ensembl_after_collapse": len(ens),
    "n_universe_protein_coding": len(universe), "n_negatives": len(neg & set(universe)),
    "truth_in_universe": len(pos_all & set(universe)), "truth_missing": missing_truth,
    "cell_types_kept": list(spec.columns),
}, open(out / "summary.json", "w"), indent=1)
print(json.load(open(out / "summary.json")))

# --- figure -------------------------------------------------------------------
from pnprank.plotting import auroc_heatmap  # noqa: E402

order_sets = ["all", "tier1", "tier2", "cat_CMT", "cat_dHMN", "cat_HSAN", "cat_HSP", "cat_SYND"]
cts = ["Schwann_M", "Schwann_N", "Satglia", "NF1", "NF2", "NF3", "PEP", "NP", "SST", "TRPM8",
       "cLTMR", "Fibroblast_Apod", "Fibroblast_Dcn", "Fibroblast_Mgp", "Vascular", "Immune"]
sub = res[res.score.str.startswith("expr_matched:") & ~res.score.str.endswith("max_specificity")].copy()
sub["cell_type"] = sub.score.str.split(":").str[1]
piv = {k: sub.pivot(index="positives", columns="cell_type", values=k).loc[order_sets, cts]
       for k in ("auroc", "ci_low", "ci_high")}
Path("results/figures").mkdir(parents=True, exist_ok=True)
auroc_heatmap(piv["auroc"], piv["ci_low"], piv["ci_high"],
              "Pilot: cell-type specificity of known neuropathy genes, expression-matched "
              "(Yang 2022 human TG)", "results/figures/pilot_specificity_auroc.png")
