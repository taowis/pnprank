"""Case study: rank genes inside published linkage intervals (config/loci.tsv).

For each locus, genes overlapping the cytoband interval (GRCh38, UCSC cytoBand; Ensembl
116 gene coordinates) are ranked by the back-test scores at the locus's cutoff year: the
model was trained on tier-1 genes known before that year only, and every interval gene is
scored out of fold. Known tier-1 genes at the cutoff are not scored (they were training
positives) and are listed separately. This is a demonstration on hand-picked loci, not a
validation; the selection of loci is described in the Stage report.

Usage: python workflow/scripts/locus_case_study.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from pnprank.backtest import oof_scores_at_cutoff, training_labels, within_fold_percentile
from pnprank.features import label_features, load_celltype_map, load_profiles

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
root = Path(paths["data_root"])
ext = root / "external" / "loci"
out = Path("results/case_studies")
out.mkdir(parents=True, exist_ok=True)
scores_dir = Path("results/backtest/scores")
SCORES = ["full", "non_atlas", "gtex_nerve", "gtex_mean", "neg_loeuf"]

loci = pd.read_csv("config/loci.tsv", sep="\t", dtype=str)
bands = pd.read_csv(ext / "cytoBand.hg38.txt.gz", sep="\t", header=None,
                    names=["chrom", "start", "end", "band", "stain"])
REST = "https://rest.ensembl.org"


def _get(url: str):
    import json
    import time
    import urllib.request

    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed: {url}")


def genes_in(locus_id: str, chrom: str, s: int, e: int, target: str) -> pd.DataFrame:
    """Protein-coding genes overlapping chrom:s-e plus the target gene, from the Ensembl
    REST API (GRCh38; release recorded in the cache), cached under the data root."""
    f = ext / f"genes_{locus_id}.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t", dtype={"chrom": str})
    rows = []
    for a in range(s, e, 5_000_000):     # the overlap endpoint allows 5 Mb per request
        b = min(a + 5_000_000 - 1, e)
        for g in _get(f"{REST}/overlap/region/human/{chrom}:{a}-{b}?feature=gene"
                      "&biotype=protein_coding"):
            rows.append((g["id"], g["seq_region_name"], g["start"], g["end"],
                         g.get("external_name", ""), True))
    t = _get(f"{REST}/lookup/symbol/homo_sapiens/{target}")
    rows.append((t["id"], t["seq_region_name"], t["start"], t["end"], target,
                 t["seq_region_name"] == chrom and t["end"] >= s and t["start"] <= e))
    release = _get(f"{REST}/info/data")["releases"][0]
    df = pd.DataFrame(rows, columns=["ensembl_gene_id", "chrom", "start", "end", "symbol",
                                     "in_interval"]).drop_duplicates("ensembl_gene_id")
    df["ensembl_release"] = release
    df.to_csv(f, sep="\t", index=False)
    return df

truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)


def interval(chrom: str, b0: str, b1: str) -> tuple[int, int]:
    """Genomic span from the first sub-band of b0 to the last sub-band of b1."""
    b = bands[bands.chrom == f"chr{chrom}"].sort_values("start")
    def hits(label):
        # sub-bands of a label: q21 -> q21.1, q21.31; q22.3 -> q22.31, q22.32
        rest = b.band.str[len(label):]
        sub = rest.str.startswith(".") | (("." in label) & rest.str.match(r"^\d"))
        m = b[(b.band == label) | (b.band.str.startswith(label) & sub)]
        if m.empty:
            raise ValueError(f"band {chrom}{label} not found")
        return m
    s = min(hits(b0).start.min(), hits(b1).start.min())
    e = max(hits(b0).end.max(), hits(b1).end.max())
    return int(s), int(e)


def scores_at(cutoff: int) -> pd.DataFrame:
    """Cached Stage 5 scores at a cutoff; computed with the same code if missing."""
    f = scores_dir / f"tier1_{cutoff}.tsv.gz"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0)
    X = pd.read_csv(root / "processed" / "features_static.tsv.gz", sep="\t", index_col=0)
    U = list(X.index)
    fc, mc = cfg["features"], cfg["model"]
    prof = load_profiles(fc["atlases"], root / "interim", load_celltype_map(fc["celltype_map"]),
                          cfg["pseudobulk"]["min_cells_per_group"], cfg["pseudobulk"]["min_donors"],
                         root / "processed" / "features_static.tsv.gz")
    non_atlas = [c for c in X.columns if c.startswith(("gnomad|", "gtex|"))]
    y = training_labels(U, truth, cutoff, {"1"})
    seed, folds = cfg["seeds"]["global"], cfg["backtest"]["scoring_folds"]
    full, _, fold = oof_scores_at_cutoff(X, y, "logistic", seed, folds, mc["inner_folds"],
                                         lambda k: label_features(prof, U, k), 5)
    non, _, _ = oof_scores_at_cutoff(X[non_atlas], y, "logistic", seed, folds,
                                     mc["inner_folds"], None, 5)
    df = pd.DataFrame({"fold": fold, "full": full, "non_atlas": non})
    df["gtex_nerve"] = X["gtex|Nerve_Tibial|logtpm"].reindex(df.index)
    df["gtex_mean"] = X["gtex|mean_all|logtpm"].reindex(df.index)
    df["neg_loeuf"] = -X["gnomad|loeuf"].reindex(df.index)
    df.to_csv(f, sep="\t", compression="gzip")
    print(f"computed scores at {cutoff}: n_pos_train={int(y.sum())}")
    return df


rows, tops = [], []
for _, L in loci.iterrows():
    cutoff = int(L.cutoff_year)
    df = scores_at(cutoff)
    pct = pd.DataFrame({s: within_fold_percentile(df[s], df["fold"]) for s in SCORES})
    s, e = interval(L.chrom, L.band_start, L.band_end)
    genes = genes_in(L.locus_id, L.chrom, s, e, L.gene)
    g = genes[genes.in_interval]
    target = genes[genes.symbol == L.gene].ensembl_gene_id.iloc[0]
    ids = list(genes.ensembl_gene_id)
    in_interval = target in set(g.ensembl_gene_id)
    known = sorted(truth.hgnc_symbol[truth.ensembl_gene_id.isin(ids) & (truth.tier == "1")
                                     & (truth.discovery_year.astype(int) < cutoff)])
    scored = [i for i in ids if i in pct.index]
    p = pct.loc[scored]
    row = {"locus_id": L.locus_id, "interval": f"{L.chrom}{L.band_start}-{L.band_end}",
           "grch38_start": s, "grch38_end": e, "grch38_span_mb": round((e - s) / 1e6, 1), "gene": L.gene, "cutoff_year": cutoff,
           "gene_in_interval": in_interval, "n_protein_coding": len(g),
           "n_scored": len(scored), "known_tier1_not_scored": ",".join(known),
           "genome_pct_full": round(float(p.at[target, "full"]), 3)}
    for sc in SCORES:
        r = p[sc].rank(ascending=False, method="min")
        row[f"rank_{sc}"] = int(r[target]) if not np.isnan(p.at[target, sc]) else None
    rows.append(row)
    sym = dict(zip(genes.ensembl_gene_id, genes.symbol))
    top = p.sort_values("full", ascending=False).head(10)
    for k, (gid, v) in enumerate(top.iterrows(), 1):
        tops.append({"locus_id": L.locus_id, "rank": k, "gene": sym.get(gid, gid),
                     "is_target": gid == target, "pct_full": round(v.full, 3),
                     "pct_non_atlas": round(v.non_atlas, 3)})

res = pd.DataFrame(rows)
res.to_csv(out / "locus_ranking.tsv", sep="\t", index=False)
pd.DataFrame(tops).to_csv(out / "locus_top10.tsv", sep="\t", index=False)
print(res.to_string(index=False))
