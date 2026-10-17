"""GRCh38 coordinates (Ensembl 116 GTF) for the genes in the released ranking, used by
`pnprank interval`. The table is small and committed so the tool works without the data
root. Coordinates are 1-based and inclusive, as in the GTF.

Usage: python workflow/scripts/gene_coordinates.py
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
h = cfg["external"]["human_gtf"]
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
gtf = Path(paths["data_root"]) / Path(h["dest"]).relative_to("data")

rows = []
with gzip.open(gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "gene":
            continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        rows.append((a["gene_id"], a.get("gene_name", ""), f[0], int(f[3]), int(f[4]), f[6],
                     a.get("gene_biotype", "")))
genes = pd.DataFrame(rows, columns=["ensembl_gene_id", "gtf_gene_name", "chrom", "start", "end",
                                    "strand", "biotype"]).set_index("ensembl_gene_id")

ranking = pd.read_csv("results/tables/pnprank_ranking.tsv.gz", sep="\t",
                      usecols=["ensembl_gene_id", "hgnc_symbol"])
out = ranking.join(genes, on="ensembl_gene_id", how="left")
missing = out[out.chrom.isna()]
out = out.dropna(subset=["chrom"]).astype({"start": int, "end": int})
out = out.sort_values(["chrom", "start"])
out.to_csv(h["coordinates_table"], sep="\t", index=False, compression="gzip")

print(f"{len(out)} of {len(ranking)} ranked genes have Ensembl 116 coordinates")
print(f"not in the Ensembl 116 GTF: {len(missing)} {list(missing.hgnc_symbol)[:20]}")
print("chromosomes:", sorted(out.chrom.unique(), key=lambda c: (len(c), c)))
print("biotypes:", out.biotype.value_counts().to_dict())
diff = out[out.gtf_gene_name.fillna("") != out.hgnc_symbol.fillna("")]
print(f"symbol differs from GTF gene_name: {len(diff)}")
