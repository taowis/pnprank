"""How many known neuropathy genes are expressed in tissues used for patient RNA-seq.

Patient RNA-seq is usually done on blood or cultured fibroblasts. A gene not expressed in
the sampled tissue cannot show aberrant expression or splicing there. This counts tier-1
genes with GTEx v10 median TPM at or above a threshold in each tissue, next to tibial nerve,
skeletal muscle and cervical spinal cord. TPM >= 1 is a common detection cut-off; 5 and 10
are shown as well.

Usage: python workflow/scripts/patient_tissue_expression.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

paths = yaml.safe_load(Path("config/paths.yaml").read_text())
gct = (Path(paths["data_root"]) / "external" / "gtex"
       / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
TISSUES = ["Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts",
           "Muscle_Skeletal", "Nerve_Tibial", "Brain_Spinal_cord_cervical_c-1"]
THRESHOLDS = [1, 5, 10]

tpm = pd.read_csv(gct, sep="\t", skiprows=2)
tpm["ensembl_gene_id"] = tpm["Name"].str.split(".").str[0]
tpm = tpm.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")[TISSUES]

truth = pd.read_csv("data/truth/neuropathy_genes.tsv", sep="\t")
tier1 = truth[truth["tier"] == 1]
missing = sorted(set(tier1["ensembl_gene_id"]) - set(tpm.index))
t = tpm.reindex(tier1["ensembl_gene_id"]).dropna()

rows = []
for tissue in TISSUES:
    row = {"tissue": tissue, "n_genes": len(t)}
    for th in THRESHOLDS:
        row[f"n_tpm_ge_{th}"] = int((t[tissue] >= th).sum())
    rows.append(row)
summary = pd.DataFrame(rows)

out = Path("results/case_studies")
out.mkdir(parents=True, exist_ok=True)
summary.to_csv(out / "patient_tissue_expression.tsv", sep="\t", index=False)

per_gene = tier1.set_index("ensembl_gene_id")[["hgnc_symbol"]].join(t, how="inner")
per_gene.round(3).to_csv(out / "patient_tissue_expression_genes.tsv", sep="\t")

print(summary.to_string(index=False))
print("tier-1 genes without a GTEx row:", missing)
blood_fib_low = per_gene[(per_gene["Whole_Blood"] < 1) & (per_gene["Cells_Cultured_fibroblasts"] < 1)]
print(f"below 1 TPM in both whole blood and fibroblasts: {len(blood_fib_low)}; of these, "
      f">= 1 TPM in tibial nerve: {(blood_fib_low['Nerve_Tibial'] >= 1).sum()}")
print(blood_fib_low.round(2).to_string())
