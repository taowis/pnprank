# Known neuropathy genes in tissues used for patient RNA-seq

**Date:** 2026-10-16. Script: `workflow/scripts/patient_tissue_expression.py`; outputs:
`results/case_studies/patient_tissue_expression.tsv` (counts) and
`patient_tissue_expression_genes.tsv` (median TPM per tier-1 gene and tissue). Input: GTEx
v10 gene median TPM (`config/config.yaml`, `external.gtex_median_tpm`).

Patient RNA-seq for inherited neuropathy is usually done on blood or cultured fibroblasts.
A gene that is not expressed in the sampled tissue cannot show aberrant expression or
splicing there. This counts how many of the 100 tier-1 genes reach a median TPM threshold
in each tissue. All 100 have a GTEx row.

| GTEx tissue | TPM >= 1 | TPM >= 5 | TPM >= 10 |
|---|---|---|---|
| Whole blood | 67 | 54 | 35 |
| EBV-transformed lymphocytes | 74 | 63 | 59 |
| Cultured fibroblasts | 84 | 72 | 62 |
| Skeletal muscle | 75 | 55 | 45 |
| Tibial nerve | 92 | 85 | 78 |
| Spinal cord (cervical C1) | 89 | 76 | 60 |

- **15 tier-1 genes are below 1 TPM in both whole blood and cultured fibroblasts.**
- **Seven of them are at least 1 TPM in tibial nerve, mostly Schwann cell genes:** GJB1, SH3TC2, PMP2, DRP2, EGR2 and KIF1A. Their top cell type in the ranking is a Schwann cell population, and GDAP1 is the seventh.
- **SCN10A, SCN11A and NTRK1 are below 1 TPM in tibial nerve as well.** They act in sensory neuron cell bodies, which tibial nerve samples do not contain. In the ranking, all three are most specific to sensory neurons in the human DRG atlas.
- **SLC5A7 is the same case for motor neurons.** It is near zero in every tissue listed, and it is most specific to motor neurons in the human spinal cord atlas.
- **So most known genes are detectable in fibroblasts (84 of 100), but a fifth are missed by blood and fibroblasts alike or are barely detectable in both.** Several of these are classic CMT genes.

This supports describing the atlases as a complement to patient RNA-seq: for genes like these, the atlases show which cell types of nerve, ganglion or spinal cord express them. It is not a test of PNPrank's ranking.

Caveats:
- **TPM >= 1 is a conventional cut-off, not a measure of whether an outlier could be detected.** Detection also depends on depth and on the method.
- **GTEx tissues are from adult post-mortem donors.** Patient samples, especially cultured fibroblasts, differ in culture conditions and passage.
