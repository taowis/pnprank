# Stage 3 — Feasibility pilot on one atlas

**Date:** 2026-10-07
**Atlas:** Yang 2022 human trigeminal ganglion (GEO GSE197289, PMID 35349784).
snRNA-seq, 38,028 nuclei, 3 donors, 16 author cell types.

**Outputs:**
- `results/pilot/`:
  - `summary.json`
  - `auroc_by_celltype.tsv`
  - `donor_consistency.tsv`
  - `pseudobulk_groups.tsv`
  - `top_truth_genes_by_celltype.tsv`
- `results/figures/pilot_specificity_auroc.png`
- Script: `workflow/scripts/pilot.py`

This is **not** a back-test. It uses all discovery years and asks only whether the
pipeline works and whether cell-type expression carries any signal about known
neuropathy genes. No multiple-testing correction is applied; treat results as
exploratory.

## What was done

1. Added a guarded download command (`pnprank download-atlas`). It enforces the 5 GB
   single-file limit and the disk budget, and writes size and SHA-256 to
   `config/atlases.yaml`. Downloaded the human snRNA raw counts (182 MB) and metadata.
2. Converted the author RDS to Matrix Market with R (added `r-base` and `r-matrix`
   to the environment; lockfile regenerated).
3. Mapped symbols to Ensembl 116 through HGNC.
4. Built cell type × donor pseudobulk (at least 20 cells per group, at least 2 donors
   per cell type) and computed EWCE-style fraction specificity and tau.
5. Scored truth genes against protein-coding background genes (all 249 candidates
   removed from negatives), using AUROC with 1,000-resample stratified bootstrap 95% CIs.
6. Added an expression-matched control: specificity ranked within deciles of mean
   expression.
7. Added unit tests for pseudobulk, specificity, AUROC/bootstrap, symbol mapping and
   the matching control (18 passed).

## What was found

**Mapping:**
- 21,003 symbols map as approved, 804 via previous symbols and 182 via aliases;
  10,225 are unmapped (almost all clone-named lncRNAs).
- 17,979 protein-coding genes are in the universe.
- 240 of the 242 truth genes are present; MT-RNR1 and MT-TL1 are non-coding.

**Expression level beats specificity (all positives, n = 240):**

| Score | AUROC [95% CI] |
|---|---|
| Mean CPM across cell types | 0.68 [0.65, 0.71] |
| Mean CPM, tier 1 only | 0.74 [0.70, 0.78] |
| Mean CPM, HSAN | 0.83 [0.76, 0.90] |
| Tau (specificity) | 0.36 [0.33, 0.39] |
| Max fraction specificity | 0.39 [0.36, 0.43] |
| Max specificity, expression-matched | 0.47 [0.43, 0.51] |
| Max specificity, expression-matched, tier 1 | 0.56 [0.50, 0.62] |

- Known neuropathy genes are more highly expressed and **less** cell-type-specific
  than background genes.
- Most of the raw neuron-specificity signal (AUROC about 0.6) is explained by
  expression level. After matching, the all-gene neuron signals drop to about 0.55.

**Category-specific signals survive expression matching, in the expected direction:**

| Category | Cell type | AUROC [95% CI] |
|---|---|---|
| CMT (n=59) | myelinating Schwann | 0.61 [0.53, 0.69] |
| HSAN (n=16) | SST neurons | 0.71 [0.58, 0.82] |
| HSAN | PEP neurons | 0.67 [0.54, 0.80] |
| dHMN (n=19) | NF1 neurons | 0.71 [0.58, 0.83] |
| dHMN | PEP neurons | 0.70 [0.56, 0.81] |
| HSP (n=24) | NP neurons | 0.69 [0.58, 0.79] |

- Every category is depleted in fibroblasts and vascular cells (AUROC 0.25–0.46).
- The most Schwann-specific truth genes are the classic myelin genes (FA2H, MAG, PRX,
  DRP2, MPZ, EGR2, GJB1, SH3TC2). Sodium channel genes SCN9A, SCN10A and SCN11A top
  the nociceptor-like subtypes. This is a sanity check that labels and mapping are right.

**Donor consistency:**
- Spearman correlation of per-donor specificity is 0.63–0.77 between donors 1 and 2.
- It is 0.45–0.66 for comparisons involving donor 3, which has only 142 neurons.
- Schwann_N is the least consistent cell type.

## Surprises and problems

- **The main pitfall is expression level, not specificity.** A model on cell-type
  specificity alone would underperform a trivial "how highly is this gene expressed"
  baseline. Stage 4/5 must therefore:
  - include expression-level baselines, both tissue-specific and tissue-agnostic, for
    example GTEx-wide expression;
  - report whether nerve-atlas features add anything beyond them.
- **Tier 2 behaves differently.** Syndromic genes are mostly broadly expressed
  metabolic or housekeeping genes, with weaker cell-type signal. This supports making
  tier 1 the primary analysis set, with tier 1+2 as a sensitivity analysis.
- **Small category sizes** (HSAN 16, dHMN 19) give wide CIs, and the many category ×
  cell type comparisons are uncorrected.
- **Cranial ganglion only:** this atlas has no motor neurons and no peripheral nerve
  trunk. The dHMN signal in NF1 neurons is therefore a proxy (large-diameter sensory
  neurons), not motor neurons.
- **Data quirk:** the GEO RDS file is gzip-compressed twice; handled in the converter.

## Feasibility verdict

The pipeline works end to end: download with checksums, author labels, gene mapping,
pseudobulk, specificity and bootstrapped evaluation. Cell-type expression carries
biologically sensible signal per neuropathy category. The overall signal is modest,
and much of it is expression level.

## Proposed next step

**Stage 4: full feature set and model.**
1. Add the remaining shortlisted atlases, in priority order:
   - Heming sural nerve controls;
   - Alkaslasi and Yadav for motor neurons;
   - Jung DRG;
   - Yim and Stefanovska mouse nerve;
   - muscle/NMJ atlases.
2. Harmonise cell types across atlases.
3. Add non-atlas features: gnomAD constraint and GTEx tissue expression, both as
   features and as baselines.
4. Add label-derived features (co-expression with known genes) computed per cutoff
   year only.
5. Extend `tests/test_leakage.py` to cover every feature before any back-test.

Decisions needed at the start of Stage 4:
- whether tier 1 is the primary positive set;
- approval for the 21 GB iPain file, if it is used.
