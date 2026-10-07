# Stage 4 — Full feature set and model

**Date:** 2026-10-07

**Outputs:**
- `results/tables/stage4/`:
  - `cv_performance.tsv`: AUROC and average precision with 95% bootstrap CIs.
  - `paired_auroc_differences.tsv`: paired bootstrap differences.
  - `oof_*.tsv.gz`: out-of-fold scores per model and positive set.
  - `cv_best_params.json`: inner-CV hyperparameters per outer fold.
  - `logistic_full_coefficients.tsv`.
  - `ranking_current_logistic_full.tsv.gz`: all 19,228 genes.
  - `summary.json`.
  - `yadav_motoneuron_qc.tsv`.
- Scripts:
  - `workflow/scripts/standardize.py`: per-atlas readers.
  - `workflow/scripts/seurat_to_mtx.R`.
  - `workflow/scripts/stage4_model.py`.
- Package: `features.py`, `model.py` and `coexpression.py` in `src/pnprank`.
- Log of the reported run: `logs/stage4_model_run3.log` (not committed).

**This is not the temporal back-test.** It uses all currently known genes and asks
one question: do single-cell atlas features add signal beyond simple non-atlas
predictors? The back-test (Stage 5) is the test of the central claim.

## What was done

1. **Standardised 10 atlases.** Each atlas was reduced to raw counts, author labels,
   Ensembl IDs and pseudobulk (cell type × donor). Mouse genes were mapped to human
   by Ensembl Compara 116 one-to-one orthologues.

   | Atlas | Species | Tissue | Cells* | Donors/replicates | Harmonised types used |
   |---|---|---|---|---|---|
   | yang2022_human_tg | human | trigeminal ganglion | 38,028 | 3 | 12 |
   | heming2025_human_sural | human | sural nerve (controls) | 33,732 | 4 | 10 |
   | sparc465_human_sympathetic | human | sympathetic ganglia | 52,950 | 8 libraries | 9 |
   | yadav2023_human_spinal_cord | human | spinal cord | 55,040 | 7 | 12 |
   | kedlian2024_human_muscle | human | skeletal muscle (nuclei) | 92,259 | 13 | 15 |
   | yang2022_mouse_tg | mouse | trigeminal ganglion (naive) | 30,821 | 7 | 12 |
   | stefanovska2025_mouse_sciatic | mouse | sciatic nerve (young) | 7,888 | 3 | 10 |
   | krauter2025_mouse_drg | mouse | DRG neurons | 4,714 | 4 | 6 |
   | alkaslasi2021_mouse_mn | mouse | spinal cholinergic neurons | 16,042 | 9 | 5 |
   | ham2025_mouse_nmj | mouse | muscle / NMJ | 36,191 | 4 | 12 |

   *Cells in standardised author groups, before cell-type exclusions and the
   20-cell / 2-donor thresholds.

2. **Harmonised author labels** to a shared vocabulary in `config/celltype_map.tsv`.
   Every author label has an explicit entry, either a harmonised type or "exclude",
   and unmapped labels raise an error.

3. **Built features for 19,228 protein-coding genes.** These are the genes seen in at
   least one atlas or in GTEx.
   - **Atlas features (label-free):**
     - log1p mean CPM and fraction specificity per atlas × cell type;
     - tau per atlas.
   - **Non-atlas features:**
     - gnomAD v4.1 LOEUF, pLI and missense Z;
     - GTEx v10 log TPM for tibial nerve, skeletal muscle and cervical spinal cord;
     - mean GTEx expression across tissues and GTEx tau.
   - **Label-derived features:** co-expression with known genes, one per atlas, with
     the gene itself left out.
   - Total: 224 static features, 234 with label-derived features.

4. **Evaluated models by nested cross-validation.**
   - 5 outer folds, 3 inner folds.
   - Label-derived features are recomputed from each outer fold's training positives
     only, and hyperparameters are tuned inside the training fold
     (`tests/test_leakage.py::test_cv_label_features_use_training_fold_only`).
   - Models:
     - L2 logistic regression (median imputation, missing indicators, scaling,
       balanced weights);
     - histogram gradient boosting.
   - Each model was fitted with non-atlas features only and with all features.
   - Single-feature baselines: GTEx mean expression, GTEx tibial nerve expression,
     and −LOEUF.

5. **Labels.**
   - Positives:
     - Primary: tier 1 (100 genes).
     - Sensitivity analysis: tier 1+2 (240 genes in the universe).
   - Negatives: protein-coding genes that are not candidates (18,981).
   - Other candidates are dropped from both sets.
   - AUROC and average precision have 1,000-resample stratified bootstrap 95% CIs.
     Model comparisons use paired bootstrap on identical resamples.

6. **Fitted a final logistic model** on all tier 1 genes. It gives the current ranking
   and the coefficients.

7. **Tests:** 24 passed and 1 skipped (the Stage 5 end-to-end back-test test).

## What was found

**Atlas features add signal beyond the non-atlas predictors**

Performance (out-of-fold):

| Positives | Score | AUROC [95% CI] | Avg. precision [95% CI] |
|---|---|---|---|
| tier 1 | GTEx mean expression | 0.672 [0.620, 0.724] | 0.010 [0.008, 0.013] |
| tier 1 | GTEx tibial nerve | 0.705 [0.654, 0.754] | 0.013 [0.010, 0.022] |
| tier 1 | −LOEUF | 0.649 [0.600, 0.701] | 0.008 [0.007, 0.011] |
| tier 1 | logistic, non-atlas | 0.739 [0.691, 0.785] | 0.022 [0.015, 0.045] |
| tier 1 | **logistic, full** | **0.791 [0.750, 0.829]** | **0.060 [0.029, 0.111]** |
| tier 1 | HGB, non-atlas | 0.685 [0.634, 0.732] | 0.011 [0.009, 0.015] |
| tier 1 | HGB, full | 0.749 [0.701, 0.795] | 0.021 [0.014, 0.036] |
| tier 1+2 | GTEx tibial nerve | 0.671 [0.643, 0.701] | 0.024 [0.020, 0.037] |
| tier 1+2 | logistic, non-atlas | 0.708 [0.679, 0.737] | 0.030 [0.025, 0.042] |
| tier 1+2 | **logistic, full** | **0.762 [0.735, 0.790]** | **0.068 [0.046, 0.100]** |
| tier 1+2 | HGB, full | 0.734 [0.705, 0.764] | 0.043 [0.032, 0.064] |

Prevalence is 0.0052 (tier 1) and 0.0125 (tier 1+2).

Paired AUROC differences:

| Positives | Comparison | ΔAUROC [95% CI] |
|---|---|---|
| tier 1 | logistic full − non-atlas | +0.052 [0.011, 0.094] |
| tier 1 | HGB full − non-atlas | +0.064 [0.017, 0.116] |
| tier 1 | logistic full − HGB full | +0.041 [0.005, 0.076] |
| tier 1+2 | logistic full − non-atlas | +0.054 [0.028, 0.084] |
| tier 1+2 | HGB full − non-atlas | +0.044 [0.010, 0.078] |
| tier 1+2 | logistic full − HGB full | +0.028 [0.005, 0.052] |

- **The atlas gain is consistent.** All six intervals exclude zero, in both model
  families and both positive sets. The gain is modest, about 0.05 AUROC. For tier 1
  its lower bound is close to zero (0.011).
- **Logistic beats gradient boosting.** With 100 positives a strongly regularised
  linear model is more stable. The inner CV chose C = 0.001 for the final model.
- **Top-k enrichment (tier 1, out-of-fold):**
  - top 100 genes (0.5% of the universe): 11 of 100 positives, about 20× the 0.5
    expected by chance;
  - top 1,000: 33 of 100.
  - Average precision is low in absolute terms (0.06). Most genes ranked high are
    not known neuropathy genes.
- **Six atlases versus ten.** An earlier run with six atlases gave logistic full
  0.793 [0.755, 0.830]. The four additional atlases add no measurable discrimination:
  - Yadav spinal cord;
  - Yang mouse TG;
  - Stefanovska sciatic nerve;
  - Krauter DRG.

  They are kept because they are needed for category-specific analyses: motor neurons
  for dHMN, mouse nerve and DRG subtypes.

**Largest standardised coefficients (final logistic model) are biologically
sensible**

- Positive:
  - gnomAD missense Z (+0.34);
  - motor-neuron specificity in human spinal cord (+0.23);
  - SST and PEP sensory-neuron expression in human TG (+0.20, +0.18);
  - myelinating-Schwann specificity in mouse TG (+0.20);
  - C-LTMR specificity in mouse DRG (+0.18);
  - Schwann specificity at the mouse NMJ (+0.16).
- Negative: specificity for fibroblasts, satellite glia, immune cells, OPCs and spinal
  interneurons.
- **Caveat:** these are coefficients of correlated features under strong
  regularisation. They show what the model uses, not independent effects.

**The current ranking still contains an expression-level confound**

- **Known genes at the top.** The top 10 includes MPZ, PMP22, PRX, DYNC1H1, NEFH and
  NEFL (all tier 1).
- **Top 100 composition:**
  - 12 tier 1 genes;
  - 8 tier 2 genes;
  - 1 unverified candidate (L1CAM).
- **Housekeeping and cytoskeletal genes also rank high:**
  - TUBA1A, TUBA1B and TUBB;
  - ACTB, ACTG1 and TMSB4X;
  - YWHAG and CALM2.
- **Interpreting the top negatives.** Many of the other top-ranked negatives are
  myelin or neuronal genes, for example MBP, UGT8, LGI4, CDH19, SCN7A and KCNA1. They
  are plausible by expression pattern, which is what the model can see. Nothing here
  is evidence that they cause neuropathy.
- **Ranking scope.** This ranking is trained on current labels and is shown for
  inspection only. It is not a validated prediction.

## Surprises and problems

- **Expression level is still a confound.** It was found in the pilot and it persists.
  The model has no feature that separates "highly expressed in nerve/neuron" from
  "broadly highly expressed".
  - Stage 5 should add an expression-matched evaluation, as in the pilot, and report
    it next to the raw AUROC.
  - A candidate model change is residualising atlas expression on GTEx mean
    expression. It would be pre-specified before the back-test, not tuned on it.
- **More atlases did not help overall discrimination.** Gains from the four added
  atlases, if any, may show up only within categories. Stage 5 should report
  per-category results (CMT, dHMN, HSAN, HSP).
- **Jung 2023 human DRG is not included.**
  - Its counts and author labels are only in one 4.2 GB Seurat RDS. Loading it needs
    an estimated 15–25 GB of RAM, more than this machine has free.
  - Option to be checked: the 240 MB cross-species object
    (GSE201654_mousegpcynohumanDRG.rds.gz) together with the per-sample matrices
    (GSE201586_RAW.tar, about 2 GB). Pseudobulk sums are additive, so the per-sample
    route would be exact.
  - Not yet verified: whether that object holds human labels, and whether its
    barcodes match the per-sample matrices.
- **Yim 2022 mouse nerves excluded.** Only numeric clusters were deposited.
- **Yadav donors 1 and 5 fail the motoneuron QC.** Donor 5's "motoneurons" look like
  glial contamination.
- **Memory pressure.** One earlier model run thrashed under memory pressure and was
  killed. The reported run (run 3) completed normally.

## Leakage status

- `tests/test_leakage.py` passes:
  - post-cutoff invariance of label-derived features;
  - the sensitivity check;
  - leave-self-out co-expression;
  - training-fold-only label features in CV.
- The end-to-end back-test test is still skipped and becomes active in Stage 5.
- PPI is not used (`features.ppi: false`), so no interaction-database release needs
  to be pinned.

## Proposed next step

**Stage 5: validation.**
1. Temporal back-test at several cutoff years:
   - train on genes discovered before the cutoff (features, labels and tuning all
     pre-cutoff);
   - score genes discovered at or after the cutoff against non-candidate negatives.
2. Baselines at every cutoff:
   - GTEx nerve expression;
   - GTEx mean expression;
   - −LOEUF;
   - non-atlas model.
3. Expression-matched evaluation and per-category results.
4. Interval simulation: rank of each newly discovered gene at the time it was found.
5. Failure analysis: well-known genes the model ranks low, and why.
6. Activate the end-to-end leakage test before reporting any back-test number.

Decisions needed at the start of Stage 5:
- **Cutoff years.** Proposal: every 2–3 years from 2005, with at least 30 pre-cutoff
  tier 1 genes. To be confirmed after counting genes per year.
- **Jung via the 240 MB object.** Whether to try it before Stage 5, or leave Jung
  out.
