# Stage 4 — Full feature set and model

> Note (2026-11-03): this report describes results before the inner cross-validation correction. Current values are in `results/` and `manuscript/revision_v3_notes.md` (section 6); see `docs/decisions.md`.

**Date:** 2026-10-07 (updated the same day after adding Jung 2023; 11 atlases)

**Outputs:**
- `results/tables/stage4/`:
  - `cv_performance.tsv`: AUROC and average precision with 95% bootstrap CIs.
  - `paired_auroc_differences.tsv`: paired bootstrap differences.
  - `oof_*.tsv.gz`: out-of-fold scores per model and positive set.
  - `cv_best_params.json`: inner-CV hyperparameters per outer fold.
  - `logistic_full_coefficients.tsv`.
  - `ranking_current_logistic_full.tsv.gz`: all 19,229 genes.
  - `summary.json`.
  - `yadav_motoneuron_qc.tsv`.
- Scripts:
  - `workflow/scripts/standardize.py`: per-atlas readers.
  - `workflow/scripts/seurat_to_mtx.R`.
  - `workflow/scripts/stage4_model.py`.
- Package: `features.py`, `model.py` and `coexpression.py` in `src/pnprank`.
- Logs (not committed):
  - reported run: `logs/stage4_model_run4.log` (11 atlases);
  - earlier run: `logs/stage4_model_run3.log` (10 atlases).

**This is not the temporal back-test.** It uses all currently known genes and asks
one question: do single-cell atlas features add signal beyond simple non-atlas
predictors? The back-test (Stage 5) is the test of the central claim.

## What was done

1. **Standardised 11 atlases.** Each atlas was reduced to raw counts, author labels,
   Ensembl IDs and pseudobulk (cell type × donor). Mouse genes were mapped to human
   by Ensembl Compara 116 one-to-one orthologues.

   | Atlas | Species | Tissue | Cells* | Donors/replicates | Harmonised types used |
   |---|---|---|---|---|---|
   | yang2022_human_tg | human | trigeminal ganglion | 38,028 | 3 | 12 |
   | heming2025_human_sural | human | sural nerve (controls) | 33,732 | 4 | 10 |
   | jung2023_human_drg | human | dorsal root ganglion | 111,925 | 7 (18 libraries) | 7 |
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

3. **Built features for 19,229 protein-coding genes.** These are the genes seen in at
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
   - Total: 239 static features, 250 with label-derived features.

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
   - Negatives: protein-coding genes that are not candidates (18,982).
   - Other candidates are dropped from both sets.
   - AUROC and average precision have 1,000-resample stratified bootstrap 95% CIs.
     Model comparisons use paired bootstrap on identical resamples.

6. **Fitted a final logistic model** on all tier 1 genes. It gives the current ranking
   and the coefficients.

7. **Tests:** 24 passed and 1 skipped (the Stage 5 end-to-end back-test test).

## What was found

**Atlas features add signal, but for the primary comparison the gain is borderline**

Performance (out-of-fold, 11 atlases):

| Positives | Score | AUROC [95% CI] | Avg. precision [95% CI] |
|---|---|---|---|
| tier 1 | GTEx mean expression | 0.672 [0.626, 0.723] | 0.010 [0.008, 0.013] |
| tier 1 | GTEx tibial nerve | 0.705 [0.655, 0.756] | 0.013 [0.010, 0.023] |
| tier 1 | −LOEUF | 0.649 [0.597, 0.699] | 0.008 [0.007, 0.011] |
| tier 1 | logistic, non-atlas | 0.744 [0.699, 0.787] | 0.021 [0.014, 0.039] |
| tier 1 | **logistic, full** | **0.792 [0.750, 0.833]** | **0.065 [0.031, 0.127]** |
| tier 1 | HGB, non-atlas | 0.694 [0.643, 0.742] | 0.010 [0.008, 0.014] |
| tier 1 | HGB, full | 0.756 [0.708, 0.804] | 0.020 [0.014, 0.033] |
| tier 1+2 | GTEx tibial nerve | 0.671 [0.642, 0.699] | 0.024 [0.020, 0.036] |
| tier 1+2 | logistic, non-atlas | 0.705 [0.674, 0.733] | 0.031 [0.025, 0.046] |
| tier 1+2 | **logistic, full** | **0.764 [0.735, 0.791]** | **0.069 [0.046, 0.099]** |
| tier 1+2 | HGB, non-atlas | 0.708 [0.674, 0.738] | 0.030 [0.025, 0.039] |
| tier 1+2 | HGB, full | 0.745 [0.715, 0.775] | 0.043 [0.033, 0.061] |

Prevalence is 0.0052 (tier 1) and 0.0125 (tier 1+2).

Paired AUROC differences (11 atlases; the 10-atlas run is shown for comparison):

| Positives | Comparison | ΔAUROC [95% CI], 11 atlases | Fraction ≤ 0 | 10 atlases |
|---|---|---|---|---|
| tier 1 | logistic full − non-atlas | +0.048 [−0.000, 0.091] | 0.026 | +0.052 [0.011, 0.094] |
| tier 1 | HGB full − non-atlas | +0.063 [0.009, 0.111] | 0.010 | +0.064 [0.017, 0.116] |
| tier 1 | logistic full − HGB full | +0.036 [0.000, 0.074] | 0.025 | +0.041 [0.005, 0.076] |
| tier 1+2 | logistic full − non-atlas | +0.059 [0.030, 0.089] | 0.000 | +0.054 [0.028, 0.084] |
| tier 1+2 | HGB full − non-atlas | +0.037 [0.002, 0.073] | 0.018 | +0.044 [0.010, 0.078] |
| tier 1+2 | logistic full − HGB full | +0.019 [−0.003, 0.041] | 0.044 | +0.028 [0.005, 0.052] |

- **The direction is consistent; the size is modest, about +0.05 AUROC.**
  - All six differences are positive in both runs.
  - For the primary comparison (tier 1, logistic), the 11-atlas 95% CI touches zero
    (lower bound rounds to 0.000; 2.6% of resamples ≤ 0).
  - The tier 1+2 gain is clearer: +0.059 [0.030, 0.089].
  - A fair summary: atlas features probably add a small amount of discrimination.
    The evidence is not strong for tier 1 alone.
- **Cross-validation split noise is of the same order as the CI margin.**
  - Adding Jung added one gene to the universe (19,228 → 19,229). That changed the
    stratified fold assignment.
  - The non-atlas logistic model, whose features did not change, moved from 0.739 to
    0.744.
  - The bootstrap CIs above do not include this split-to-split variance. Stage 5
    should repeat CV over several seeds and report it.
- **With all features, logistic beats gradient boosting** in both positive sets, though
  the tier 1+2 difference is not significant. With 100 positives a strongly regularised linear
  model is more stable. The inner CV chose C = 0.001 for the final model.
- **Top-k enrichment (tier 1, out-of-fold, 11 atlases):**
  - top 100 genes (0.5% of the universe): 11 of 100 positives, about 20× the 0.5
    expected by chance;
  - top 1,000: 34 of 100.
  - Average precision is low in absolute terms (0.065). Most genes ranked high are
    not known neuropathy genes.
- **More atlases did not improve overall discrimination.**

  | Atlases | logistic full, tier 1 |
  |---|---|
  | 6 | 0.793 [0.755, 0.830] |
  | 10 | 0.791 [0.750, 0.829] |
  | 11 (with Jung) | 0.792 [0.750, 0.833] |

  The added atlases are kept because category-specific analyses need them: motor
  neurons for dHMN, DRG and nerve cell types for CMT and HSAN.

**Largest standardised coefficients (final logistic model) are biologically
sensible**

- Positive:
  - gnomAD missense Z (+0.35);
  - motor-neuron specificity in human spinal cord (+0.19);
  - myelinating-Schwann specificity in mouse TG (+0.19);
  - C-LTMR specificity in mouse DRG (+0.18);
  - SST sensory-neuron expression in human TG (+0.17);
  - DRG-neuron specificity in human DRG (Jung, +0.17);
  - Schwann specificity at the mouse NMJ (+0.15).
- Negative: specificity for fibroblasts, OPCs, immune cells, endothelium (Jung) and
  cholinergic interneurons.
- **Caveat:** these are coefficients of correlated features under strong
  regularisation. They show what the model uses, not independent effects.

**The current ranking still contains an expression-level confound**

- **Top 30:** MPZ, TUBA1A, PRX, PMP22, L1CAM, NEFH, DYNC1H1, TUBA1B, NEFL, IL16, MBP,
  SLC36A2, FA2H, TARDBP, TUBB, BCAS1, TUBB2A, ACTG1, PRNP, CADM3, KCNA2, ATP1A1,
  YWHAG, HSPA12A, UGT8, ANXA2, KIF1A, TUBB4A, MAP1B, KCNA1.
- **Top 100 composition:**
  - 12 tier 1 genes;
  - 8 tier 2 genes;
  - 1 unverified candidate (L1CAM).
- **Housekeeping and cytoskeletal genes rank high:**
  - TUBA1A, TUBA1B, TUBB and TUBB2A;
  - ACTG1;
  - YWHAG.
- **Other top-ranked negatives.** Many are myelin or neuronal genes, for example MBP,
  UGT8, BCAS1 and KCNA1. They are plausible by expression pattern, which is what the
  model can see. Nothing here is evidence that they cause neuropathy.
- **Ranking scope.** This ranking is trained on current labels and is shown for
  inspection only. It is not a validated prediction.

## Surprises and problems

- **The atlas gain is fragile for tier 1.** It depends on CV split noise of about
  0.005 AUROC, and its CI lower bound sits near zero. Stage 5 must report:
  - repeated-CV variance;
  - the back-test, rather than relying on one CV split.
- **Expression level is still a confound.** It was found in the pilot and it persists.
  The model has no feature that separates "highly expressed in nerve/neuron" from
  "broadly highly expressed".
  - Stage 5 should add an expression-matched evaluation and report it next to the raw
    AUROC.
  - A candidate model change is residualising atlas expression on GTEx mean
    expression. It would be pre-specified before the back-test, not tuned on it.
- **More atlases did not help overall discrimination.** Gains, if any, may show up
  only within categories. Stage 5 should report per-category results (CMT, dHMN,
  HSAN, HSP).
- **Jung 2023 needed less memory than feared.**
  - The 9.6 GB serialised object loaded with about 3 GB resident memory. Extraction
    took under 4 minutes.
  - The deposited human object has a single "DRG Neuron" class, with no neuron
    subtypes. Jung therefore adds human DRG non-neuronal types and one pooled neuron
    profile.
  - Donors were taken from the GEO series matrix (18 libraries, 7 donors).
- **Yim 2022 mouse nerves excluded.** Only numeric clusters were deposited.
- **Yadav donors 1 and 5 fail the motoneuron QC.** Donor 5's "motoneurons" look like
  glial contamination.
- **Memory pressure.** One earlier model run thrashed under memory pressure and was
  killed. The reported runs completed normally.

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
3. Repeated CV (several seeds) for the current-label comparison, to quantify split
   variance.
4. Expression-matched evaluation and per-category results.
5. Interval simulation: rank of each newly discovered gene at the time it was found.
6. Failure analysis: well-known genes the model ranks low, and why.
7. Activate the end-to-end leakage test before reporting any back-test number.

Decision needed at the start of Stage 5:
- **Cutoff years.** Proposal: every 2–3 years from 2005, with at least 30 pre-cutoff
  tier 1 genes. To be confirmed after counting genes per year.
