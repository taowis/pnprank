# Stage 5 — Validation: temporal back-test, baselines, interval simulation, failure analysis

**Date:** 2026-10-07

**Outputs (`results/backtest/`):**
- `backtest_performance.tsv`: AUROC (95% CI) and average precision for every:
  - cutoff;
  - positive set;
  - evaluation window (`all_after`, `within_5y`);
  - score transform;
  - score.
- `backtest_paired_differences.tsv`: the full model against each comparator (paired bootstrap).
- `interval_simulation.tsv`, `interval_summary.tsv`, `interval_paired.tsv`.
- `failure_analysis.tsv`, `failure_feature_contrast.tsv`.
- `repeated_cv_tier1.tsv`.
- Figures: `results/figures/backtest_auroc_tier1.png`, `backtest_auroc_tier1_2.png` and
  `interval_simulation_percentiles.png`.
- Scripts: `workflow/scripts/stage5_backtest.py`, `stage5_summarise.py` and
  `stage5_repeated_cv.py`.
- Functions in `src/pnprank/backtest.py`: `training_labels`, `oof_scores_at_cutoff`,
  `within_fold_percentile` and `evaluation_labels`.

**Leakage tests:** `tests/test_leakage.py` passes in full (10 leakage tests; 29 tests in the whole suite at the time), including the
end-to-end back-test test, which is no longer skipped. The end-to-end test checks three
things:
- Removing, moving or adding post-cutoff truth genes, or changing their tier, leaves
  every training label and every score unchanged.
- Pre-cutoff positives are never scored.
- Scores are identical under three Python hash seeds.

## Design (pre-specified, with one recorded deviation)

The design was committed before any back-test was run (`docs/decisions.md`, commit
e9740cf).

**Cutoffs and labels**
- Cutoffs: 2005, 2008, 2011, 2014, 2017 and 2020.
- Training labels at cutoff c:
  - 1 = tier-1 genes discovered before c.
  - Tier-2 genes known before c are dropped.
  - 0 = every other gene, including future positives and today's unverified
    candidates.
- Label-derived co-expression features use pre-cutoff genes only. C is tuned by inner
  CV on training genes only.

**Out-of-fold scoring**
- Genes labelled 0 are split into 5 folds.
- Each fold is scored by a model that did not see it.

**Evaluation**
- Positives: tier-1 genes discovered at or after c.
- Negatives: non-candidate protein-coding genes.
- The primary model is logistic regression on all features. Comparators:
  - the non-atlas model;
  - GTEx tibial nerve expression;
  - GTEx mean expression;
  - −LOEUF.

**Deviation (found while running).**
- *What went wrong:* the 5 fold models choose their own C, so their probabilities are
  on different scales, and pooling raw scores across folds mixes those scales. At
  cutoff 2008, one fold chose C = 1.0. The median score per fold ranged from 0.000 to
  0.114. The pooled AUROC (0.675) was below four of the five within-fold AUROCs, which
  ranged from 0.675 to 0.854.
- *Primary analysis now:* each score is converted to a percentile within the fold model
  that produced it. This is applied to every score, to the expression-matched analysis
  and to the interval simulation.
- *Pre-specified analysis:* the pooled analysis is still reported (`transform = pooled`).

**Reproducibility bug (fixed).**
- *What went wrong:* training rows were assembled from a Python set, so inner-CV folds,
  and with them the chosen C, changed between runs with the same seed.
- *Fix:* rows now follow label order. A test checks that scores are identical under
  three hash seeds, and re-scoring cutoff 2008 in a fresh process gave identical C and
  scores (maximum difference 1e-16).
- All results below come from the fixed code.

## What was found

### 1. Temporal back-test: PNPrank is the best or tied-best score at every cutoff

> Correction (2026-11-02): at 2020 the non-atlas model is 0.002 higher (0.755 vs 0.753), and 0.013 higher after expression matching (0.681 vs 0.668). "Best or tied-best" holds against the dated existing approaches, not against this ablation.

Tier 1, all genes discovered after the cutoff, within-fold percentiles (primary):

| Cutoff | Train / test genes | PNPrank | Non-atlas model | GTEx nerve | GTEx mean | −LOEUF |
|---|---|---|---|---|---|---|
| 2005 | 30 / 70 | **0.759** [0.710, 0.808] | 0.732 | 0.714 | 0.674 | 0.662 |
| 2008 | 38 / 62 | **0.773** [0.720, 0.828] | 0.746 | 0.714 | 0.678 | 0.650 |
| 2011 | 46 / 54 | **0.766** [0.702, 0.826] | 0.751 | 0.720 | 0.686 | 0.647 |
| 2014 | 68 / 32 | **0.781** [0.696, 0.856] | 0.710 | 0.726 | 0.690 | 0.578 |
| 2017 | 81 / 19 | **0.762** [0.650, 0.858] | 0.737 | 0.755 | 0.735 | 0.657 |
| 2020 | 87 / 13 | 0.753 [0.635, 0.866] | 0.755 | 0.748 | 0.712 | 0.618 |

Paired differences, PNPrank minus comparator. Each cell is ΔAUROC [95% CI].

| Cutoff | vs non-atlas | vs GTEx nerve | vs GTEx nerve, expression-matched |
|---|---|---|---|
| 2005 | +0.027 [−0.028, 0.084] | +0.045 [−0.016, 0.110] | **+0.092 [0.010, 0.179]** |
| 2008 | +0.027 [−0.033, 0.087] | +0.059 [−0.004, 0.119] | **+0.114 [0.038, 0.188]** |
| 2011 | +0.014 [−0.053, 0.077] | +0.046 [−0.026, 0.124] | **+0.117 [0.018, 0.219]** |
| 2014 | +0.071 [−0.002, 0.144] | +0.055 [−0.042, 0.152] | +0.098 [−0.032, 0.222] |
| 2017 | +0.024 [−0.039, 0.094] | +0.007 [−0.100, 0.113] | +0.105 [−0.064, 0.270] |
| 2020 | −0.003 [−0.081, 0.073] | +0.004 [−0.130, 0.124] | +0.080 [−0.137, 0.291] |

**Tier 1 summary**
- **Better than every simple baseline.** PNPrank scores genes discovered after the
  cutoff above every simple baseline at all six cutoffs, with AUROC 0.75–0.78. The CIs
  of single cutoffs are wide, and for raw AUROC none of the differences from GTEx
  tibial nerve expression excludes zero.
- **Advantage grows after expression matching.** With expression matched:
  - GTEx nerve drops to 0.58–0.63;
  - PNPrank stays at 0.67–0.73;
  - the difference excludes zero at 2005, 2008 and 2011.

  The atlas-based model ranks future genes better than "is it highly expressed in
  nerve" once the generic expression-level signal is removed.
- **Atlas gain over the non-atlas model is small and not individually significant.**
  - The difference is positive at 5 of 6 cutoffs, about +0.01 to +0.07.
  - The CI includes zero at every cutoff; the closest is 2014, +0.071 [−0.002, 0.144].
  - The cutoffs share test genes, so they are not independent replicates. They cannot
    be combined into a single significance test.

**Tier 1+2 sensitivity analysis (24–164 test genes)**
- **Stronger and consistent.** PNPrank reaches AUROC 0.73–0.81.
- **Against the non-atlas model:**
  - positive at all six cutoffs, +0.040 to +0.062;
  - the CI excludes zero at 2005, 2011 and 2014.
- **Against GTEx nerve:**
  - +0.065 to +0.103, with the CI excluding zero at every cutoff;
  - after expression matching, +0.10 to +0.17.

**Other windows and transforms**
- **5-year window (tier 1, genes discovered in [c, c+5)).** This uses 10–33 test genes.
  PNPrank scores 0.75–0.82, against 0.70–0.78 for GTEx nerve and 0.71–0.76 for the
  non-atlas model.
- **Pre-specified pooled analysis.** PNPrank tier-1 AUROC is 0.737, 0.732, 0.767, 0.780,
  0.762 and 0.755 across the six cutoffs, close to the primary analysis.

### 2. Interval simulation: each gene ranked in the year it was discovered

70 tier-1 genes discovered from 2005 to 2026. Each was ranked by a model trained on the
genes known before its discovery year, among about 19,000 scored genes (within-fold
percentile; 1 = top).

| Score | Median percentile [95% CI] | Top 1% | Top 5% | Top 10% |
|---|---|---|---|---|
| PNPrank | 0.830 [0.758, 0.897] | 8 | 16 | 26 (37%) |
| Non-atlas model | 0.833 [0.710, 0.870] | 5 | 17 | 24 (34%) |
| GTEx tibial nerve | 0.811 [0.714, 0.849] | 0 | 11 | 19 (27%) |
| GTEx mean expression | 0.740 [0.626, 0.801] | 1 | 4 | 13 (19%) |
| −LOEUF (n = 64) | 0.724 [0.608, 0.794] | 1 | 7 | 13 (20%) |

Paired, same genes (PNPrank minus comparator):

| Comparator | Mean percentile difference [95% CI] | Δ genes in top 10% [95% CI] |
|---|---|---|
| Non-atlas model | +0.033 [−0.018, 0.080] | +2 [−6, 10] |
| GTEx tibial nerve | +0.056 [−0.003, 0.114] | +7 [−2, 16] |
| GTEx mean expression | +0.096 [0.035, 0.164] | +13 [3, 23] |
| −LOEUF | +0.110 [0.038, 0.184] | +11 [2, 20] |

- **Top of the ranking.** 8 of the 70 genes were in the top 1% in their discovery year
  (about 190 genes). The tibial-nerve baseline put none there, and the non-atlas model
  put 5.
- **Medians.** The median percentile is similar for PNPrank, the non-atlas model and
  GTEx nerve. Most of the gain is concentrated at the top of the ranking.
- **By category** (small groups, wide CIs):

  | Category | Genes | PNPrank median | GTEx nerve median | Non-atlas median |
  |---|---|---|---|---|
  | CMT | 41 | 0.855 | 0.822 | 0.848 |
  | HSAN | 10 | 0.867 | 0.684 | 0.812 |
  | dHMN | 14 | 0.755 | 0.820 | 0.748 |

  - HSAN is where PNPrank departs most from GTEx nerve (5 of 10 HSAN genes in the top
    10%, against 0 for GTEx nerve).
  - dHMN is the weakest category.

### 3. Repeated CV on current labels (Stage 4 comparison, 10 seeds)

| | Mean (SD) | Range |
|---|---|---|
| Full model AUROC | 0.792 (0.009) | 0.782–0.810 |
| Non-atlas model AUROC | 0.736 (0.008) | 0.723–0.749 |
| Difference | **+0.057 (0.008)** | +0.048 to +0.073 |

- The atlas gain on current labels is positive under every split, at least +0.048.
- The borderline Stage 4 CI reflected the bootstrap uncertainty of a single split.
  Split-to-split SD of AUROC is about 0.009.

### 4. Failure analysis

Eight of the 70 genes ranked in the bottom half in their discovery year
(`failure_analysis.tsv`):

| Gene | Year | Category | Percentile |
|---|---|---|---|
| FBLN5 | 2011 | CMT | 0.20 |
| CFAP276 | 2019 | CMT | 0.22 |
| SLC5A7 | 2012 | dHMN | 0.23 |
| PRPS1 | 2007 | CMT | 0.27 |
| VWA1 | 2021 | dHMN | 0.34 |
| RCC1 | 2025 | CMT | 0.38 |
| MORC2 | 2015 | CMT | 0.41 |
| SYT2 | 2014 | dHMN | 0.47 |

Four patterns appear among them:

- **Cell-type-specific but lowly expressed motor-neuron genes.** SLC5A7 and SYT2 are
  most specific to spinal motor neurons, which is correct biology. But SLC5A7 is almost
  absent from bulk nerve (GTEx log TPM 0.02), and the model weighs overall expression
  heavily. dHMN genes from rare cell types are penalised.
- **Missing mouse features.** PRPS1 and MORC2 have no one-to-one mouse orthologue in
  Ensembl 116, so all five mouse atlases are missing for them. The missing-value
  handling (median imputation plus indicators) does not recover this.
- **Genes acting outside neurons and glia.**
  - FBLN5 is an extracellular-matrix gene, high in nerve but most specific to
    fibroblast-like cells.
  - VWA1 is highly expressed in nerve but tolerant to loss of function (LOEUF 1.51).
- **Weak prior evidence.** CFAP276 (cilia-related, LOEUF 1.49, low nerve expression)
  has no expression or constraint signal pointing to neuropathy.

Genes ranked in the top 10% had higher nerve expression (median log TPM 4.0 vs 2.8 for
failures) and similar constraint (LOEUF 0.55 vs 0.58).

## Surprises and problems

- **Pooling OOF scores across separately tuned fold models is invalid.** It produced a
  spurious drop at 2008. The fix is principled and was applied to all scores, but it is
  a deviation from the pre-specified design. Both analyses are reported, and they agree
  closely after the reproducibility fix.
- **C selection is unstable with 30–40 positives.**
  - It changed with nothing more than row order. That is how the hash-order bug showed
    up.
  - Within-fold percentiles remove the scale problem but not this variance.
  - A fixed C, or a more stable tuning criterion, is worth considering in Stage 6 for
    the released model. It would not be tuned on these results.
- **The atlas gain over the non-atlas model is consistent in direction but not
  significant at any single tier-1 cutoff.**
  - The clearer evidence comes from:
    - tier 1+2;
    - the repeated CV on current labels;
    - the comparison with GTEx nerve after expression matching.
  - The honest summary: atlas features help modestly, mostly by moving a handful of
    genes into the top 1%, and by separating nerve-relevant from generically
    highly-expressed genes.
- **Expression level remains the dominant confound.** Every score loses AUROC after
  expression matching. In tier 1, GTEx nerve loses most (−0.10 to −0.17) and PNPrank
  least (−0.04 to −0.09).
- **Later cutoffs have few test genes.** 2017 has 19 and 2020 has 13, so their CIs span
  about 0.2. They cannot distinguish the scores.

## Is the central claim supported?

- **"A model trained only on genes known before a cutoff ranks later genes near the
  top":**
  - **Partly supported.** Future genes sit at a median percentile of about 0.83, and
    37% fall in the top 10% at discovery.
  - **Weaker at the very top.** Tier-1 average precision stays low (0.005–0.016), so most
    top-ranked genes are not (yet) neuropathy genes.
- **"…and beats simple baselines":**
  - **Supported against GTEx mean expression and −LOEUF.**
  - **Supported against GTEx nerve expression after expression matching.** Without
    matching, the advantage is consistent in direction but not significant at single
    cutoffs.
  - **Not significantly better than the non-atlas model** at single tier-1 cutoffs,
    although consistent in direction. It is significant in tier 1+2.

## Proposed next step

**Stage 6: command-line tool.** Package the model as `pnprank rank` and `pnprank explain
<gene>`:
- `pnprank rank`: current-label model, gene ranking with percentile and feature
  contributions.
- `pnprank explain <gene>`: the atlases and cell types driving a gene's score.

The tool also needs:
- documentation of the back-test results and the limitations above;
- tests for the CLI.

Decisions for Stage 6:
- **Regularisation in the released model.** Proposal: keep inner-CV tuning, but report
  the chosen C and its stability.
- **Which outputs to expose.** Proposal: the expression-matched percentile alongside the
  raw score, so users can see when a high rank is driven by expression level.
