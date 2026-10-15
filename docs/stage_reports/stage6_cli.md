# Stage 6 — Command-line tool

**Date:** 2026-10-07

**Outputs:**
- `results/tables/pnprank_ranking.tsv.gz`: one row per protein-coding gene (19,229), about
  0.8 MB, committed so lookups need no raw data.
- `results/tables/pnprank_ranking_meta.json`: label counts, seeds, the C chosen by every
  fold model, and seed-stability statistics.
- `results/tables/stage6/release_checks.tsv`: post hoc checks of the ranking.
- Code:
  - `src/pnprank/ranking.py`: `build_ranking`, `resolve`, `top_genes`;
  - `src/pnprank/cli.py`: `build-ranking`, `query`, `top`;
  - `src/pnprank/backtest.py`: `fit_fold_models`, factored out of `oof_scores_at_cutoff`;
  - `workflow/rules/tool.smk` and `workflow/scripts/stage6_checks.py`.

**Tests:** 35 pass, including the 10 leakage tests. Six new tests in `tests/test_cli.py`
cover:
- gene lookup: Ensembl ID with version, case-insensitive symbol, alias, an ambiguous
  alias, the approved symbol taking precedence over an alias;
- `top` filtering and ordering, and the CLI exit codes;
- the ranking build: statuses, every non-tier-1 gene scored with a unique rank, the
  percentile inside its seed range, reproducibility, and a planted positive-like gene
  ranking above the 90th percentile.

## What the tool does

```bash
pnprank query SPTLC1 KIF5A ENSG00000158887   # symbols, aliases or Ensembl IDs
pnprank query --file genes.txt --tsv         # one gene per line; '-' reads stdin
pnprank top -n 50                            # best-ranked genes not yet known
pnprank top -n 50 --expr-matched             # order by expression-matched percentile
pnprank build-ranking --n-jobs 5             # rebuild from the processed data
```

Every lookup prints a note on stderr: percentiles prioritise genes for follow-up; they are
not evidence of causality. Unmatched and ambiguous queries are reported on stderr. The
exit code is 1 only when nothing matched.

## How the ranking is built

The ranking applies the procedure validated in the Stage 5 back-test to today's labels; the
release itself is not back-tested.
- **Training labels:** the 100 tier-1 genes are positives. The 140 tier-2 genes are left
  out of training. The other 18,989 genes are labelled 0.
- **Scoring:** the 0-labelled genes are split into 5 folds. Each fold is scored by a model
  trained on the positives and the other four folds, with C tuned by inner CV.
- **Percentile:** each gene's percentile is taken within the fold model that scored it.
- **Tier-2 genes:** they are scored by all five fold models and given the mean of their
  percentiles.
- **Tier-1 genes:** they were training positives, so they have no percentile or rank. The
  tool reports them as known.
- **Shared code:** the release calls the same `fit_fold_models` as the back-test. After the
  refactor, the cutoff-2005 back-test scores were reproduced to 1e-16.

The two Stage 6 decisions you agreed to (`docs/decisions.md`):
- **Regularisation:** C is tuned by inner CV, as in the back-test. The reported percentile
  comes from one seed. Nine more seeds give each gene's percentile range, plus the C
  chosen by every fold model.
- **Columns shown beside the percentile:** the expression-matched percentile (within
  GTEx mean-expression decile) and the non-atlas model's percentile.

The Stage 4 ranking (`ranking_current_logistic_full.tsv.gz`) is superseded. It scored
genes in-sample with a model that had seen them as negatives.

## What was found

**Stability across 10 seeds**
- Every fold model chose C = 0.001 for the full model.
- Spearman correlation of percentiles with the primary seed: 0.999.
- Each gene's percentile range across seeds: median 0.024 overall, 0.002 in the top 100.
- The top 100 shares 80 to 91 genes with the primary seed's top 100.

With 100 positives, C selection is stable, unlike the 30-40-positive back-test cutoffs.

**C sits at the bottom of the grid.** Because C = 0.001 is the smallest grid value, I
tested a wider grid on the release training set.
- The optimum was C = 3e-4.
- Its inner-CV AUROC beat C = 0.001 by only 0.001-0.006 (three seeds).
- The grid was kept as back-tested.

**Held-out tier-2 genes (post hoc, not pre-specified).** Tier-2 genes never enter
training, so their rank against unlabelled genes is a held-out check:

| Score | AUROC (95% CI), 140 tier-2 vs 18,982 unlabelled |
|---|---|
| percentile (full model) | 0.713 (0.674-0.752) |
| percentile, expression-matched | 0.660 (0.614-0.704) |
| percentile, non-atlas model | 0.647 (0.605-0.689) |
| full − non-atlas (paired) | +0.067 (+0.036 to +0.101) |

- 42 of 140 tier-2 genes are in the top 10%.
- Atlas features help here as they did in Stage 4 CV and in the tier 1+2 back-test.
- This is current knowledge, not a temporal test.
- Tier 2 mixes established HSP genes with secondary or atypical neuropathy genes, so it is
  a weaker and noisier positive set than tier 1.

**Top of the list**

| | Top 100 | Top 500 |
|---|---|---|
| Unlabelled | 91 | 478 |
| Tier 2 | 8 | 21 |
| Listed but not in the truth set | 1 (L1CAM) | 1 |
| In the top GTEx expression decile | 40% | 31% |
| In deciles 1-5 | 28% | 26% |

- **Highly expressed genes.** Some of the top unlabelled genes are highly expressed
  cytoskeletal or housekeeping genes: TUBA1A, TUBA1B, TUBB, ACTG1, CLTC. Within the top
  decile, expression matching does not demote them. This confound is in
  `docs/limitations.md`.
- **Neural genes.** Others are glial or neuronal genes whose percentile rises with atlas
  features: KCNA1, MLIP, SCN7A, SLC36A2, MTCL3. Their non-atlas percentile is 0.49-0.69.
  These are the kind of candidates the atlas features are meant to surface.
- **L1CAM** ranks 2nd. It is on a source gene list but was not curated into the truth
  set, so it was a training 0. Its high rank is a consistency check, not a discovery.
- **No claim of novelty.** I have not checked whether any top unlabelled gene already has
  a published neuropathy link. That needs a literature search with a recorded source for
  each gene (hard rule 1).

## Surprises

- `query` printed ranks as `2.0`: the row-by-row lookup lost the integer dtype. Fixed and
  covered by a test.
- `docs/decisions.md` still said 19,228 genes from before Jung 2023 was added. The correct
  count is 19,229, and the line now says so.

## Proposed next step

All six stages are done. Possible follow-ups, none started:
1. Compare PNPrank with existing gene-prioritisation tools in the same temporal back-test.
   Reviewers are most likely to ask for this.
2. Curate the top ~50 unlabelled genes against the literature, with a source for each,
   to separate already-reported links from genuinely new candidates.
3. List top candidates that have clear C. elegans orthologues, for the collaborator.
4. Write a two-page summary for collaborators and draft the paper.

## Addendum (2026-10-16): intended use alongside sequencing

- **The ranking does not find variants.** It orders genes for follow-up once sequencing has narrowed a locus.
  - **WGS of a linked family.** An interval can still hold many rare variants and several structural variants, and the ranking gives an order for checking them.
  - **Structural or regulatory variants.** Where the cause acts through one of these, the ranking can miss it (CMTX3), or it only helps after the variant has pointed to the gene (DHMN1, UBE3C 11 of 146).
- **Atlases add what patient RNA-seq on blood or fibroblasts cannot see.** In GTEx v10:
  - 33 of 100 tier-1 genes are below 1 TPM in whole blood;
  - 16 are below 1 TPM in cultured fibroblasts;
  - 15 are below 1 TPM in both, including GJB1, SH3TC2, PMP2 and DRP2, which are well expressed in tibial nerve.

  The atlases show which nerve, ganglion or spinal cord cell types express such genes (`docs/case_studies/patient_tissue_expression.md`).
- **Best setting:** loci that remain unsolved after sequencing.
- **Scoring a user-supplied interval.** Use GRCh38 coordinates, lifting over from hg19 or T2T-CHM13 and checking against flanking markers. Unpublished intervals are filtered outside the repository. A prospective test fixes the ranking before the intervals are seen (`docs/decisions.md`, 2026-10-16).
