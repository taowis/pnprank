# pnprank

**Status:** Stage 6 (Command-line tool) — complete. All stages done. On 2026-11-03 the inner cross-validation was corrected so that label-derived features are recomputed inside each inner split, and every analysis was rerun (`docs/decisions.md`); the released ranking did not change.

| Stage | Description | Status |
|---|---|---|
| 0 | Scaffold | done |
| 1 | Truth set with discovery years | done (242 genes) |
| 2 | Locate and assess atlases | done (13 shortlisted of 133) |
| 3 | Feasibility pilot on one atlas | done (Yang 2022 human TG) |
| 4 | Full feature set and model | done (11 atlases; nested-CV AUROC 0.79 tier 1) |
| 5 | Validation: temporal back-test, baselines, interval simulation, failure analysis | done (back-test AUROC 0.75–0.79 tier 1 after the 2026-11-03 rerun; see report) |
| 6 | Command-line tool | done (`pnprank query`, `pnprank top`, `pnprank interval`; see report) |

## What this is

PNPrank ranks human genes by their likelihood of causing inherited peripheral
neuropathy (CMT, dHMN, HSAN and related HSP) using cell-type-resolved
expression from public single-cell atlases.

The central claim to support or refute: a model trained only on genes known
before a cutoff year ranks genes discovered after that year near the top, and
beats simple baselines. Negative results will be reported as such.

## Setup

```bash
conda env create -f environment.yml
conda activate pnprank
pip install -e . --no-deps
cp config/paths.example.yaml config/paths.yaml   # then edit for this machine
pytest
```

## Using the ranking

The ranking of every protein-coding gene ships in
`results/tables/pnprank_ranking.tsv.gz`, so lookups need no downloads:

```bash
pnprank query SPTLC1 KIF5A ENSG00000158887        # symbols (aliases ok) or Ensembl IDs
pnprank query --file my_candidates.txt --tsv > hits.tsv   # one gene per line
pnprank top -n 50                                  # best-ranked genes not yet known
pnprank top -n 50 --expr-matched                   # order by expression-matched percentile
pnprank interval DHMN1=chr7:138,500,001-155,200,000  # 7q34-q36.2 (UCSC hg38 bands), ranked
pnprank interval --bed loci.bed --pad-kb 500 --tsv > loci.tsv  # several intervals, padded
```

`pnprank interval` takes GRCh38 coordinates, 1-based and inclusive as shown in Ensembl and
the UCSC browser (BED files are read as 0-based, as BED defines). Gene coordinates are
from Ensembl 116 (`results/tables/gene_coordinates_e116.tsv.gz`). A gene is included if any
part of it overlaps the interval. `interval_rank` orders the scored genes inside the
interval; known tier-1 genes have no score and are listed last. Intervals on hg19 or
T2T-CHM13 must be lifted over first, and checked against their flanking markers. Eight ranked genes
have no Ensembl 116 primary-assembly coordinates and never appear in interval output
(listed in `docs/decisions.md`).

Columns:

| Column | Meaning |
|---|---|
| `percentile` | rank among scored genes (1 = top), from a model that never saw the gene as a training example |
| `percentile_seed_min`, `percentile_seed_max` | range over 10 random fold splits; a wide range means an unstable rank |
| `percentile_expr_matched` | rank among genes of similar overall expression (GTEx decile); low here but high in `percentile` means the rank is driven mostly by expression level |
| `percentile_non_atlas` | the same model without single-cell atlas features (GTEx and gnomAD only) |
| `status` | `known_tier1` (training positive, not ranked), `known_tier2`, `listed_not_in_truth` (on a source gene list but not curated as established), `unlabelled` |
| `top_cell_type` | the atlas cell type with the highest specificity for the gene (descriptive) |

What the percentiles mean: in the back-test, tier-1 genes ranked in the year they were
reported had a median percentile of about 0.83, but most top-ranked
genes are not neuropathy genes. Use the ranking to prioritise follow-up, not as
evidence of causality. See `docs/stage_reports/stage5_validation.md` and
`docs/limitations.md`.

To rebuild the ranking (Stages 1–4 must have been run; if `interim/` has been removed, the atlas profiles are rebuilt from `data/processed/features_static.tsv.gz`):
`pnprank build-ranking --n-jobs 5`, or `snakemake -c5 results/tables/pnprank_ranking.tsv.gz`.

## Layout

- `config/` — run settings, atlas registry, machine paths template
- `data/` — not in git except `data/truth/` (curated labels) and `data/README.md`, which lists the source of every other file
- `src/pnprank/` — the importable package and `pnprank` CLI
- `workflow/` — Snakemake rules and thin wrapper scripts; every step is listed in `docs/workflow.md`
- `results/` — small tables and figures (committed); models and matrices are not
- `docs/` — decisions, data sources, limitations, per-stage reports
- `tests/` — including `test_leakage.py`, which must pass before any back-test result is reported

## Data policy

Public data only. No patient or lab-internal data enters this repository.

## Licence

Code: MIT (see `LICENSE`). Data: each source under its own licence, see `docs/data_sources.md`.
