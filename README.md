# pnprank

**Status:** Stage 6 (Command-line tool) — complete. All stages done; awaiting review.

| Stage | Description | Status |
|---|---|---|
| 0 | Scaffold | done |
| 1 | Truth set with discovery years | done (242 genes) |
| 2 | Locate and assess atlases | done (13 shortlisted of 133) |
| 3 | Feasibility pilot on one atlas | done (Yang 2022 human TG) |
| 4 | Full feature set and model | done (11 atlases; nested-CV AUROC 0.79 tier 1) |
| 5 | Validation: temporal back-test, baselines, interval simulation, failure analysis | done (back-test AUROC 0.75–0.78 tier 1; see report) |
| 6 | Command-line tool | done (`pnprank query`, `pnprank top`; see report) |

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
```

Columns:

| Column | Meaning |
|---|---|
| `percentile` | rank among scored genes (1 = top), from a model that never saw the gene as a training example |
| `percentile_seed_min`, `percentile_seed_max` | range over 10 random fold splits; a wide range means an unstable rank |
| `percentile_expr_matched` | rank among genes of similar overall expression (GTEx decile); low here but high in `percentile` means the rank is driven mostly by expression level |
| `percentile_non_atlas` | the same model without single-cell atlas features (GTEx and gnomAD only) |
| `status` | `known_tier1` (training positive, not ranked), `known_tier2`, `listed_not_in_truth` (on a source gene list but not curated as established), `unlabelled` |
| `top_cell_type` | the atlas cell type with the highest specificity for the gene (descriptive) |

What the percentiles mean: in the temporal back-test, genes discovered after a cutoff
had a median percentile of about 0.83 at the time of discovery, but most top-ranked
genes are not neuropathy genes. Use the ranking to prioritise follow-up, not as
evidence of causality. See `docs/stage_reports/stage5_validation.md` and
`docs/limitations.md`.

To rebuild the ranking from the raw data (Stages 1–4 must have been run):
`pnprank build-ranking --n-jobs 5`, or `snakemake -c5 results/tables/pnprank_ranking.tsv.gz`.

## Layout

- `config/` — run settings, atlas registry, machine paths template
- `data/` — not in git except `data/truth/` (curated labels) and `data/README.md`, which lists the source of every other file
- `src/pnprank/` — the importable package and `pnprank` CLI
- `workflow/` — Snakemake rules and thin wrapper scripts; every step is listed in `docs/workflow.md`
- `results/` — small tables and figures (committed); models and matrices are not
- `docs/` — decisions, data sources, limitations, per-stage reports
- `manuscript/` — paper draft and collaborator summary (Word, with an EndNote library)
- `presentations/` — talks, one folder per talk, with the scripts that draw their animations
- `tests/` — including `test_leakage.py`, which must pass before any back-test result is reported

## Data policy

Public data only. No patient or lab-internal data enters this repository.

## Licence

Code: MIT (see `LICENSE`). Data: each source under its own licence, see `docs/data_sources.md`.
