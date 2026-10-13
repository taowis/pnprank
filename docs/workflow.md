# Analysis workflow

This page lists every step of the analysis in the order it runs, with what each step
reads and writes and where its results are discussed. The steps are Snakemake rules
(`Snakefile`, `workflow/rules/*.smk`) that call thin scripts in `workflow/scripts/`;
the analysis code itself is in `src/pnprank/`. A few steps are manual or call the
`pnprank` CLI directly; they are marked below.

Settings live in `config/`: `config.yaml` (cutoff years, seeds, model and bootstrap
settings, reference file URLs and checksums), `atlases.yaml` (atlas registry with file
URLs and SHA-256), `celltype_map.tsv` (author cell-type labels to shared classes),
`comparators.yaml` (dated GO and BioGRID releases), `loci.tsv` (linkage intervals) and
the machine-specific `paths.yaml` (data root, disk budget). Each non-obvious choice is
in `docs/decisions.md`; data provenance is in `docs/data_sources.md` and `data/README.md`.

## Running

```
conda activate pnprank
snakemake -n -c5               # list what would run
snakemake -c5                  # everything except the live web queries
snakemake -c1 live_queries     # literature and C. elegans queries (results change over time)
pytest                         # unit tests, including the leakage tests
```

On the current data root `snakemake -n` reports nothing to do, which means the
committed results match the current code and inputs as far as file times can tell.
Snakemake decides what to rerun from modification times. Reference files that are fixed
by release and checked by checksum are marked `ancient()`, so fetching them again does
not force the atlases to be reprocessed. Starting from an empty data root, the workflow
downloads about 13 GB of atlas files (the largest single file is 4.2 GB) and checks
the disk budget before each download.

`tests/test_leakage.py` is a rule input (`leakage_tests`) for every step that
summarises back-test results, so those steps cannot run unless the leakage tests pass.

## Steps

### Stage 1: truth set

Not rules, because the curation is partly manual and the source snapshots cannot be
fetched again in their original versions.

1. `pnprank build-candidates` snapshots PanelApp panels 85, 846 and 165, the ClinGen
   gene-validity export and the HGNC complete set into `data/external/truth_sources/`
   and writes the candidate list `data/truth/candidates.tsv`.
2. Discovery papers were curated per gene under `docs/truth_protocol.md` and
   `docs/truth_protocol_round2.md`, giving `data/truth/curation.json` and
   `data/truth/review.tsv`.
3. `pnprank assemble-truth` checks every PMID against PubMed, computes the discovery
   year and writes `data/truth/neuropathy_genes.tsv`.
4. Rule `validate_truth` checks the schema.

Report: `docs/stage_reports/stage1_truth_set.md`; build notes in `data/truth/sources.md`.

### Stage 2: atlas search

Manual search under `docs/atlas_search_protocol.md`. Results: `config/atlases.yaml`
(shortlist) and `results/tables/atlas_candidates.tsv` (all 133 candidates).
Report: `docs/stage_reports/stage2_atlases.md`.

### Downloads (`workflow/rules/external.smk`)

| Rule | Script | Writes |
|---|---|---|
| `fetch_external` | `fetch_external.py` | gnomAD v4.1 constraint, GTEx v10 median TPM, Ensembl 116 mouse GTF and gene table, Compara 116 mouse homologies and the human-mouse one-to-one table, UCSC hg38 cytobands |
| `download_atlas_file` | `download_atlas.py` | one atlas file under `raw/<atlas>/`, checked against its SHA-256 in `config/atlases.yaml` |

### Stage 3: pilot (`features.smk`)

| Rule | Script | Reads | Writes |
|---|---|---|---|
| `stage3_pilot` | `pilot.py` | Yang 2022 human trigeminal ganglion counts | `results/pilot/`, `results/figures/pilot_specificity_auroc.png` |

Report: `docs/stage_reports/stage3_pilot.md`.

### Stage 4: features and model (`features.smk`)

| Rule | Script | Reads | Writes |
|---|---|---|---|
| `standardize_atlas` (one job per atlas, 11 atlases) | `standardize.py`, which calls `rds_to_mtx.R` / `seurat_to_mtx.R` for RDS files | raw atlas files, HGNC, mouse gene table, orthology table | `interim/<atlas>/cells.h5ad`, pseudobulk sums and group sizes |
| `stage4_model` | `stage4_model.py` | pseudobulk of the 11 atlases, `celltype_map.tsv`, gnomAD, GTEx, truth set | `processed/features_static.tsv.gz`, `results/tables/stage4/` (nested-CV performance, paired differences, coefficients, current ranking) |

Report: `docs/stage_reports/stage4_features_model.md`.

### Stage 5: validation (`features.smk`)

| Rule | Script | Reads | Writes |
|---|---|---|---|
| `stage5_backtest` | `stage5_backtest.py` | feature matrix, truth set | `results/backtest/backtest_performance.tsv`, `backtest_paired_differences.tsv`, `interval_simulation.tsv`; per-cutoff scores in `results/backtest/scores/` (not committed) |
| `stage5_repeated_cv` | `stage5_repeated_cv.py` | feature matrix, truth set | `results/backtest/repeated_cv_tier1.tsv` |
| `leakage_tests` | `pytest tests/test_leakage.py` | | `logs/leakage_tests_passed.txt` |
| `stage5_summarise` | `stage5_summarise.py` | back-test and interval tables, feature matrix | `interval_summary.tsv`, `interval_paired.tsv`, `failure_analysis.tsv`, `failure_feature_contrast.tsv`, back-test figures in `results/figures/` |

Report: `docs/stage_reports/stage5_validation.md`.

### Comparators in the back-test (`comparators.smk`)

| Rule | Script | Reads | Writes |
|---|---|---|---|
| `comparators_fetch` | `comparators_fetch.py` | `config/comparators.yaml` | GOA, GO, BioGRID and UniProt files under `external/comparators/`, with `manifest.json` |
| `comparators_backtest` | `comparators_backtest.py` | dated releases, cached back-test scores, feature matrix | `results/backtest/comparators/` performance, paired differences, coverage and random-seed null tables |
| `comparators_summarise` | `comparators_summarise.py` | the tables above | `comparators_auroc.png`, `comparator_headline_differences.tsv` |

Report: `docs/stage_reports/comparators.md`.

### Stage 6: released ranking (`tool.smk`)

| Rule | Command or script | Reads | Writes |
|---|---|---|---|
| `stage6_ranking` | `pnprank build-ranking` | feature matrix, truth set | `results/tables/pnprank_ranking.tsv.gz` and `_meta.json` |
| `stage6_checks` | `stage6_checks.py` | the ranking | `results/tables/stage6/release_checks.tsv` |

The `pnprank query` and `pnprank top` commands read the committed ranking. Report:
`docs/stage_reports/stage6_cli.md`.

### Case studies (`case_studies.smk`)

| Rule | Script | Reads | Writes |
|---|---|---|---|
| `locus_case_study` | `locus_case_study.py` | `config/loci.tsv`, cytobands, cached back-test scores; interval genes from the Ensembl REST API, cached under `external/loci/` | `results/case_studies/locus_ranking.tsv`, `locus_top10.tsv` |
| `top50_literature_queries` (live) | `top50_literature_queries.py` | the ranking; PubMed, HGNC, PanelApp and ClinGen queried live | `results/tables/top50_literature_queries.tsv`, `top50_literature_curated.tsv` |
| `celegans_orthologues` (live) | `celegans_orthologues.py` | the ranking; HGNC, Ensembl and Alliance APIs, cached under `external/celegans/` | `results/tables/celegans_shortlist.tsv`, `tier1_celegans_orthologues.tsv` |

The A to E literature class of each top-50 gene (`results/tables/top50_literature.tsv`)
was assigned by reading the search hits; the script only reproduces the searches.
Write-ups: `docs/case_studies/linkage_loci.md`, `top50_literature.md`,
`celegans_shortlist.md`.

### Manuscript

`manuscript/` holds the paper draft and the collaborator summary (Word files with an
EndNote library). They are written by hand from the reports above and are not built by
the workflow. `manuscript/figures/plain_figures.py` redraws the summary figures from the
committed result tables.

### Presentations

`presentations/<yyyy-mm>_<audience>/` holds one talk each. The animations are drawn by a
script in the same folder that reads `results/` and the data root, and the deck is
assembled by a second script; see the folder's README.

## Dependency outline

```
truth set (manual + CLI) ─┐
fetch_external ───────────┼─> standardize_atlas x11 ─> stage4_model ─┬─> stage5_backtest ─┬─> stage5_summarise
download_atlas_file ──────┘          │                               │                    ├─> comparators_backtest ─> comparators_summarise
                                     └─> stage3_pilot                │                    └─> locus_case_study
                                                                     ├─> stage5_repeated_cv
                                                                     └─> stage6_ranking ─┬─> stage6_checks
                                                                                         ├─> top50_literature_queries (live)
                                                                                         └─> celegans_orthologues (live)
```

`snakemake --rulegraph | dot -Tpng > rulegraph.png` draws the exact graph.
