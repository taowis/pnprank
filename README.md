# pnprank

**Status:** Stage 2 (Atlases) — complete, awaiting review before Stage 3.

| Stage | Description | Status |
|---|---|---|
| 0 | Scaffold | done |
| 1 | Truth set with discovery years | done (242 genes) |
| 2 | Locate and assess atlases | done (13 shortlisted of 133) |
| 3 | Feasibility pilot on one atlas | not started |
| 4 | Full feature set and model | not started |
| 5 | Validation: temporal back-test, baselines, interval simulation, failure analysis | not started |
| 6 | Command-line tool | not started |

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

## Layout

- `config/` — run settings, atlas registry, machine paths template
- `data/` — not in git except `data/truth/` (curated labels) and `data/README.md`
- `src/pnprank/` — the importable package and `pnprank` CLI
- `workflow/` — Snakemake rules (one file per stage) and thin wrapper scripts
- `results/` — small tables and figures (committed); models and matrices are not
- `docs/` — decisions, data sources, limitations, per-stage reports
- `tests/` — including `test_leakage.py`, which must pass before any back-test result is reported

## Data policy

Public data only. No patient or lab-internal data enters this repository.

## Licence

Code: MIT (see `LICENSE`). Data: each source under its own licence, see `docs/data_sources.md`.
