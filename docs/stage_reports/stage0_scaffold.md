# Stage 0 — Scaffold

**Date:** 2026-10-06

## Done
- Directory layout per project spec; git repository initialised.
- `environment.yml` + osx-arm64 explicit lockfile; `pyproject.toml` with `pnprank` CLI.
- Config templates: `config.yaml`, `paths.example.yaml`, `atlases.yaml` (empty, to be verified in Stage 2).
- Truth TSV schema and validator (`pnprank validate-truth`), with tests.
- Snakemake chosen as workflow engine (see `docs/decisions.md`).

## Open items
- `disk_limit_gb` in `config/paths.yaml` is a placeholder (200 GB) pending user confirmation.
- Ensembl/GENCODE release, cutoff years and model type are TODOs for later stages.

## Next step
Stage 1: build the truth set with verified discovery years and PMIDs.
