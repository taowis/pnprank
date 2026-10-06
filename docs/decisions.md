# Decisions

Each non-obvious choice, with its reason. Newest last.

## 2026-10-06 — Stage 0

- **Workflow engine: Snakemake** (not Make). File-based DAG with per-rule
  conda/config support and wildcards per atlas and cutoff year, which fits
  one-rules-file-per-stage under `workflow/rules/`. Uses `snakemake-minimal`
  to keep the environment small.
- **Packaging:** `pyproject.toml` with setuptools, `src/` layout; dependencies
  are managed by `environment.yml`, so `pip install -e . --no-deps`.
- **CLI uses argparse** (stdlib) to avoid an extra dependency at this stage.
- **Python 3.11** pinned for broad compatibility with scanpy/anndata on conda-forge.
- **Lockfile:** `conda-lock.osx-arm64.txt`, an explicit spec exported from the
  created environment (`conda list --explicit --md5`). It is platform-specific;
  other platforms rebuild from `environment.yml`.
- **Truth TSV schema:** `hgnc_symbol, ensembl_gene_id, disease_category,
  inheritance, discovery_year, pmid, source_url, notes`. Ensembl ID is the key;
  every row must carry a PMID and a verification URL.
- **Back-test cutoff split is strict:** a gene is "known before" year Y only if
  `discovery_year < Y`.
- **Random seed:** `20261006` (the scaffold date; arbitrary but fixed).
- **Private notes:** `private_notes/` and the root `CLAUDE.md` are gitignored and
  never pushed. The root `CLAUDE.md` is a symlink to `private_notes/CLAUDE.md`
  so the assistant loads it from the repo root.
