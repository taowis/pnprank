# data/

Only this file and `truth/` are tracked in git. Everything else is regenerated.

| Directory | Contents | How to regenerate |
|---|---|---|
| `truth/` | Curated neuropathy gene list with discovery years and PMIDs (committed) | Curated by hand, see `truth/sources.md` |
| `raw/` | Atlas downloads, untouched | `snakemake download` (Stage 2+), checksums in `config/atlases.yaml` |
| `external/` | gnomAD constraint, GTEx, interaction DB, gene annotation | `snakemake external` (Stage 4) |
| `interim/` | Pseudobulk matrices, per-atlas intermediates | `snakemake pseudobulk` |
| `processed/` | Final gene-by-feature matrix per cutoff year | `snakemake features` |

The data root can be moved off the repo via `config/paths.yaml`.
