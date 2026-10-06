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
- **Disk budget:** `disk_limit_gb: 200` under the data root, confirmed by the user;
  raised to 300 GB by the user on 2026-10-07. Single-download limit stays 5 GB.
- **Remote:** private GitHub repository; `private_notes/` and root `CLAUDE.md` are never pushed.

## 2026-10-06 — Stage 1 (truth set)

- **Scope (user-confirmed):** core set = PanelApp Green genes on the hereditary
  neuropathy panels ∪ ClinGen CMT GCEP curations at Moderate/Strong/Definitive.
  HSP genes are included only if Green on the HSP panel (165) *and* present on a
  neuropathy panel at any confidence, and only if a neuropathy phenotype is
  verified in the literature; flagged `also_hsp`.
- **Two PanelApp neuropathy panels:** 85 (Hereditary neuropathy, 100K legacy,
  v1.513) and 846 (Hereditary neuropathy or pain disorder, NHS GMS R78, v8.37).
  846 is the actively maintained clinical panel and contains recently discovered
  genes (e.g. CADM3, ITPR3, COQ7) that 85 lacks; omitting it would bias the
  back-test against recent discoveries.
- **Tiers:** 846 also contains syndromes where neuropathy is a minor feature. Each
  gene gets `tier` 1 (neuropathy is the primary phenotype) or 2 (syndromic). The
  primary analysis can use tier 1 with tier 1+2 as a sensitivity analysis
  (decided in Stage 5).
- **Category `SYND`** added for syndromic disorders with neuropathy.
- **Discovery year (user-confirmed):** year of the earliest PubMed record
  reporting pathogenic variants in the gene in human patients with a peripheral
  neuropathy phenotype. Linkage-only, animal-only and association studies do not
  count. A gene first reported for a non-neuropathy disease takes the year of its
  first neuropathy report.
- **Year from PubMed date:** min(print year, electronic year) taken from PubMed
  esummary by code, never from the curator; a discovery is public from whichever
  date came first.
- **Symbols and Ensembl IDs from HGNC**, not from PanelApp (PanelApp symbols lag,
  e.g. AARS vs AARS1).
- **Reference:** Ensembl 116 / GENCODE 50 (GRCh38.p14), the current release on
  2026-10-06.
- **Literature curation method:** per-gene search by LLM subagents under a fixed
  protocol (`docs/truth_protocol.md`), followed by a scripted PubMed re-check of
  every PMID (existence, date, gene mentioned) and manual review of flagged rows.
- **Stricter discovery rule (user-confirmed, after round 1):** neuropathy must be
  documented in the discovery paper's own genotyped patients (abstract or full
  text). "Recognised core feature of the disorder" no longer suffices. Round 1
  subagents applied the looser rule inconsistently (one batch already used the
  strict one), so ~30 genes are re-curated in round 2 (`docs/truth_protocol_round2.md`).
- **Secondary or atypical neuropathies kept as tier 2 (user-confirmed):** acute
  porphyrias (HMBS, CPOX, PPOX), ADA2, AGXT, SMN1, PDYN and PRNP are not
  excluded; they are forced to tier 2.
- **Independent earlier-paper screen:** every verified pick is re-searched on
  PubMed (gene names AND neuropathy terms AND mutation/variant, published before
  the pick's year). 130/244 picks had hits; each hit is read and accepted or
  rejected with a reason in round 2.
- **Unverified and excluded candidates** are removed from the positive set and
  must also be removed from the negative set in modelling (they are neither
  confirmed positives nor clean negatives).

## 2026-10-07 — Stage 2 (atlases)

- **Search method:** CELLxGENE curation API plus three scoped literature/GEO searches
  by LLM subagents (peripheral nerve; sensory and autonomic ganglia; spinal motor
  neurons and muscle) under `docs/atlas_search_protocol.md`. Subagents were told not
  to download data or circumvent any access control. 133 candidate datasets are in
  `results/tables/atlas_candidates.tsv`; raw subagent records are kept (not in git)
  under `data/external/atlas_search/subagent_results/`.
- **Independent check:** every shortlisted GEO accession was re-queried via E-utilities
  (title, taxon, linked PMID), and CELLxGENE, SPARC and Zenodo records via their APIs.
- **Shortlist (13 entries in `config/atlases.yaml`, status `shortlisted`):** chosen for
  healthy/control cells, raw counts, deposited author annotations and coverage of
  Schwann cells, nerve fibroblasts, sensory neuron subtypes, satellite glia, motor
  neurons and NMJ. Not yet downloaded.
- **Excluded from the shortlist:** Dropbox-only releases (e.g. the Bhuiyan 2025 human
  DRG reference) because they have no version, checksum or licence; controlled-access
  data (Tavares-Ferreira 2022, dbGaP); embargoed data (SPARC 464 until 2027-04-29);
  spot-level spatial data; datasets without healthy controls.
- **Disease datasets:** only control/naive cells enter features (e.g. Heming 2025 uses
  its 4 control donors; iPain uses naive animals). Patient cells are not used.
- **Licences:** GEO pages state no licence; recorded as "not stated". CELLxGENE data is
  CC-BY 4.0 per its publication terms (verified on the CELLxGENE docs page).
- **Orthology:** Ensembl Compara 116 protein_default homologies, one-to-one only. BioMart
  is no longer reachable, and the homo_sapiens dump omits mouse in release 116, so the
  mus_musculus dump is used. The parser accepts either gene order (Compara README).
- **Large files:** iPain DRG all-cells h5ad (21 GB) exceeds the 5 GB single-download
  limit and needs approval. ATAC fragment files (SPARC 465 17.8 GB, Krauter 8 GB) are
  not needed and will not be downloaded.
