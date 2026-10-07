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

## 2026-10-07 — Stage 3 (pilot)

- **Pilot atlas (user-confirmed):** Yang 2022 human trigeminal ganglion (GSE197289),
  snRNA-seq raw counts and author labels (`subtype`, 16 cell types, 3 donors).
- **R for RDS input:** `r-base` and `r-matrix` added to the conda environment; the
  RDS is converted to Matrix Market by `workflow/scripts/rds_to_mtx.R`. The GEO file
  is gzip-compressed twice; the script strips one layer if `readRDS` fails.
- **Gene IDs:** atlas symbols mapped to Ensembl 116 via HGNC (approved, then unique
  previous symbol, then unique alias); ambiguous names dropped; rows mapping to the
  same Ensembl ID summed. 10,225 of 32,214 symbols are unmapped, almost all
  clone-based lncRNA names absent from HGNC.
- **Pseudobulk:** raw counts summed per cell type x donor, groups with < 20 cells
  dropped, CPM per group, mean over donors; cell types must appear in >= 2 donors.
- **Specificity:** EWCE-style fraction of a gene's expression in each cell type, plus
  tau on log1p(CPM).
- **Gene universe and negatives:** HGNC protein-coding genes present in the atlas.
  Negatives exclude all 249 candidates (including unverified/excluded ones).
- **Expression-matched control:** specificity is ranked within deciles of mean
  expression before computing AUROC, because known disease genes are more highly
  expressed than average and raw specificity partly tracks expression.
- **Pilot is not a back-test:** it uses all discovery years and reports no
  temporal metric. No multiple-testing correction; results are exploratory.

## 2026-10-07 — Stage 4 (features and model)

- **Primary positive set (user-confirmed):** tier 1 (100 genes); tier 1+2 is a
  sensitivity analysis.
- **iPain not downloaded (user-confirmed):** mouse DRG is covered by Krauter 2025;
  the 21 GB iPain file is not needed.
- **Yim 2022 excluded from features:** the deposited Seurat object holds raw counts but
  only numeric cluster IDs (res.0.6 / active identity); author cell-type names are not
  deposited. Re-annotating clusters ourselves would add an undocumented judgement, and
  mouse sciatic nerve is already covered by Stefanovska 2025 with author labels.
- **R isolation:** the conda R must not load packages from the user R library
  (binary-incompatible, segfaults); `R_LIBS_USER` is pointed to a non-existent path in the
  environment (`environment.yml` `variables`).
- **Yadav motoneuron QC (pre-specified rule):** a donor's "Motoneurons" are kept only if
  pseudobulk CPM of both SLC5A7 and PRPH is >= 10x the interneuron median. Donor 5
  (SLC5A7 6, PRPH 0, high MBP/PLP1/GFAP: glial contamination or mislabel) and donor 1
  (SLC5A7 ~3x) fail; 1,082 motoneurons from 5 donors remain
  (`results/tables/stage4/yadav_motoneuron_qc.tsv`).
- **Yang 2022 mouse:** naive animals only; "Control" and "PBS" are procedural controls
  for the headache models and are excluded. "Injured" neurons are excluded.
- **Kedlian 2024:** nuclei only (cells and nuclei are mixed in the object).
- **SPARC 465:** library used as replicate unit (library-to-donor mapping not deposited).
- **Stefanovska 2025:** raw counts from GEO (the Zenodo h5ad holds only normalised
  values); young (3-month) samples only; GEO files matched to h5ad labels by barcode
  overlap (each file matched exactly one sample suffix, 100% of its labelled cells).
- **Jung 2023 human DRG (added after the first Stage 4 report):** the deposited human
  object (GSE201586_humanallcelltypes.rds.gz, gzip applied twice; 9.6 GB serialised)
  loaded with about 3 GB resident memory, less than estimated. Raw integer counts
  (32,174 genes x 111,925 nuclei) and the author `celltype` column were exported. The
  object has one "DRG Neuron" class (no neuron subtypes), mapped to `neuron_sensory`.
  Donors come from the GEO series matrix (sample IDs in the per-sample file names; 18
  libraries from 7 donors); libraries of one donor are summed. FACS and density-gradient
  nuclei isolations are pooled. GEO lists hg19 for alignment; genes are matched by
  symbol through HGNC, so the assembly does not affect the Ensembl 116 keys.
- **Gene universe for modelling:** HGNC protein-coding genes seen in at least one atlas
  or in GTEx (19,228 before Jung 2023 was added; 19,229 after); genes with no expression
  data at all are dropped.
- **Stage 4 evaluation is nested CV on current labels, not a back-test:** it asks only
  whether atlas features add signal over non-atlas predictors (paired bootstrap on
  identical resamples).
- **Lockfile regenerated** with `conda list --explicit --md5` after adding `r-seuratobject`.

## 2026-10-07 — Stage 5 (validation), pre-specified before any back-test was run

- **Cutoff years:** 2005, 2008, 2011, 2014, 2017, 2020 (every 3 years from the first year
  with >= 30 pre-cutoff tier-1 genes; user-confirmed rule). Post-cutoff tier-1 genes:
  70, 62, 54, 32, 19, 13.
- **Training labels at cutoff c (no future information):** 1 = tier-1 genes discovered
  before c. Tier-2 genes discovered before c are dropped (known before c). Every other
  gene is 0, including genes discovered at or after c and unverified/excluded candidates,
  because at time c they were not known neuropathy genes and their current status is
  future information.
- **Out-of-fold scoring:** genes labelled 0 are split into 5 folds; each fold is scored by
  a model trained on the pre-cutoff positives and the other folds, so no evaluated gene
  is scored by a model that saw it as a training negative. Hyperparameters are tuned by
  inner CV on the training genes only. Label-derived co-expression uses only pre-cutoff
  tier-1 genes.
- **Evaluation at cutoff c:** positives = tier-1 genes discovered in or after c;
  negatives = protein-coding genes that are not candidates (current knowledge is allowed
  for evaluation). Primary: all post-cutoff genes. Secondary: genes discovered in
  [c, c + 5).
- **Primary model:** logistic regression with all features (chosen in Stage 4). Gradient
  boosting is not back-tested (worse in Stage 4 and much slower).
- **Comparators at every cutoff:** non-atlas logistic model; GTEx tibial nerve
  expression; GTEx mean expression; -LOEUF. Paired bootstrap differences on identical
  resamples.
- **Expression-matched evaluation:** every score is converted to its percentile within
  GTEx mean-expression deciles (computed over the universe) before AUROC, as in the pilot.
- **Interval simulation:** for each year y from 2005 to 2026 with tier-1 discoveries, a
  model is trained with cutoff y and each gene discovered in y is ranked among all genes
  scored at that cutoff (percentile). Per-category summaries (CMT, dHMN, HSAN) use these
  percentiles. Failure analysis lists the genes ranked lowest at discovery.
- **Repeated CV:** the Stage 4 current-label comparison (logistic full vs non-atlas) is
  repeated over 10 seeds to measure split-to-split variance.
- **Tier 1+2** is a sensitivity analysis with the same procedure (positives tier 1+2).
- **Deviation from the pre-specified design (found while running, before results for all
  cutoffs were seen):** pooling raw out-of-fold scores across the 5 scoring folds is
  invalid when the fold models choose different regularisation, because their predicted
  probabilities are on different scales. At cutoff 2008 one fold chose C = 1.0 and the
  others 0.001-0.01; median scores per fold ranged 0.000-0.114, and the pooled AUROC of
  the full model (0.675) was below every one of its within-fold AUROCs except the lowest
  (0.675-0.854). Primary analysis is therefore changed to percentiles within the fold
  model that scored each gene (each fold is a random fifth of the scored genes), applied
  to every score including baselines and to the expression-matched analysis (fold x
  expression decile) and the interval simulation. The pre-specified pooled analysis is
  still computed and reported (`transform = pooled`). The Stage 4 CV comparison is barely
  affected (tier 1 gain +0.048 pooled vs +0.044 within-fold; tier 1+2 +0.059 both),
  because the full model chose C = 0.001 in every fold there.
- **Reproducibility fix:** the training rows at a cutoff were assembled from a Python set,
  whose iteration order changes with string-hash randomisation between processes; inner
  CV folds, and so the chosen C, changed between runs with the same seed (e.g. cutoff
  2005 tier 1: pooled AUROC 0.762 vs 0.726). Rows now follow the label order, and a test
  runs the scoring under three hash seeds. The episode also shows that C selection is
  unstable with 30-40 positives; within-fold percentiles remove the scale problem, not
  this variance. All Stage 5 results were regenerated after the fix.

## 2026-10-07 — Stage 6 (command-line tool)

- **The released ranking reuses the back-tested procedure, not the Stage 4 final model.**
  The Stage 4 ranking (`ranking_current_logistic_full.tsv.gz`) scores genes with a model
  that saw them as training negatives, which pushes their scores down in-sample. The
  release instead applies the Stage 5 out-of-fold scoring to today's labels: tier-1 genes
  are positives, every other gene is scored by a fold model that did not train on it,
  and scores are reported as within-fold percentiles. `fit_fold_models` was factored out
  of `oof_scores_at_cutoff` so both use one code path (cutoff 2005 scores reproduced to
  1e-16).
- **Tier-2 genes** are left out of training (as in the back-test) and scored by all five
  fold models; each gets the mean of its percentiles against the genes each model scored.
  **Tier-1 genes** were training positives in every fold model, so they get no
  percentile and no rank; the CLI reports them as known.
- **Regularisation (user-agreed proposal):** C stays tuned by inner CV in each fold model,
  as back-tested. The reported percentile comes from the first seed (the validated
  procedure, one seed). Nine further seeds are run only to report the range of each
  gene's percentile and the chosen C per fold, so users can see how much a rank depends
  on the random split. Averaging over seeds was not adopted because it was not
  back-tested.
- **Outputs (user-agreed proposal):** each gene's percentile is shown next to its
  expression-matched percentile (within fold x GTEx mean-expression decile) and the
  non-atlas model's percentile, so a high rank driven mainly by expression level, or not
  helped by atlas features, is visible.
- **Gene lookup** accepts Ensembl IDs (version ignored), HGNC approved symbols, then
  HGNC previous and alias symbols; an alias shared by several genes is reported as
  ambiguous rather than guessed. Aliases come from the HGNC complete set snapshot used
  for the truth set.
- **`top` excludes known tier-1/2 genes by default** (they are already known);
  `--include-known` shows them. Unverified candidate-list genes are shown with status
  `listed_not_in_truth`.
- **The ranking table is committed** (`results/tables/pnprank_ranking.tsv.gz`, one row
  per protein-coding gene, about 1 MB) so `pnprank query` works without the raw data.
- **C at the grid boundary (checked, grid kept):** with today's 100 tier-1 positives every
  fold model and seed chose C = 0.001, the smallest value in the grid. An extended inner-CV
  grid on the full release training set peaked at C = 3e-4, with inner-CV AUROC higher
  than at 0.001 by 0.001-0.006 (three seeds). The grid is kept as back-tested; the gain
  from extending it is within split-to-split noise.
- **Post hoc release checks** (`workflow/scripts/stage6_checks.py`), not pre-specified:
  tier-2 genes never enter training, so their percentiles are a held-out check against
  unlabelled genes; and the share of highly expressed genes at the top of the list.
