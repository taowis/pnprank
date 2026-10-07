# Comparators in the temporal back-test

**Date:** 2026-10-08 (after Stage 6)

**Question:** inside the Stage 5 back-test, how does PNPrank compare with existing
gene-prioritisation approaches? The same setup is used throughout:
- cutoffs 2005–2020;
- tier 1 as the primary set, tier 1+2 as the sensitivity set;
- the same evaluation labels;
- within-fold-percentile and fold × expression-decile transforms;
- the paired bootstrap from `pnprank.backtest`.

**Outputs (`results/backtest/comparators/`):**
- `comparator_performance.tsv`
- `comparator_paired_differences.tsv`
- `comparator_headline_differences.tsv`
- `comparator_coverage.tsv`
- `comparator_data_coverage.tsv`
- `comparator_random_seed_null.tsv`
- figure `comparators_auroc.png`

**Code:**
- `src/pnprank/comparators.py`
- `workflow/scripts/comparators_fetch.py`, `comparators_backtest.py` and
  `comparators_summarise.py`
- release choices in `config/comparators.yaml`

**Tests:**
- `tests/test_comparators.py`
- three new leakage tests in `tests/test_leakage.py`

`tests/test_leakage.py` passes in full (13 tests), as does the whole suite (47 tests).

## Step 1: feasibility

The core problem is leakage. Most tools draw on resources that today record the very
gene–disease links being predicted. A comparator counts as leakage-free only if both of
these hold:
- (a) it uses no disease labels, or uses only pre-cutoff known genes as seeds;
- (b) its data come from a release dated before the cutoff, or carry no disease
  information at all.

A subagent verified the tool facts and the PMIDs against PubMed and the full texts on
2026-10-08. I looked up the release archives myself, as listed under "Data" below.

| Comparator | Disease labels? | Dated pre-cutoff data? | Verdict | Source |
|---|---|---|---|---|
| **GO guilt-by-association** (ToppGene/Endeavour-style functional similarity) | Only through seeds; we use pre-cutoff genes only | **Yes**: EBI GOA human archive (releases from 2001-11-02); GO ontology archive (monthly from 2004-03-01) | **Run, all 6 cutoffs** | ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ ; release.geneontology.org |
| **Network propagation (RWR)** on BioGRID | Only through seeds | **Yes from 2006**: earliest BioGRID archive 2.0.17 (2006-03-28); IntAct archive from 2006-01-01; HPRD download gave HTTP 503 | **Run at 2008–2020; infeasible at 2005** | downloads.thebiogrid.org/BioGRID/Release-Archive/ ; ftp.ebi.ac.uk/pub/databases/intact/ |
| **gnomAD constraint** (−LOEUF, pLI, missense Z) | No: observed vs expected variation (PMID 32461654); disease genes used only for validation | gnomAD v4.1 is post-cutoff but label-free (same assumption as for the expression atlases) | **Run** (−LOEUF was already a Stage 5 baseline) | PMC7334197 |
| ToppGene as distributed (PMID 19465376) | Its features include "disease-gene associations", "Phenotype (human disease…)" and literature co-citation | No dated archive found | **Infeasible as-is**; its GO component is reimplemented above | PMC2703978 |
| Endeavour (PMIDs 16680138, 27131783) | 75 current data sources, including disease-derived ones | No archived data releases found | **Infeasible as-is**; same reimplementation | PMC4987917 |
| Exomiser / hiPHIVE / PhenIX (PMIDs 26562621, 23660285, 25078397) | Yes: OMIM/Orphanet gene–disease links and HPO annotations | Data bundles only from `1711` (2017) | **Infeasible** (leaks; no pre-2017 data) | data.monarchinitiative.org/exomiser/data/ |
| Phen2Gene (PMID 32500119) | Yes: HPO, OMIM, ClinVar, Orphanet, GeneReviews | Releases from 2019-04-23 | **Infeasible** | PMC7252576 |
| AMELIE (PMID 32434849) | Yes: trained on OMIM/HGMD/ClinVar | Not verified | **Infeasible** | PMC9366928 |
| GADO (PMID 31253775) | Yes: HPO-annotated genes as training sets | HPO release not verified; HPO GitHub assets only from 2022-02-14 | **Infeasible** | PMC6599066 |
| Episcore (PMID 29849042) | Yes: ClinGen HI score-3 genes and earlier HI gene sets as positives | Single static 2018 release | **Infeasible** (label-trained) | PMC5976622 |
| pHaplo / pTriplo (PMID 35917817) | Partly: gold standard includes DDG2P/ClinGen HI genes; fitted to case–control CNV data | Single release 2022-03-11 (Zenodo 6347673) | **Infeasible** | PMC9742861 |
| HIPred (PMID 28137713), GHIS (PMID 26001969) | Yes: known HI genes (Dang et al.) as positives | Static 2015/2017 | **Infeasible** | PMC5581952, PMC4551909 |
| GeneBayes s_het (PMID 38977852) | No disease labels, but its prior uses current GO and PPI features, which indirectly encode disease knowledge | Single release (Zenodo 10.5281/zenodo.7939767) | **Not run**: −LOEUF and pLI already cover the constraint signal leakage-free | zenodo.org/records/7939767 |
| GeneMANIA (PMID 20576703) | Partly: OMIM/phenotype networks can be weighted in | Archives 2010-09-23 to 2021-04-19 | **Not run**: no pre-2008 archive; the network channels would need auditing; BioGRID RWR covers the PPI idea | genemania.org/data/archive/ |
| STRING | Partly: text mining and curated databases | Archives from v5.1 (2004-07-30) | **Not run**: needs the experimental channel only; a possible extension that would also cover PPI at 2005 | string-db.org/cgi/access |

## Step 2: implemented comparators (`src/pnprank/comparators.py`)

Methods were fixed before any result was seen (`docs/decisions.md`). They were not tuned.

**Data, one release per cutoff.** Each is the last release dated before 1 January of the
cutoff year; a test checks every date.

| Cutoff | GOA human release (date) | GO ontology | BioGRID release (releaseDate) |
|---|---|---|---|
| 2005 | 25 (2004-12-14) | 2004-12-01 | none (earliest 2006-03-28) |
| 2008 | 58 (2007-12-18) | 2007-12-01 | 2.0.36 (2007-12-15) |
| 2011 | 92 (2010-12-13) | 2010-12-01 | 3.1.72 (2010-12-25) |
| 2014 | 128 (2013-12-09) | 2013-12-01 | 3.2.108 (2013-12-26) |
| 2017 | 163 (2016-11-28) | 2016-12-01 | 3.4.144 (2016-12-30) |
| 2020 | 197 (2019-12-17) | 2019-12-09 | 3.5.180 (2019-12-31) |

**Where the dates come from:**
- GOA: the EBI archive listing.
- GO: the release directory name.
- BioGRID: the `<source releaseDate>` in each release's PSI-MI 2.5 file.

URLs and checksums are in `docs/data_sources.md`. Total download is about 400 MB, with no
single file over 76 MB. The data root is 19 GB of the 300 GB limit.

**Comparators:**
- **`go_similarity`.** GO terms are propagated over is_a and part_of using the ontology of
  the same era. Terms with 3–500 genes are kept, with IDF weights. The score is the cosine
  similarity to the centroid of the seed genes.
- **`ppi_rwr`.** Random walk with restart (restart probability 0.5) on human
  physical-interaction networks:

  | Cutoff | Edges | Genes |
  |---|---|---|
  | 2008 | 23k | 7.6k |
  | 2011 | 33k | 8.7k |
  | 2014 | 125k | 14.5k |
  | 2017 | 214k | 15.3k |
  | 2020 | 339k | 16.5k |

- **Label-free:** `pli` and `mis_z` (gnomAD v4.1). They sit alongside the Stage 5
  baselines `neg_loeuf` and `gtex_nerve`.
- **Seed-free study-bias controls:** `go_term_count` (number of GO terms per gene) and
  `ppi_degree`.

**Seeds** are exactly the PNPrank training positives at that cutoff (`training_labels`).
Genes that cannot be scored (unannotated, or outside the network) score 0, which puts them
in a tie at the bottom.

Coverage of post-cutoff tier-1 genes:

| Cutoff | GO | PPI |
|---|---|---|
| 2005 | 43/70 | — |
| 2008 | 51/62 | 40/62 |
| 2011 | 45/54 | 34/54 |
| 2014 | 29/32 | 26/32 |
| 2017 | 18/19 | 17/19 |
| 2020 | 13/13 | 13/13 |

**Leakage tests** (`tests/test_leakage.py`):
- Comparator scores are unchanged when post-cutoff truth is removed, moved later,
  extended or re-tiered.
- The seed set passed to the scoring function is exactly the pre-cutoff tier-1 genes.
- Seeds are never scored.
- A sanity test shows that the seeds do change the scores.
- The GAF reader drops annotations dated on or after the cutoff. No release contained
  any.

## Step 3: results

### Tier 1, all genes after the cutoff, within-fold percentile (primary)

AUROC [95% CI]. The PNPrank column reproduces the Stage 5 point estimates exactly.

| Cutoff | PNPrank | GO similarity | PPI RWR | Missense Z | pLI | GTEx nerve | Non-atlas |
|---|---|---|---|---|---|---|---|
| 2005 | **0.759** [0.710, 0.808] | 0.653 [0.588, 0.720] | — | 0.695 [0.621, 0.763] | 0.540 [0.453, 0.627] | 0.714 | 0.732 |
| 2008 | **0.773** [0.720, 0.828] | 0.706 [0.646, 0.762] | 0.643 [0.567, 0.713] | 0.685 [0.608, 0.756] | 0.546 [0.456, 0.630] | 0.714 | 0.746 |
| 2011 | **0.766** [0.702, 0.826] | 0.671 [0.594, 0.743] | 0.646 [0.573, 0.715] | 0.670 [0.583, 0.749] | 0.559 [0.469, 0.647] | 0.720 | 0.751 |
| 2014 | **0.781** [0.696, 0.856] | 0.716 [0.630, 0.791] | 0.622 [0.508, 0.727] | 0.604 [0.494, 0.710] | 0.499 [0.382, 0.610] | 0.726 | 0.710 |
| 2017 | **0.762** [0.650, 0.858] | 0.727 [0.622, 0.819] | 0.700 [0.569, 0.828] | 0.660 [0.516, 0.789] | 0.554 [0.401, 0.711] | 0.755 | 0.737 |
| 2020 | 0.753 [0.635, 0.866] | 0.720 [0.626, 0.800] | 0.710 [0.583, 0.830] | 0.630 [0.443, 0.804] | 0.507 [0.346, 0.665] | 0.748 | **0.755** |

Paired differences, PNPrank minus comparator (ΔAUROC [95% CI]). **Bold** marks a CI
that excludes 0.

| Cutoff | vs GO similarity | vs PPI RWR | vs missense Z | vs pLI |
|---|---|---|---|---|
| 2005 | **+0.106 [+0.034, +0.178]** | — | +0.064 [−0.014, +0.142] | **+0.219 [+0.123, +0.313]** |
| 2008 | +0.067 [−0.002, +0.137] | **+0.130 [+0.048, +0.209]** | +0.088 [+0.000, +0.171] | **+0.227 [+0.124, +0.335]** |
| 2011 | **+0.095 [+0.008, +0.186]** | **+0.120 [+0.045, +0.197]** | **+0.096 [+0.003, +0.192]** | **+0.207 [+0.100, +0.314]** |
| 2014 | +0.064 [−0.018, +0.148] | **+0.159 [+0.058, +0.265]** | **+0.177 [+0.058, +0.283]** | **+0.282 [+0.133, +0.422]** |
| 2017 | +0.035 [−0.079, +0.146] | +0.062 [−0.059, +0.189] | +0.102 [−0.021, +0.220] | **+0.207 [+0.046, +0.348]** |
| 2020 | +0.032 [−0.116, +0.167] | +0.043 [−0.125, +0.206] | +0.123 [−0.039, +0.286] | **+0.246 [+0.044, +0.420]** |

The 2008 missense-Z interval rounds to +0.000 at its lower end, so it is not marked as
excluding 0.

**Expression-matched (tier 1)**
- PNPrank: 0.67–0.73.
- GO similarity: 0.59–0.68, and it keeps its AUROC better than GTEx nerve (0.58–0.63).
- PPI RWR: 0.54–0.64.
- PNPrank minus GO, matched: +0.102 [+0.010, +0.188] at 2005; +0.002 to +0.076 at the
  other cutoffs, with CIs including 0.
- PNPrank minus PPI, matched: +0.15 to +0.19 with CIs excluding 0 at 2008–2014;
  +0.03 to +0.07 at 2017–2020, with CIs including 0.

### Tier 1+2 sensitivity (24–164 test genes)

| Comparator | AUROC range | PNPrank minus comparator | Cutoffs where the CI excludes 0 |
|---|---|---|---|
| PNPrank | 0.735–0.811 | — | — |
| GO similarity | 0.638–0.724 | +0.07 to +0.10 | 4 of 6 (2005–2014) |
| PPI RWR | 0.564–0.644 | +0.15 to +0.21 | all 5 |
| Missense Z | 0.56–0.60 | — | all 6 |
| pLI | 0.47–0.52 | — | all 6 |

Expression-matched: the GO CI excludes 0 at 3 of 6 cutoffs; PPI, missense Z and pLI
at every cutoff.

### Is the comparator signal specific to the seeds? (study-bias controls)

These are tier-1, raw AUROCs.

| Cutoff | GO similarity | GO, random seeds: mean [2.5–97.5%] | GO term count | PPI RWR | PPI, random seeds | PPI degree |
|---|---|---|---|---|---|---|
| 2005 | 0.640 | 0.606 [0.554, 0.647] | 0.681 | — | — | — |
| 2008 | 0.697 | 0.598 [0.552, 0.635] | 0.679 | 0.640 | 0.619 [0.611, 0.630] | 0.619 |
| 2011 | 0.677 | 0.588 [0.551, 0.619] | 0.627 | 0.619 | 0.603 [0.588, 0.621] | 0.638 |
| 2014 | 0.713 | 0.596 [0.525, 0.662] | 0.648 | 0.619 | 0.606 [0.590, 0.630] | 0.608 |
| 2017 | 0.729 | 0.606 [0.515, 0.678] | 0.638 | 0.696 | 0.682 [0.656, 0.715] | 0.694 |
| 2020 | 0.719 | 0.620 [0.521, 0.678] | 0.674 | 0.710 | 0.701 [0.657, 0.755] | 0.694 |

The GO-similarity and PPI-RWR columns are raw AUROCs, so they differ slightly from the
within-fold values in the main table. The random-seed null uses 20 draws. GO term count and
PPI degree are within-fold AUROCs.

## What was found

1. **PNPrank is better than or level with every leakage-free comparator at every cutoff.**
   No comparator beats it anywhere.
   - In tier 1 the difference excludes 0:
     - against GO similarity at 2 of 6 cutoffs;
     - against PPI propagation at 3 of 5;
     - against missense Z at 2 of 6;
     - against pLI at 6 of 6.
   - In tier 1+2 it excludes 0 at 4 of 6 against GO and at every cutoff against the others.
2. **GO guilt-by-association is the strongest existing approach.**
   - It catches up with PNPrank at the late cutoffs: 0.72–0.73 at 2017–2020, against
     0.75–0.76 for PNPrank.
   - After expression matching it equals PNPrank at 2017–2020 (0.68 vs 0.69; 0.67 vs
     0.67).
   - Its real seeds beat random seeds at 2008–2020, so part of its signal is genuine
     functional similarity to known neuropathy genes.
   - At 2005 it is no better than annotation count (0.653 vs 0.681) or random seeds.
     The 2004 GOA release is sparse: about 10.9k genes are annotated, mostly by IEA.
3. **Network propagation on dated BioGRID adds almost nothing beyond node degree.**
   - Seeded RWR is within, or barely above, the random-seed null at every cutoff.
   - It is close to plain degree: 0.64 vs 0.62 at 2008, 0.71 vs 0.69 at 2020.
   - Its apparent rise at late cutoffs follows network growth (23k to 339k edges), not
     neuropathy-specific information.
   - After expression matching it falls to 0.54–0.64.
4. **Label-free constraint scores are weak for this disease group.**
   - pLI is near chance (0.50–0.56), which fits the many dominant missense and recessive
     neuropathy genes.
   - Missense Z (0.60–0.70) is the better constraint score in tier 1, but worse in tier 1+2.

## Surprises and problems

- **GO similarity at late cutoffs.** It becomes competitive with PNPrank at 2017–2020,
  both raw and expression-matched.
  - The test sets are small (19 and 13 genes).
  - Annotation of the remaining "undiscovered" genes is more complete by then.
  - The honest reading: PNPrank's lead over GO-based ranking is clear at early cutoffs
    and not demonstrable at late ones.
  - Whether GO and atlas features are complementary (for example, adding GO similarity as
    a feature) was not tested. It would be a new model and would need its own pre-specified
    back-test.
- **The PPI result looked like a plausible baseline at first (0.70–0.71 at 2017–2020).**
  The random-seed and degree controls show this is study bias. The controls were added
  after the first run, because RWR is known to track degree. They are reported for every
  cutoff, as recorded in `docs/decisions.md`.
- **Identifier loss in old GO releases.**
  - 2004 and 2007 GOA files key proteins by UniProt accessions, many now deleted. 70%
    (2005) and 61% (2008) of GOA protein IDs map to a current gene.
  - Mapping goes through UniProt idmapping and secondary accessions. These are identifier
    maps only.
  - This handicaps GO at early cutoffs. Unmapped proteins are mostly deleted TrEMBL
    entries, but the bias cannot be measured exactly.
- **No pre-2006 interaction release in BioGRID or IntAct, and HPRD is down.** PPI is
  therefore missing at 2005. STRING v5.1 (2004-07-30) exists, but only its experimental
  channel would be leakage-free; this was not attempted.
- **Nothing looked too good.**
  - No comparator exceeded PNPrank.
  - PNPrank point estimates match Stage 5 exactly.
  - The scored gene sets were asserted to be identical to the cached back-test scores.

## Proposed next step

For the paper:
- Report the tier-1 table and the paired-difference table above, with the feasibility
  table as a supplement.
- State that phenotype-driven tools cannot be back-tested without leakage.

Optional:
- A pre-specified back-test of PNPrank + dated GO similarity.
- A STRING v5.1 experimental-channel network to cover PPI at 2005.
