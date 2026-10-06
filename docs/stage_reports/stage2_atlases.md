# Stage 2 — Locate and assess atlases

**Date:** 2026-10-07
**Outputs:**
- `config/atlases.yaml`: 13 shortlisted atlases, status `shortlisted`, nothing downloaded.
- `results/tables/atlas_candidates.tsv`: all 133 candidates with suitability ratings.
- `docs/data_sources.md`: verified accessions and licences.
- `docs/atlas_search_protocol.md`: the search protocol.
- `data/external/orthology/human_mouse_one2one_e116.tsv`: not in git; regenerable.

## What was done

1. Queried the CELLxGENE curation API for nerve, ganglion, spinal cord and muscle tissues.
2. Ran three scoped searches by LLM subagents under a written protocol. They used GEO,
   PubMed, Europe PMC, Zenodo, SPARC, ArrayExpress and CELLxGENE. Every accession, PMID,
   file size and licence had to be confirmed on the official page. No data was downloaded,
   and no access control was circumvented.
3. Independently re-queried every shortlisted accession (GEO E-utilities, CELLxGENE,
   Pennsieve, Zenodo) and checked the linked PMIDs on PubMed.
4. Set up mouse-to-human orthology from Ensembl Compara 116, using one-to-one
   orthologues only.

## What was found

133 candidate datasets: 21 high, 43 medium, 69 low suitability. Shortlist of 13:

| Cell types | Human | Mouse |
|---|---|---|
| Schwann (myelinating / Remak), nerve fibroblasts, nerve vascular and immune | Heming 2025 sural nerve, 4 control donors, ~34k nuclei | Yim 2022 (4 nerves, healthy); Stefanovska 2025 (young sciatic) |
| Sensory neuron subtypes, satellite glia | Jung 2023 DRG; Yang 2022 TG; SPARC 465 sympathetic ganglia | Krauter 2025 DRG (22 neuron subtypes); Yang 2022 TG; iPain (naive) |
| Motor neurons | Yadav 2023 lumbar cord (1,331 labelled "Motoneurons", needs QC) | Alkaslasi 2021 (5,054 skeletal MNs, alpha/gamma) |
| Muscle, NMJ, terminal/perineurial cells | Kedlian 2024 | Ham 2025 (NMJ myonuclei) |

- The 13 shortlisted entries total about 89 GB as listed. About 26 GB of that is ATAC
  fragment files that are not needed, so the planned download is about 63 GB. That is
  within the 200 GB budget.
- One file exceeds the 5 GB single-download limit: the iPain DRG all-cells h5ad (21 GB).
  It needs approval before download.
- Orthology: 16,335 human–mouse one-to-one pairs (15,696 high confidence).
  - 232 of the 242 truth genes have a one-to-one mouse orthologue.
  - The 10 without one are ADA2, CD59, DEGS1, FDXR, MORC2, MT-RNR1, MT-TL1, PRPS1,
    SLC52A2 and SMN1. These are mostly paralog expansions in mouse, or mitochondrial
    non-coding genes. In mouse atlases they will have missing values, not zeros.

## Surprises and problems

- **Human peripheral nerve is barely covered.**
  - The only usable healthy human nerve single-cell data is 4 sural-nerve control donors,
    all from one centre (Heming 2025).
  - There is no healthy human single-cell data for mixed motor nerves (sciatic, tibial).
  - Tabula Sapiens, Tabula Muris and the GTEx single-nucleus data contain no peripheral nerve.
  - GTEx tibial nerve exists only as bulk RNA-seq.
- **Human motor neurons are scarce.**
  - Yadav 2023 is the only healthy human dataset with deposited motor neuron labels.
  - Its count (1,331) is far higher than other human datasets (around 100–650), so the
    cluster must be quality-checked before use.
  - No human dataset separates alpha from gamma motor neurons in deposited labels.
- **Key resources are not versioned.**
  - The largest human DRG reference (Bhuiyan 2025, 126 donors) is distributed only as a
    Dropbox Seurat object, with no licence or checksum. It is excluded from the shortlist
    for reproducibility.
  - The Zeisel 2018 mousebrain.org downloads return 404.
- **Metadata errors:**
  - GEO GSE167597 links the wrong PMID.
  - GEO GSE241386 appears to link a different paper.
  - Some author metadata files have column headers shifted by one (GSE124312, Kupari 2021).
- **Some data is controlled or embargoed:**
  - Tavares-Ferreira 2022 human DRG (dbGaP).
  - SPARC 464 human DRG multiome (embargoed to 2027-04-29).
  - A Biogen ALS spinal cord atlas (not yet public).
- **Ensembl changes:**
  - BioMart no longer responds.
  - The release-116 homo_sapiens homology dump contains no Mus musculus rows; mouse
    pairs are only in the mus_musculus dump. Using the human dump silently returned
    zero orthologues. A test now covers both gene orders.
- **Not yet confirmed** for some shortlisted objects (Yim 2022 RDS, Jung 2023 RDS):
  - whether they hold raw counts or normalised values;
  - exact per-cell-type counts.

  These are checked on download.
- **Licences:** GEO pages state no licence; these are recorded as "not stated". Only CELLxGENE
  (CC-BY 4.0), Zenodo and SPARC records state one.

## Proposed next step

**Stage 3: feasibility pilot on one atlas.** Recommended: **Yang 2022 human trigeminal
ganglion** (GSE197289, about 0.7 GB, 3 donors). It is small, has deposited raw counts and
verified labels, and covers both mechanisms in one dataset: neurons (axonal CMT, HSAN)
and myelinating and non-myelinating Schwann cells (demyelinating CMT). The alternative is
Heming 2025 sural nerve controls (about 1.5 GB), which has better Schwann and fibroblast
resolution but no neurons.

The pilot would:
1. Download the atlas and record checksums.
2. Build per-cell-type pseudobulk profiles.
3. Compute per-gene specificity.
4. Run an exploratory check that known neuropathy genes are enriched in the expected
   cell types.

No back-test numbers will be reported before `tests/test_leakage.py` covers features.
