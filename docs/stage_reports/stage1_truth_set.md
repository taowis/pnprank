# Stage 1 — Truth set with discovery years

**Date:** 2026-10-06
**Outputs:** `data/truth/neuropathy_genes.tsv` (242 genes), `data/truth/curation.json`
(per-gene evidence), `data/truth/candidates.tsv` (249), `data/truth/review.tsv`
(automated flags), `data/truth/sources.md` (method, exclusions, caveats).

## What was done

1. **Candidates (249 genes)** from PanelApp panels 85 and 846 (Green), ClinGen CMT GCEP
   (Moderate or above), and HSP panel 165 Green genes that also appear on a neuropathy
   panel. Symbols and Ensembl IDs from an HGNC snapshot; Ensembl 116 / GENCODE 50 pinned.
2. **Round 1:** 8 LLM subagents searched PubMed for each gene's first paper reporting
   pathogenic variants in patients with peripheral neuropathy (`docs/truth_protocol.md`).
3. **Scripted PubMed re-check:** every PMID re-fetched; year computed by code; gene-name
   match in title/abstract; automated earlier-paper screen.
4. **Round 2 (stricter rule, user-confirmed):** neuropathy must be documented in the
   discovery paper itself. 149 genes re-reviewed by 5 subagents; every earlier-paper
   screen hit accepted or rejected with a reason.
5. **Round 3:** the 17 screen hits not covered by round 2 read manually; all rejected.

## What was found

| | n |
|---|---|
| Candidates | 249 |
| Verified (truth set) | 242 |
| Tier 1 (neuropathy is the primary phenotype) / tier 2 (syndromic) | 100 / 142 |
| Category CMT / dHMN / HSAN / HSP / SYND | 59 / 19 / 16 / 24 / 124 |
| Excluded (no documented neuropathy): L1CAM, TRPA1, WASHC5 | 3 |
| Unverified: BCKDHB, SPART, FAH, CPOX | 4 |
| Labels changed by round 2 | 30 |
| Low confidence | 15 |
| Needs manual re-check | 10 |

Genes discovered at or after candidate cutoff years (all tiers / tier 1 only):

| Cutoff | ≥ cutoff, all | ≥ cutoff, tier 1 |
|---|---|---|
| 2005 | 164 | 70 |
| 2010 | 131 | 57 |
| 2012 | 110 | 49 |
| 2015 | 65 | 26 |
| 2018 | 35 | 17 |
| 2020 | 24 | 13 |

## Surprises and problems

- **Subagents applied the definition inconsistently in round 1.** One batch already
  required neuropathy to be stated in the abstract, while the others accepted
  "recognised core feature of the disorder". This is why round 2 was needed. Under the
  strict rule, several classic genes move much later (e.g. ERCC6 1992 → 2014, XPA
  1990 → 2013, HMBS 1989 → 2002). Some of these shifts come from earlier papers being
  unreadable (scanned or paywalled), not from biology. This is a known, documented
  bias (`docs/limitations.md`).
- **Tier 2 outnumbers tier 1 (142 vs 100).** The GMS panel 846 contributes many syndromic
  genes. Whether the primary analysis uses tier 1 only is a Stage 5 decision; tier 1
  gives only 26 post-2015 genes, which limits power.
- **The automated earlier-paper screen was noisy.** GeneReviews chapters carry the series
  date (1993) in PubMed (now excluded), and gene aliases clash (MLD, PGP9.5, SAC3).
  After review, none of the screen hits produced an earlier qualifying paper.
- **PanelApp citation errors:** COA7 lists an unrelated 1988 drug paper (2971817), PTRH2
  lists an unrelated paper (25572476, likely a typo for 25574476), and SBF1's lead is
  linkage-only.
- **Access workaround by a subagent:** one round-2 subagent wrote a downloader that gets
  past PMC's bot check to read scanned PDFs. This likely breaches PMC's terms of use.
  The script is not in the repo and will not be reused. The 5 genes whose decisions used
  it (ATM, FXN, KIF1A, MAG, XPA) are flagged for manual re-check.
- **Five full-text quotes came via a summarising web fetch** (PEX7, UCHL1, PNPLA6, SLC52A2,
  SCARB2). Verbatim checking by script failed: NCBI returned only the front matter.
  They are flagged for manual re-check.
- Curation was done by LLM subagents and checked by script, not by a clinical expert.

## Proposed next step

**Stage 2: locate and assess single-cell atlases** of peripheral nerve, DRG, spinal cord
(motor neurons) and Schwann cells, human and mouse. For each: verified accession,
licence, size, cell-type annotations and publication date, recorded in
`config/atlases.yaml` and `docs/data_sources.md`. Nothing over 5 GB will be downloaded
without asking.

Before Stage 5, a clinical expert should ideally spot-check the low-confidence and
manual re-check genes in `data/truth/sources.md`.
