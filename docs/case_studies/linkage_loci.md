# Case study — ranking genes inside published linkage intervals

**Date:** 2026-10-08. Script: `workflow/scripts/locus_case_study.py`; loci and sources:
`config/loci.tsv`; outputs: `results/case_studies/locus_ranking.tsv`, `locus_top10.tsv`.
Design and locus selection: `docs/decisions.md` (2026-10-08).

For each locus, the gene later reported as the cause is ranked among the protein-coding
genes in the linkage interval, using scores from a model trained only on tier-1 genes known
before the gene was reported (the Stage 5 back-test scores at that year). Rank 1 = top.

| Locus | Interval (GRCh38 span) | Gene (cutoff) | Genes ranked | PNPrank | Non-atlas model | GTEx nerve |
|---|---|---|---|---|---|---|
| HSN1 | 9q22.1-q22.3 (12.0 Mb) | SPTLC1 (2001) | 78 | **12** | 30 | 31 |
| DI-CMTB | 19p13.2-p12 (17.3 Mb) | DNM2 (2005) | 398 | 51 | **36** | 53 |
| X-linked dHMN | Xq13.1-q21 (30.6 Mb) | ATP7A (2010) | 103 | **36** | 77 | 48 |
| CMTX4 | Xq24-q26 (21.5 Mb) | AIFM1 (2012) | 130 | 46 | 60 | 52 |
| CMTX6 | Xp22.11 (3.0 Mb) | PDK3 (2013) | 19 | **5** | 9 | 8 |
| CMTX3 | Xq26.3-q27.1 (6.7 Mb) | SOX3 (2026) | 50 | 34 | **21** | 25 |
| DHMN1 | 7q34-q36.2 (16.7 Mb) + UBE3C | UBE3C (2022) | 146 | **11** | 24 | 22 |

- **Relative rank** (rank / genes ranked; lower is better):
  - PNPrank: median 0.26; the gene is in the top half of its interval at 6 of 7 loci.
  - Non-atlas model: median 0.42.
  - GTEx nerve: median 0.40.
- **PNPrank vs the non-atlas model:** better at 5 of 7 loci, worse at DI-CMTB and CMTX3.
- **The gene was never ranked first.** The best placements:
  - UBE3C: 11 of 146 (top 8%);
  - SPTLC1: 12 of 78 (top 15%), with a model trained on only 12 known genes;
  - DNM2: 51 of 398 (top 13%).
- **CMTX3 / SOX3 is the clearest failure** (34 of 50, genome-wide percentile 0.11). SOX3
  is a developmental transcription factor with little expression in adult nerve, DRG or
  muscle atlases. The CMTX3 mechanism, altered SOX3 expression caused by a nearby
  insertion, is outside what an adult expression-based model can see.
- **Caveats:**
  - There are only 7 loci, all from one group's publications.
  - The intervals are whole cytobands, wider than the marker-defined intervals.
  - The CMTX4 and CMTX3 intervals overlap.
  - LOEUF is missing for the X-linked target genes.
