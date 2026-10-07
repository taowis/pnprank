# Case study — candidates with C. elegans orthologues

**Date:** 2026-10-08. Script: `workflow/scripts/celegans_orthologues.py`. Outputs:
`results/tables/celegans_shortlist.tsv` (74 genes) and
`results/tables/tier1_celegans_orthologues.tsv` (100 genes). Rule and reasons:
`docs/decisions.md` (2026-10-08). Sources and releases: `docs/data_sources.md`.

This is a hand-off list for a C. elegans collaborator. It is not a validation of the model,
and a high rank is not evidence that a gene causes neuropathy.

## Selection

- **Candidates:** the top 200 genes of the release ranking, with known tier-1/2 genes
  excluded (`pnprank top -n 200`).
- **Orthology:** each gene was queried in two sources:
  - Ensembl REST pan-taxonomic Compara (release 116);
  - Alliance of Genome Resources 9.1.0, which reports how many of 9 prediction methods
    support each human-worm pair.
- **Clear orthologue:** a pair counts as clear if either
  - Ensembl calls it one-to-one, or
  - Alliance marks it as the best match in both directions and at least 5 of the 9
    methods support it.
- **Less testable in worms:** genes whose top cell type is a Schwann-cell type or
  oligodendrocyte are flagged (`worm_testable = False`). Worms have no myelin.

## Counts

**Tier-1 genes (established neuropathy genes, n = 100):**
- 55 have a clear worm orthologue;
- 12 more have a worm orthologue supported by at least 5 methods that is not the
  reciprocal best match. Examples: MFN2 / fzo-1, PDK3 / pdhk-2, PRPS1 / prps-1,
  SYT2 / snt-1. In each of these, Alliance does not mark the pair as the best match in
  both directions, and the worm gene is also linked to human paralogues (MFN1; PDK1/2/4;
  PRPS1L1/PRPS2; SYT1/SYT5);
- 23 have only weaker calls;
- 10 have no call in either source: ARHGAP19, CADM3, FBLN5, FBXO38, GJB1, LRSAM1, MPZ,
  NGF, RETREG1, SH3TC2.

So "55 of 100 tier-1 genes have a clear worm orthologue (67 with at least 5 supporting
methods)" is the sourced figure. For 5 of the 10 genes with no call (ARHGAP19, GJB1,
MPZ, RETREG1, SH3TC2), the top cell type in the ranking table is a Schwann-cell type.

**Top 200 candidates:**

| Orthology class | Genes |
|---|---|
| Clear | 74 (66 distinct worm genes) |
| Supported (≥ 5 methods, not reciprocal best) | 20 |
| Weak | 61 |
| None | 45 |

- 22 of the 74 clear genes have a Schwann-cell or oligodendrocyte top cell type.
- That leaves **52 testable genes** (45 distinct worm genes). By top cell type: 41
  neuronal, 10 non-neural, 1 other glia.

The 20 "supported" candidates are not in the shortlist file. They include KIF1B / unc-104,
KIF21A / klp-12, TUBB4A / tbb-4, SPTBN1 / unc-70 and DPYSL2 / dhp-1. All 200 genes, with
their orthology class, are in `<data_root>/external/celegans/top200_all_orthology.tsv`.

## Top 20 testable candidates

Columns:
- **Rank:** rank among unlabelled genes.
- **Methods:** Alliance methods supporting the pair (of 9).
- **Expr-matched pct / Non-atlas pct:** the expression-matched percentile and the
  non-atlas model's percentile.
- **Shares worm gene with:** other human genes that Alliance links to the same worm gene
  with at least 5 methods.

| Rank | Human | Worm gene (WormBase) | Ensembl type | Methods | Expr-matched pct | Non-atlas pct | Top cell type | Shares worm gene with |
|---|---|---|---|---|---|---|---|---|
| 1 | TUBA1A | mec-12 (WBGene00003175) | not called | 5/9 | 1.000 | 1.000 | neuron_sensory | TUBA1B;TUBA1C;TUBA3C;TUBA3D;TUBA4A |
| 6 | TUBA1B | mec-12 (WBGene00003175) | not called | 5/9 | 0.997 | 0.999 | neuron_motor | TUBA1A;TUBA1C;TUBA3C;TUBA3D;TUBA4A |
| 7 | ACTG1 | act-2 (WBGene00000064) | many2many | 8/9 | 1.000 | 0.999 | neuron_motor | ACTB |
| 8 | TUBB | tbb-4 (WBGene00006538) | not called | 6/9 | 1.000 | 1.000 | neuron_motor | TUBB2A;TUBB2B;TUBB4A;TUBB4B |
| 12 | TARDBP | tdp-1 (WBGene00006514) | one2one | 9/9 | 0.995 | 0.979 | myonucleus_mtj | - |
| 13 | TUBB2A | tbb-4 (WBGene00006538) | not called | 6/9 | 1.000 | 0.986 | neuron_motor | TUBB;TUBB2B;TUBB4A;TUBB4B |
| 17 | CLTC | chc-1 (WBGene00011867) | one2many | 9/9 | 1.000 | 0.999 | neuron_motor | CLTCL1 |
| 26 | CALM2 | cmd-1 (WBGene00000552) | not called | 6/9 | 0.997 | 0.998 | myonucleus_type2 | CALM1;CALM3 |
| 42 | BLTP1 | lpd-3 (WBGene00003060) | one2one | 8/9 | 1.000 | 0.996 | neuron_cholinergic_interneuron | - |
| 46 | ATP5F1B | atp-2 (WBGene00000229) | one2one | 9/9 | 0.992 | 0.988 | neuron_motor | - |
| 51 | TUBB4B | tbb-4 (WBGene00006538) | not called | 6/9 | 0.989 | 0.997 | neuron_sensory | TUBB;TUBB2A;TUBB2B;TUBB4A |
| 58 | PKM | pyk-1 (WBGene00009126) | many2many | 8/9 | 0.987 | 0.967 | neuron_motor | PKLR |
| 60 | ACTB | act-2 (WBGene00000064) | many2many | 8/9 | 0.992 | 0.999 | vascular_mural | ACTG1 |
| 78 | NDRG4 | ndrr-1 (WBGene00013020) | one2many | 9/9 | 0.981 | 0.830 | neuron_sensory | NDRG1;NDRG2;NDRG3 |
| 79 | KLC1 | klc-2 (WBGene00002215) | one2one | 8/9 | 1.000 | 0.948 | neuron_motor | KLC2;KLC4 |
| 83 | DYNLT3 | dylt-1 (WBGene00008764) | not called | 7/9 | 0.995 | 0.811 | neuron_sensory | DYNLT1 |
| 87 | ATP5F1A | atp-1 (WBGene00010419) | one2one | 9/9 | 0.989 | 0.981 | neuron_motor | - |
| 89 | NAA10 | daf-31 (WBGene00000923) | many2many | 9/9 | 0.983 | 0.872 | myonucleus_mtj | NAA11 |
| 98 | SAT1 | D2023.4 (WBGene00008408) | one2many | 9/9 | 0.973 | 0.996 | immune | SAT2;SATL1 |
| 100 | PRKCA | pkc-2 (WBGene00004033) | many2many | 9/9 | 0.992 | 0.942 | neuron_sensory_np | PRKCB;PRKCG |

**Flagged as less testable** (clear orthologue, but the top cell type is Schwann-cell or
oligodendrocyte):

L1CAM, SLC36A2, KCNA1, KIF5C, EHBP1, QKI, STARD13, AK5, AKR1B1, SHC4, FUT8, SEC14L5, LCTL,
CHL1, PTPDC1, NRXN1, CRYAB, CTSD, SLC9A3, CPEB1, CNP, HSPG2.

## Caveats

- **The top of the list has two problems for worm work.**
  - **Highly expressed housekeeping and cytoskeletal genes.** Tubulins, actins,
    calmodulins, CLTC and ATP synthase subunits rank high partly because they are highly
    expressed (`docs/limitations.md`). A strong worm phenotype in these genes would not
    show that they are neuropathy genes.
  - **Shared worm genes.** Several human paralogues map to one worm gene, so a worm
    mutant tests the gene family, not the specific human gene:
    - TUBA1A and TUBA1B both map to mec-12;
    - TUBB, TUBB2A and TUBB4B map to tbb-4;
    - CALM1, CALM2 and CALM3 map to cmd-1;
    - ACTB and ACTG1 map to act-2.

    The column `human_genes_sharing_worm_gene` shows this for each gene. Genes with no
    shared paralogue and an Ensembl one-to-one call are the cleanest tests. Nine of the
    52 testable genes qualify: TARDBP, BLTP1, ATP5F1B, ATP5F1A, VPS35, MARCHF6, EPRS1,
    PET100 and SUPT6H.
- **The two sources sometimes disagree.**
  - For tubulins and calmodulins, Ensembl pan-Compara does not list the Alliance-chosen
    worm gene, so the Ensembl type is empty.
  - CNP / nonu-1 is clear only through the Ensembl one-to-one call. Only 2 Alliance
    methods support it.
- **Top cell type is a crude proxy.** It is the single highest-specificity atlas cell
  type. A gene flagged as glial may also act in neurons, and the reverse.
  - KIF5C and NRXN1 are among the genes flagged by this rule. Their expression in other
    cell types was not checked for this note.
  - TARDBP's top cell type is a muscle nucleus type.
- **Enrichment columns are record counts from Alliance, not curated findings.**
  - `worm_n_alleles` counts Alliance entries with `alleleCategory=allele`.
  - `worm_expr_nervous_system` counts expression annotations under nervous system. A zero
    means no annotation, not evidence that the gene is absent from neurons.
  - All 52 testable genes have at least one listed allele; 35 have a nervous-system
    expression annotation.

  Whether a given allele is a null, viable, or available from the CGC has not been
  checked.
- **Novelty is not checked.** Whether any candidate already has a published neuropathy
  link or worm neurodegeneration phenotype has not been checked; this needs a literature
  search with sources.
- **Live APIs.** The Ensembl pan-Compara endpoint returned intermittent HTTP 500/503
  errors during the run. Requests were retried and the responses cached. A rerun after
  the next Ensembl or Alliance release can change the calls.
