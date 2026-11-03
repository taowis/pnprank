# Case study — literature check of the top 50 candidates

**Date searched:** 2026-10-08. **Table:** `results/tables/top50_literature.tsv` (one row per
gene). **Searches:** `workflow/scripts/top50_literature_queries.py` re-runs the scripted
queries and writes `results/tables/top50_literature_queries.tsv` (PubMed hit counts) and
`results/tables/top50_literature_curated.tsv` (PanelApp and ClinGen entries, HGNC aliases).
The classification is manual. Design choices are recorded in `docs/decisions.md`
(2026-10-08).

**Question:** for the 50 best-ranked genes of the release ranking that are not known tier-1/2
genes (`pnprank top -n 50`, ranks 1-58), which have already been linked to neuropathy, and
which are new candidates?

## Method

- **Gene list:** `pnprank top -n 50 --tsv` on `results/tables/pnprank_ranking.tsv.gz`. It
  contains 49 unlabelled genes and L1CAM (`listed_not_in_truth`).
- **PubMed:** five fixed queries per gene, searched in title/abstract. Each uses the HGNC
  symbol plus HGNC alias and previous symbols. Aliases are dropped if they are shorter than
  4 characters, non-ASCII, generic words (Beta2, beta-5), or another gene's approved symbol
  (e.g. TUBB as an alias of TUBB2A). The exact strings are in the queries table.
  - Q1: symbol AND neuropathy terms (neuropathy, CMT, HMN/dHMN, HSAN/HSN, spastic paraplegia,
    motor neuron disease, SMA).
  - Q2: as Q1, adding the aliases.
  - Q3: Q2 restricted to human genetic reports (variant terms AND patient/family terms).
  - Q4: related neurological diseases (ALS, ataxia, leukodystrophy, myopathy, epilepsy,
    brain malformation and so on) AND variant terms AND patient terms.
  - Q5: functional or model-organism work (Schwann, myelin, axon, sensory or motor neuron,
    DRG, peripheral nerve AND knockout, knockdown, mutant, zebrafish and so on). Q5 also
    searches the HGNC approved name and a common protein name (e.g. osteopontin).
- **Ad hoc queries:** where a class needed checking, I ran targeted queries, e.g. TARDBP
  with sensory neuronopathy terms or KCNA1 with neuromyotonia. Each is listed with its hit
  count in the `additional_queries` column.
- **Curated sources:**
  - Genomics England PanelApp API: all panels, with confidence level.
  - ClinGen gene-disease validity: the full CSV, created 2026-10-07.
- **Reading:** I read the titles of every Q3 record and of the relevant Q4/Q5 records, and
  the abstract of every PMID cited in the table.
- **No claims from memory:** every PMID, disease link and classification in the table rests
  on a record retrieved on 2026-10-08.

**Classes** (each gene gets the highest class its evidence supports, in the order A > E):
- **A:** established gene for an inherited peripheral neuropathy or HSP. This needs curated
  gene-disease validity (ClinGen, or PanelApp Green) or several unrelated families.
- **B:** reported in inherited neuropathy or HSP, but the evidence is limited (one family,
  case reports, or cohort variant screening), or neuropathy appears as part of a broader
  syndrome.
- **C:** Mendelian link to a related neurological disease (ALS, ataxia, CNS myelin disorder,
  brain malformation, epilepsy, myopathy), but not to peripheral neuropathy.
- **D:** only functional or model-organism evidence about axons, myelin, Schwann cells, or
  sensory or motor neurons. Expression alone (a marker gene, an atlas) does not count.
- **E:** nothing relevant found.
- A and B must rest on human genetic evidence.

## Counts

| Class | n | Genes |
|---|---|---|
| A | 3 | L1CAM, TUBB4A, GLDN |
| B | 4 | TUBA1A, TARDBP, TUBB2A, SCN7A |
| C | 14 | MBP, ACTG1, TUBB, KCNA1, MLIP, CLTC, YWHAG, MAP1B, EPB41L3, BLTP1, ANK3, ATP5F1B, RELN, KIF5C |
| D | 14 | PEX5L, BCAS1, SECISBP2L, ANXA2, UGT8, S100B, SNCG, LGALS1, GFRA3, S100A10, NCMAP, CALCA, SPP1, PKM |
| E | 15 | IL16, SLC36A2, TUBA1B, MTCL3, CALM2, SEPTIN3, HSPA12A, GAS2L3, CDH19, YWHAH, VAT1L, TMEM255A, FXYD7, TUBB4B, EHBP1 |

- **Already linked to a human Mendelian neurological disease (A-C):** 21 of 50.
- **Already reported in neuropathy or HSP (A-B):** 7 of 50.
- **Some published PNS biology but no human link (D):** 14 of 50.
- **New on this search (E):** 15 of 50.

## Class A: a possible gap in the truth set (truth set not edited)

| Gene | Rank | Evidence | Gap under the project's own rules? |
|---|---|---|---|
| GLDN | 57 | Lethal congenital contracture syndrome 11, ClinGen **Definitive**. Biallelic variants in 4 unrelated families. Patient sciatic nerve shows lengthened nodes of Ranvier. The authors call it an inherited nodopathy "among peripheral neuropathies" (PMID 27616481); confirmed in an AMC cohort (33820833). | **Yes, probably.** It is the same kind of disease as CNTNAP1 (tier 1, SYND), ADCY6 (tier 1) and ADGRG6 (tier 2), congenital neuropathies presenting as arthrogryposis. It was never a candidate because PanelApp lists it only on the arthrogryposis/fetal panels, not on neuropathy panel 85 or 846, and it has no ClinGen CMT GCEP curation. |
| L1CAM | 2 | SPG1/MASA/X-linked hydrocephalus. ClinGen L1 syndrome Definitive; PanelApp HSP panel 165 Green (PMID 7920659). | **No.** Stage 1 already reviewed and excluded it because no peripheral neuropathy is reported. This search also found none (7 hits, none relevant). |
| TUBB4A | 31 | HLD6/H-ABC and DYT4. PanelApp Green on the HSP panels (165, 567, 568), and several families with complicated HSP and hypomyelination (25772097, 37867417). | **No.** It is on no neuropathy panel, and no TUBB4A neuropathy report was found (0 hits). The truth set takes HSP genes only with a verified neuropathy. |

- **Proposed action (for the user to decide):** curate GLDN under the truth protocol, then
  decide whether L1CAM and TUBB4A should be counted as "known" when the top list is
  presented.
- **Effect on the model:** GLDN was a training 0 in the release model. Adding it would
  change the labels, so it would need a rebuild of the ranking.

## Class B: reported, but the evidence is limited

- **TUBA1A (rank 1):** one patient with spastic paraplegia and ataxia (37435044). The gene's
  main disease is tubulinopathy.
- **TUBB2A (13):** one patient with spastic paraplegia, sensorimotor polyneuropathy and
  ataxia (29547997).
- **TARDBP (12):** an ALS gene, with several case reports of sensory and motor neuronopathy /
  FOSMN (21294910, 30586030, 38841627). Borderline B/C.
- **SCN7A (15):** listed among sodium-channel variants in painful diabetic and small-fibre
  neuropathy cohorts, without segregation (37175987). Also has DRG functional data
  (23026072). Borderline B/D.

## Class E: new on this search

| Gene | Rank | GTEx decile | Non-atlas pct | Note |
|---|---|---|---|---|
| IL16 | 4 | 6 | 0.78 | cytokine studies only |
| SLC36A2 | 5 | 3 | 0.69 | Mendelian disease is renal iminoglycinuria |
| TUBA1B | 6 | 10 | 1.00 | screened and excluded as a dHMN II positional candidate in 2002 (12090300) |
| MTCL3 | 22 | 3 | 0.57 | one PubMed record in total |
| CALM2 | 26 | 10 | 1.00 | Mendelian disease is cardiac (LQT/CPVT) |
| SEPTIN3 | 27 | 8 | 0.93 | knockout has no CNS neuronal phenotype (17564677) |
| HSPA12A | 30 | 5 | 0.98 | expression in spinal anterior horn only |
| GAS2L3 | 33 | 3 | 0.62 | zebrafish brain morphogenesis only |
| CDH19 | 36 | 5 | 0.90 | Schwann cell precursor marker (expression), enteric neural crest migration |
| YWHAH | 41 | 10 | 0.95 | nothing relevant |
| VAT1L | 44 | 7 | 0.71 | nothing relevant |
| TMEM255A | 47 | 4 | 0.88 | nothing relevant |
| FXYD7 | 50 | 8 | 0.77 | brain Na,K-ATPase regulator |
| TUBB4B | 51 | 10 | 1.00 | Mendelian disease is a ciliopathy (PCD, LCA with deafness) |
| EHBP1 | 56 | 8 | 0.95 | dendrite formation in hippocampal neurons |

- **Most open candidates:** MTCL3, VAT1L, TMEM255A, GAS2L3, SLC36A2, CDH19 and FXYD7 have
  almost no published biology.
- **Driven by the atlas features:** for most of them the non-atlas percentile is low
  (0.57-0.90). They are high mainly because of their glial or neuronal expression pattern.
- **Interpretation:** these are new candidates on this search, not candidates with support.

## How the classes relate to the percentile columns

| Class | n | Median rank | Top GTEx decile (n) | Median non-atlas pct | Non-atlas pct < 0.9 (n) |
|---|---|---|---|---|---|
| A | 3 | 31 | 1 | 0.92 | 1 |
| B | 4 | 12.5 | 2 | 0.98 | 1 |
| C | 14 | 20.5 | 6 | 0.98 | 2 |
| D | 14 | 33 | 7 | 0.96 | 4 |
| E | 15 | 33 | 4 | 0.90 | 8 |

- **Highly expressed genes are not mostly class E.** Of the 20 genes in the top GTEx
  expression decile, 4 are E (TUBA1B, CALM2, YWHAH, TUBB4B); 7 are D and 9 are A-C.
  - The tubulins and cytoskeletal genes at the top (TUBA1A, TUBB, TUBB2A, TUBB4A, ACTG1,
    CLTC, MAP1B, KIF5C) mostly have Mendelian CNS disease (C), sometimes with HSP or
    neuropathy features (A/B).
  - So the ranking puts broadly "neuronal disease" genes near the top. This is consistent
    with the expression confound in `docs/limitations.md`: the model partly learns
    "neural, highly expressed, disease-prone" rather than "PNS-specific".
- **Class E is enriched for atlas-driven genes.** 8 of 15 E genes have a non-atlas
  percentile below 0.9, against 8 of 35 in other classes (Fisher exact p = 0.049; post hoc,
  not corrected).
  - Genes that rank high only because of the atlas features are the least studied. That
    fits genuine novelty, and it equally fits noise from the atlas features.
  - This search cannot tell the two apart.
- **Classes A-C, against D-E:** genes with a human Mendelian neurological link (A-C) rank
  higher (median rank 17 vs 33).
  - Their non-atlas percentile is rarely low (4 of 21 below 0.9, against 12 of 29 for D-E).
  - Mendelian disease genes tend to be constrained and highly expressed, which is what the
    non-atlas features capture.

## Limitations

- **Search coverage:** PubMed title/abstract only. Full texts, OMIM, conference abstracts
  and preprints were not searched.
  - MAP1B: a PanelApp Red entry ("giant axonal neuropathy-like phenotype") cites only a
    conference abstract that I could not retrieve. It is recorded as unverified.
- **Alias matching:**
  - Aliases shorter than 4 characters were dropped, so reports that use only such an
    alias would be missed. Protein-name synonyms were added only in Q5.
  - Some hit counts are inflated by aliases that are common words or other proteins'
    names. CALCA via "CGRP" (Q2 = 235) is the main case, and these hits are mostly
    acquired neuropathic-pain studies.
- **Borderline calls:** they are noted in the table (TARDBP B/C, SCN7A B/D).
  - The rule that expression alone does not make a gene D affects CDH19 and HSPA12A.
- **One searcher:** the classification was made by one person, without a second
  independent reviewer.
- **Counts drift:** PubMed counts will change. The queries table records the counts on
  2026-10-08.
