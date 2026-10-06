# Truth set sources

Generated 2026-10-06. `neuropathy_genes.tsv` is built by `pnprank assemble-truth` from
`curation.json` (the full per-gene evidence record) and `candidates.tsv`.

## How it was built

1. Candidates (249): PanelApp Green genes on panels 85 and 846, ClinGen CMT GCEP
   curations at Moderate or above, and HSP-panel (165) Green genes also on a neuropathy panel.
   Snapshots and checksums: `data/external/truth_sources/manifest.json` (not in git).
2. Round 1: per-gene discovery-paper search under `docs/truth_protocol.md`.
3. Round 2: stricter rule (neuropathy documented in the paper itself) and review of every
   hit from an automated earlier-paper PubMed screen (`docs/truth_protocol_round2.md`).
4. Round 3: remaining unreviewed screen hits read manually (`round3_screen_review`).
5. `discovery_year` is computed by code from PubMed esummary: min(print, electronic) year.

## Not in the truth set (exclude from positives and negatives)

| Gene | Status | Reason |
|---|---|---|
| BCKDHB | unverified | No PubMed abstract reports BCKDHB variants in a patient with peripheral neuropathy. MSUD neuropathy leads: Kleopa 2001 Muscle Nerve (11180212, acute axonal neuropathy in MSUD; it says PN was 'not identified previously as |
| CPOX | unverified | Accepted only because acute motor neuropathy is a recognised feature of acute hepatic porphyria attacks (HCP); the abstract does not document neuropathy in this patient. Martasek et al. 1994 HMG. The first explicit HCP n |
| FAH | unverified | Paper states it is the first identification of a causal FAH mutation. PN qualifies as a recognised feature of HT1 (acute porphyria-like neurologic crises with axonal peripheral neuropathy; Mitchell 1990 NEJM PMID 2153931 |
| L1CAM | exclude | No PubMed record found reporting L1CAM variants in patients with documented peripheral neuropathy; peripheral neuropathy is not a recognised core feature of L1 syndrome (hydrocephalus/MASA/SPG1/CRASH). Gene appears via t |
| SPART | unverified | No abstract documents peripheral neuropathy in SPART-mutation patients. Troyer syndrome's cardinal features are spastic paraparesis, pseudobulbar palsy and distal amyotrophy (lower motor neuron). Neuropathy is not stated |
| TRPA1 | exclude | Only Mendelian TRPA1 disease is familial episodic pain syndrome 1 (Kremeyer et al., Neuron 2010 Jun 10, PMID 20547126: N855S gain-of-function; 'Quantitative sensory testing demonstrated normal baseline sensory thresholds |
| WASHC5 | exclude | Discovery paper Valdmanis 2007 AJHG (17160902) describes 'relatively severe, pure spastic paraplegia' - no PN. No SPG8 abstract (2013-2021 case reports/series checked: 23455931, 23881105, 25454649, 26572744, 28181327, 29 |

## Needs manual re-check

| Gene | PMID | Why |
|---|---|---|
| ATM | 8755918 | decision used PMC PDFs fetched by an automated downloader; re-check by opening the PMC article |
| FXN | 9150176 | decision used PMC PDFs fetched by an automated downloader; re-check by opening the PMC article |
| KIF1A | 21820098 | decision used PMC PDFs fetched by an automated downloader; re-check by opening the PMC article |
| MAG | 26179919 | decision used PMC PDFs fetched by an automated downloader; re-check by opening the PMC article |
| PEX7 | 12522768 | full-text quote obtained via a summarising web fetch; verify verbatim |
| PNPLA6 | 18313024 | full-text quote obtained via a summarising web fetch; verify verbatim |
| SCARB2 | 18308289 | full-text quote obtained via a summarising web fetch; verify verbatim |
| SLC52A2 | 22740598 | full-text quote obtained via a summarising web fetch; verify verbatim |
| UCHL1 | 28007905 | full-text quote obtained via a summarising web fetch; verify verbatim |
| XPA | 24252196 | decision used PMC PDFs fetched by an automated downloader; re-check by opening the PMC article |

## Low-confidence labels

| Gene | Year | PMID |
|---|---|---|
| ADA2 | 2014 | 25278816 |
| AGXT | 2000 | 10693261 |
| ERCC6 | 2014 | 25376329 |
| GLA | 1999 | 10069717 |
| HMBS | 2002 | 11940707 |
| MT-TL1 | 1992 | 1436526 |
| MTTP | 2000 | 10679949 |
| PPOX | 2001 | 11286631 |
| SAMD9L | 2018 | 30046003 |
| SCN9A | 2007 | 17294067 |
| SPG21 | 2021 | 34492745 |
| TTC19 | 2011 | 21278747 |
| VAPB | 2006 | 16967488 |
| VPS13A | 2012 | 22246890 |
| XPA | 2013 | 24252196 |

## Labels changed by the stricter rule or round-2 review

| Gene | Round 1 PMID (date) | Final PMID | Final year |
|---|---|---|---|
| ABCD1 | 8441467 (1993 Feb 25) | 8297373 | 1994 |
| ADA2 | 24552285 (2014 Mar 6) | 25278816 | 2014 |
| ALDH18A1 | 11092761 (2000 Nov 22) | 15517380 | 2004 |
| APTX | 11586299 (2001 Oct) | 11586300 | 2001 |
| ARSA | 1670590 (1991 Jan 3) | 8334794 | 1993 |
| ATM | 7792600 (1995 Jun 23) | 8755918 | 1996 |
| CYP27A1 | 2019602 (1991 Apr 25) | 8309576 | 1994 |
| ERCC6 | 1339317 (1992 Dec 11) | 25376329 | 2014 |
| ERCC8 | 7664335 (1995 Aug 25) | 17084038 | 2006 |
| FXN | 8596916 (1996 Mar 8) | 9150176 | 1997 |
| GALC | 8297359 (1994 Jan 28) | 9236794 | 1997 |
| GLA | 2539398 (1989 Apr) | 10069717 | 1999 |
| HADHA | 7811722 (1994 Dec 08) | 9003853 | 1997 |
| HEXA | 2522679 (1989 Mar 17) | 7898712 | 1995 |
| HMBS | 2563167 (1989 Jan) | 11940707 | 2002 |
| LYST | 8717042 (1996 Jul 18) | 8896560 | 1996 |
| MTTP | 8361539 (1993 Sep 2) | 10679949 | 2000 |
| PPOX | 8852667 (1996 Mar) | 11286631 | 2001 |
| SACS | 10655055 (2000 Feb) | 12873855 | 2003 |
| SAMD9L | 31053103 (2019 May 3) | 30046003 | 2018 |
| SCN11A | 24036948 (2013 Nov) | 24776970 | 2014 |
| SCN9A | 21698661 (2012 Jan) | 17294067 | 2007 |
| SPG11 | 17322883 (2007 Mar) | 18067136 | 2007 |
| SPG21 | 35111129 (2022 Jan 17 (eCollection 2021)) | 34492745 | 2021 |
| TTPA | 7719340 (1995 Feb) | 8602747 | 1996 |
| TYMP | 9924029 (1999 Jan 29) | 10852545 | 2000 |
| UCHL1 | 35986737 (2022 Oct) | 28007905 | 2017 |
| VPS13A | 11381253 (2001 Jun) | 22246890 | 2012 |
| XK | 8004674 (1994 Jun 17) | 11761473 | 2001 |
| XPA | 1702221 (1990 Dec) | 24252196 | 2013 |
