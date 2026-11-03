# Limitations

- Expression atlases published after a back-test cutoff are used as inputs. This
  is assumed acceptable because expression atlases do not encode disease labels.
- Discovery years require neuropathy to be documented in the discovery paper
  (abstract or accessible full text). Where an earlier paper's full text could not
  be read (paywalled, scanned only), the year falls to the first verifiable paper.
  This can shift some labels later than the true first report; such papers are
  listed as `inaccessible_earlier` in `data/truth/curation.json`.
- Discovery-paper curation was performed by LLM subagents under a written
  protocol and checked by script (PMID existence, PubMed dates, gene mention,
  earlier-paper screen). It was not independently curated by a clinical expert.
- Expression level is a confound: known neuropathy genes are more highly expressed than
  average, and the current model also ranks broadly expressed housekeeping and
  cytoskeletal genes highly. AUROC partly reflects expression level.
- Heming 2025 sural nerve counts are CellBender-corrected (as deposited), not raw
  CellRanger counts.
- Atlases differ in protocol (cells vs nuclei), species and depth; mouse genes enter only
  through one-to-one orthologues, so genes without one are missing from mouse features.
- Yim 2022 (mouse nerves) is excluded because only numeric clusters were deposited. Jung
  2023 (human DRG) contributes one pooled DRG-neuron profile; its human object has no
  neuron subtypes.
- SPARC 465 uses library, not donor, as the replicate unit.
- Stage 4 model comparisons use a single 5-fold CV split; split-to-split variance is not
  included in the bootstrap CIs. Repeating the split over 10 seeds (Stage 5,
  `results/backtest/repeated_cv_tier1.tsv`) gave an SD of 0.009 for the full-model AUROC
  and 0.008 for the atlas gain.
- Back-test cutoffs share test genes, so per-cutoff results are not independent and are
  not combined into one significance test. Late cutoffs (2017, 2020) have 19 and 13 test
  genes and wide CIs.
- Hyperparameter selection is unstable with few pre-cutoff positives (30-40); back-test
  scores are compared within the fold model that produced them, a deviation from the
  pre-specified pooled analysis (both are reported).
- Genes without a one-to-one mouse orthologue have no features from the five mouse
  atlases (e.g. PRPS1, MORC2), which lowers their rank.
- Only protein-coding genes are ranked. All 100 tier-1 training genes and 240 of the 242
  truth-set genes are protein-coding (the other two, MT-RNR1 and MT-TL1, are mitochondrial
  RNA genes), so there are no non-coding examples to learn from. A non-coding RNA gene
  that causes neuropathy would be missed; a non-coding variant acting on a nearby
  protein-coding gene can still be reached through that gene's rank.
- Genes acting in rare cell types with low bulk expression (e.g. SLC5A7 in motor neurons)
  are ranked low because overall expression carries much of the model's weight.
- The released ranking (Stage 6) is not itself back-tested: it applies the back-tested
  procedure to today's labels. Its percentiles mean "rank among genes scored by the same
  fold model", not a probability of being a neuropathy gene.
- Expression matching uses GTEx mean-expression deciles; within the top decile, very highly
  expressed genes (tubulins, actin, clathrin) still rank at the top, so the
  expression-matched percentile does not remove this confound at the extreme. 40 of the
  top 100 ranked genes are in the top expression decile.
- Tier-1 genes are training positives and receive no percentile in the release.
- Comparators (`docs/stage_reports/comparators.md`) are re-implementations of the GO
  similarity and network-propagation ideas behind existing tools on dated releases, not
  the tools themselves; phenotype-driven tools that rely on disease annotations cannot be
  back-tested without leakage and are not compared. The PPI comparator has no pre-2005
  release and is not run at 2005. Old GO releases lose 30-40% of proteins to identifiers
  that cannot be mapped today.
- Co-expression features (the only label-derived features) were at first computed once
  from all positives of a training set, before the inner cross-validation that chooses C,
  so inner-validation positives were seeds while C was chosen (outer test and evaluation
  genes never were). Corrected on 2026-11-03: `fit_tuned` recomputes them inside each
  inner split, and every affected analysis was rerun (`docs/decisions.md`).
- The back-test restricts the timing of labels only. Atlases (published 2021-2025),
  GTEx v10 and gnomAD v4.1 are current releases at every cutoff, while the GO and
  BioGRID comparators use dated releases; the comparison favours PNPrank to an unknown
  degree.
- Evaluation background genes are unlabelled, not known negatives, and are chosen with
  the current candidate lists. Widening the background to all scored genes changed
  tier-1 AUROC by at most 0.002 (`results/tables/revision_v3/sensitivity_eval_only.tsv`).
- Bootstrap CIs hold the fitted models fixed; they do not include variation from
  refitting, fold assignment, truth-set curation or discovery-year assignment.
