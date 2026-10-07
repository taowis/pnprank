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
- Yim 2022 (mouse nerves) is excluded because only numeric clusters were deposited; Jung
  2023 (human DRG) is not yet included (memory limit).
- SPARC 465 uses library, not donor, as the replicate unit.
