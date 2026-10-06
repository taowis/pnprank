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
