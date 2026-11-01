# Manual check of the truth set

`truth_manual_check.xlsx` (and its TSV twin) lists the 242 genes in
`data/truth/neuropathy_genes.tsv` with the PubMed record of each discovery paper, for
the author to check by hand against the rule in Methods. Built by
`workflow/scripts/truth_manual_check_sheet.py` on 2026-11-02 from PubMed esummary.

Rows are ordered by priority:

| Priority | Genes | Reason |
|---|---|---|
| A | 10 | flagged "Needs manual re-check" in `data/truth/sources.md` |
| B | 14 | other low-confidence labels in `sources.md` |
| C | 16 | other labels changed by the round-2 review |
| D | 202 | the rest |

The columns from `check_result` onward are filled in by the checker; the "how to check"
sheet gives the rule and the steps. When the check is finished, the filled sheet is
committed here, changed entries are applied to `neuropathy_genes.tsv` with the reason in
`sources.md`, and the truth set, back-test and leakage tests are rerun.
