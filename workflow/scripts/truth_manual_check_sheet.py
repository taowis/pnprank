"""Build the sheet for the author's manual check of the truth set.

One row per gene in data/truth/neuropathy_genes.tsv, with the PubMed record of the
discovery paper and empty columns for the check. Genes flagged in data/truth/sources.md
come first. Writes data/truth/manual_check/truth_manual_check.xlsx and a TSV twin.

Usage: python workflow/scripts/truth_manual_check_sheet.py
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.worksheet.datavalidation import DataValidation

ROOT = Path(__file__).resolve().parents[2]
TRUTH = ROOT / "data/truth/neuropathy_genes.tsv"
SOURCES = ROOT / "data/truth/sources.md"
OUT = ROOT / "data/truth/manual_check"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def section_genes(md: str, heading: str) -> list[str]:
    """Gene symbols in the first column of the table under a '## heading'."""
    m = re.search(rf"^## {re.escape(heading)}\n(.*?)(?=^## |\Z)", md, re.S | re.M)
    if not m:
        raise SystemExit(f"section not found in sources.md: {heading}")
    rows = [r for r in m.group(1).splitlines() if r.startswith("| ")]
    return [r.split("|")[1].strip() for r in rows[1:]]  # rows[0] is the header


def pubmed_summaries(pmids: list[str]) -> dict:
    out = {}
    for i in range(0, len(pmids), 150):
        q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids[i:i + 150]),
                                    "retmode": "json"})
        with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
            res = json.load(r)["result"]
        for p in res.get("uids", []):
            out[p] = res[p]
        time.sleep(0.4)
    return out


def main():
    genes = list(csv.DictReader(TRUTH.open(), delimiter="\t"))
    md = SOURCES.read_text()
    recheck = section_genes(md, "Needs manual re-check")
    lowconf = section_genes(md, "Low-confidence labels")
    changed = section_genes(md, "Labels changed by the stricter rule or round-2 review")
    summ = pubmed_summaries(sorted({g["pmid"] for g in genes}))

    def priority(sym):
        if sym in recheck:
            return "A", "flagged for manual re-check in sources.md"
        if sym in lowconf:
            return "B", "low-confidence label in sources.md"
        if sym in changed:
            return "C", "label changed in round 2"
        return "D", ""

    rows = []
    for g in genes:
        s = summ.get(g["pmid"], {})
        ids = {a["idtype"]: a["value"] for a in s.get("articleids", [])}
        pr, why = priority(g["hgnc_symbol"])
        rows.append({
            "priority": pr, "why_flagged": why,
            "hgnc_symbol": g["hgnc_symbol"], "ensembl_gene_id": g["ensembl_gene_id"],
            "tier": g["tier"], "disease_category": g["disease_category"],
            "inheritance": g["inheritance"], "discovery_year": g["discovery_year"],
            "pmid": g["pmid"], "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{g['pmid']}/",
            "pmc_url": (f"https://www.ncbi.nlm.nih.gov/pmc/articles/{ids['pmc']}/"
                        if "pmc" in ids else ""),
            "first_author": (s.get("authors") or [{}])[0].get("name", ""),
            "journal": s.get("source", ""), "pubdate": s.get("pubdate", ""),
            "epubdate": s.get("epubdate", ""), "title": s.get("title", ""),
            "curation_note": g["notes"],
            "check_result": "", "correct_pmid": "", "correct_year": "",
            "earlier_report_found": "", "evidence": "", "checked_on": "",
        })
    order = {"A": 0, "B": 1, "C": 2, "D": 3}
    rows.sort(key=lambda r: (order[r["priority"]], r["tier"], r["hgnc_symbol"]))
    cols = list(rows[0])

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "truth_manual_check.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    wb = Workbook()
    ws = wb.active
    ws.title = "check"
    ws.append(cols)
    for r in rows:
        ws.append([r[c] for c in cols])
    head = Font(bold=True)
    fill_check = PatternFill("solid", fgColor="FFF2CC")
    for c in ws[1]:
        c.font = head
    first_check = cols.index("check_result") + 1
    for col in range(first_check, len(cols) + 1):
        ws.cell(1, col).fill = fill_check
    widths = {"why_flagged": 28, "title": 60, "curation_note": 60, "evidence": 50,
              "pubmed_url": 38, "pmc_url": 45}
    for i, c in enumerate(cols, 1):
        ws.column_dimensions[ws.cell(1, i).column_letter].width = widths.get(c, 14)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.alignment = Alignment(vertical="top", wrap_text=cols[c.column - 1] in widths)
        for c in (row[cols.index("pubmed_url")], row[cols.index("pmc_url")]):
            if c.value:
                c.hyperlink = c.value
                c.style = "Hyperlink"
    letter = ws.cell(1, first_check).column_letter
    dv = DataValidation(type="list", formula1='"agree,change,exclude,unsure"', allow_blank=True)
    ws.add_data_validation(dv)
    dv.add(f"{letter}2:{letter}{len(rows) + 1}")
    ews = cols.index("earlier_report_found") + 1
    el = ws.cell(1, ews).column_letter
    dv2 = DataValidation(type="list", formula1='"no,yes"', allow_blank=True)
    ws.add_data_validation(dv2)
    dv2.add(f"{el}2:{el}{len(rows) + 1}")
    ws.freeze_panes = "D2"
    ws.auto_filter.ref = ws.dimensions

    hw = wb.create_sheet("how to check")
    lines = [
        "Manual check of the PNPrank truth set (242 genes)",
        "",
        "Rule. The discovery year is the year of the earliest PubMed record reporting "
        "pathogenic variants in the gene in patients with a documented peripheral neuropathy. "
        "Neuropathy must be documented in the discovery paper itself. The year is the earlier "
        "of the print and electronic publication years (pubdate, epubdate).",
        "Tier 1: neuropathy is the primary phenotype. Tier 2: syndromic, or neuropathy "
        "secondary or atypical.",
        "",
        "For each row:",
        "1. Open the paper (pubmed_url, or pmc_url for full text).",
        "2. Confirm it reports pathogenic variants in this gene in patients, and that "
        "peripheral neuropathy is documented in those patients.",
        "3. Search PubMed (and Google Scholar if useful) for an earlier paper that would meet "
        "the rule. Record yes/no in earlier_report_found.",
        "4. check_result: agree / change (fill correct_pmid and correct_year) / exclude "
        "(gene does not meet the rule) / unsure.",
        "5. evidence: the sentence or section that supports the decision (short quote and "
        "where it is, e.g. 'Results, para 2' or 'Table 1').",
        "6. checked_on: date of the check (YYYY-MM-DD).",
        "",
        "Order. Priority A (genes flagged for manual re-check in sources.md) first, then B "
        "(other low-confidence labels), then C (other labels changed in round 2), then D. "
        "A gene on more than one list is placed in the highest group.",
        "Do not edit the columns left of check_result; they record what the curation found.",
    ]
    for ln in lines:
        hw.append([ln])
    hw.column_dimensions["A"].width = 120
    for c in hw["A"]:
        c.alignment = Alignment(wrap_text=True, vertical="top")
    hw["A1"].font = head
    wb.save(OUT / "truth_manual_check.xlsx")
    counts = {p: sum(r["priority"] == p for r in rows) for p in "ABCD"}
    missing = [r["pmid"] for r in rows if not r["title"]]
    print("rows", len(rows), counts, "pmc", sum(bool(r["pmc_url"]) for r in rows),
          "missing_pubmed", missing)


if __name__ == "__main__":
    main()
