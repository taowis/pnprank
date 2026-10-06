"""Build and validate the truth set of known inherited peripheral neuropathy genes."""

from __future__ import annotations

import csv
from pathlib import Path

REQUIRED_COLUMNS = [
    "hgnc_symbol",
    "ensembl_gene_id",
    "disease_category",
    "inheritance",
    "discovery_year",
    "pmid",
    "source_url",
    "notes",
]
DISEASE_CATEGORIES = {"CMT", "dHMN", "HSAN", "HSP"}


def load_truth(path: str | Path) -> list[dict[str, str]]:
    """Read the truth TSV as a list of row dicts."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def validate_truth(path: str | Path) -> list[str]:
    """Return a list of problems in the truth TSV; empty means valid."""
    problems: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            return [f"missing columns: {missing}"]
        seen: set[str] = set()
        for lineno, row in enumerate(reader, start=2):
            symbol = row["hgnc_symbol"].strip()
            if not symbol:
                problems.append(f"line {lineno}: empty hgnc_symbol")
            elif symbol in seen:
                problems.append(f"line {lineno}: duplicate symbol {symbol}")
            seen.add(symbol)
            if not row["ensembl_gene_id"].startswith("ENSG"):
                problems.append(f"line {lineno}: bad ensembl_gene_id {row['ensembl_gene_id']!r}")
            if row["disease_category"] not in DISEASE_CATEGORIES:
                problems.append(f"line {lineno}: unknown category {row['disease_category']!r}")
            year = row["discovery_year"].strip()
            if not (year.isdigit() and len(year) == 4):
                problems.append(f"line {lineno}: bad discovery_year {year!r}")
            if not row["pmid"].strip().isdigit():
                problems.append(f"line {lineno}: missing or non-numeric pmid")
            if not row["source_url"].strip().startswith("http"):
                problems.append(f"line {lineno}: missing source_url")
    return problems


def known_before(rows: list[dict[str, str]], cutoff_year: int) -> set[str]:
    """Ensembl IDs of genes with discovery_year strictly before cutoff_year."""
    return {r["ensembl_gene_id"] for r in rows if int(r["discovery_year"]) < cutoff_year}
