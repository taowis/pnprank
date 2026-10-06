"""Build and validate the truth set of known inherited peripheral neuropathy genes."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import date
from pathlib import Path

REQUIRED_COLUMNS = [
    "hgnc_symbol",
    "ensembl_gene_id",
    "disease_category",
    "tier",
    "also_hsp",
    "inheritance",
    "discovery_year",
    "pmid",
    "source_url",
    "notes",
]
# SYND: syndromic disorder in which peripheral neuropathy is a reported feature.
DISEASE_CATEGORIES = {"CMT", "dHMN", "HSAN", "HSP", "SYND"}
TIERS = {"1", "2"}  # 1: neuropathy is the primary phenotype; 2: syndromic

PANELAPP_API = "https://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/"
CLINGEN_CSV = "https://search.clinicalgenome.org/kb/gene-validity/download"
HGNC_TSV = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)
CLINGEN_MIN_CLASSIFICATIONS = {"Definitive", "Strong", "Moderate"}


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
            if row["tier"] not in TIERS:
                problems.append(f"line {lineno}: bad tier {row['tier']!r}")
            if row["also_hsp"] not in {"yes", "no"}:
                problems.append(f"line {lineno}: bad also_hsp {row['also_hsp']!r}")
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


# --- Candidate list from curated public sources -------------------------------


def download(url: str, dest: Path) -> dict[str, str]:
    """Download url to dest with curl (PanelApp rejects urllib's user agent).

    Returns a manifest entry with the URL, retrieval date and sha256.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-sSfL", "-m", "600", "-o", str(dest), url], check=True)
    sha = hashlib.sha256(dest.read_bytes()).hexdigest()
    return {"file": dest.name, "url": url, "retrieved": date.today().isoformat(), "sha256": sha}


def fetch_sources(outdir: Path, panel_ids: dict[str, int]) -> list[dict[str, str]]:
    """Snapshot PanelApp panels, the ClinGen gene-validity CSV and HGNC into outdir."""
    manifest = []
    for name, pid in panel_ids.items():
        url = PANELAPP_API.format(panel_id=pid)
        manifest.append({"source": f"panelapp:{name}", **download(url, outdir / f"panelapp_{pid}.json")})
    manifest.append({"source": "clingen", **download(CLINGEN_CSV, outdir / "clingen_gene_validity.csv")})
    manifest.append({"source": "hgnc", **download(HGNC_TSV, outdir / "hgnc_complete_set.txt")})
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def panelapp_genes(path: Path, min_confidence: int = 3) -> dict[str, dict]:
    """Genes on a PanelApp panel snapshot at or above min_confidence (3 = Green), by HGNC ID."""
    panel = json.loads(path.read_text())
    return {
        g["gene_data"]["hgnc_id"]: g
        for g in panel["genes"]
        if int(g["confidence_level"]) >= min_confidence
    }


def clingen_genes(path: Path, gcep_prefix: str) -> dict[str, list[dict[str, str]]]:
    """ClinGen curations from one expert panel at Moderate or stronger, by HGNC ID."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    header_idx = next(i for i, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
    header = rows[header_idx]
    out: dict[str, list[dict[str, str]]] = {}
    for r in rows[header_idx + 1 :]:
        if not r or r[0].startswith("+"):
            continue
        rec = dict(zip(header, r))
        if rec["GCEP"].startswith(gcep_prefix) and rec["CLASSIFICATION"] in CLINGEN_MIN_CLASSIFICATIONS:
            out.setdefault(rec["GENE ID (HGNC)"], []).append(rec)
    return out


def load_hgnc(path: Path) -> dict[str, dict[str, str]]:
    """HGNC complete set keyed by HGNC ID."""
    with open(path, newline="") as fh:
        return {r["hgnc_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def build_candidates(srcdir: Path, panel_ids: dict[str, int]) -> list[dict[str, str]]:
    """Union of neuropathy-panel Green genes and ClinGen CMT GCEP (>= Moderate) genes,
    plus HSP-panel Green genes that appear on a neuropathy panel at any confidence.

    Symbols and Ensembl IDs come from the HGNC snapshot, not from the panels.
    """
    hgnc = load_hgnc(srcdir / "hgnc_complete_set.txt")
    neuro = {k: v for k, v in panel_ids.items() if k != "hsp"}
    sources: dict[str, set[str]] = {}
    for pid in neuro.values():
        for h in panelapp_genes(srcdir / f"panelapp_{pid}.json"):
            sources.setdefault(h, set()).add(f"panelapp_{pid}_green")
    for h in clingen_genes(srcdir / "clingen_gene_validity.csv", "Charcot-Marie-Tooth"):
        sources.setdefault(h, set()).add("clingen_cmt_gcep")
    any_neuro = set()
    for pid in neuro.values():
        any_neuro |= set(panelapp_genes(srcdir / f"panelapp_{pid}.json", min_confidence=0))
    hsp = set(panelapp_genes(srcdir / f"panelapp_{panel_ids['hsp']}.json"))
    for h in hsp & any_neuro:
        sources.setdefault(h, set()).add(f"panelapp_{panel_ids['hsp']}_green+neuropathy_panel")

    rows = []
    for h, srcs in sources.items():
        rec = hgnc.get(h)
        if rec is None or rec["status"] != "Approved":
            raise ValueError(f"{h} not an approved HGNC ID in snapshot")
        rows.append(
            {
                "hgnc_id": h,
                "hgnc_symbol": rec["symbol"],
                "ensembl_gene_id": rec["ensembl_gene_id"],
                "aliases": "|".join(
                    s for s in (rec["prev_symbol"] + "|" + rec["alias_symbol"]).split("|") if s
                ),
                "gene_names": "|".join(n for n in [rec["name"], *rec["prev_name"].split("|")] if n),
                "also_hsp": "yes" if h in hsp else "no",
                "candidate_sources": ";".join(sorted(srcs)),
            }
        )
    return sorted(rows, key=lambda r: r["hgnc_symbol"])


# --- PubMed verification ------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _get(url: str) -> bytes:
    import time

    time.sleep(0.4)  # NCBI limit without an API key: 3 requests/second
    return subprocess.run(["curl", "-sSfL", "-m", "120", url], check=True, capture_output=True).stdout


def pubmed_records(pmids: list[str]) -> dict[str, dict]:
    """Fetch summary and abstract text for PMIDs from PubMed."""
    out: dict[str, dict] = {}
    for i in range(0, len(pmids), 100):
        chunk = ",".join(pmids[i : i + 100])
        summ = json.loads(_get(f"{EUTILS}esummary.fcgi?db=pubmed&retmode=json&id={chunk}"))
        text = _get(f"{EUTILS}efetch.fcgi?db=pubmed&rettype=abstract&retmode=text&id={chunk}")
        abstracts = _split_abstracts(text.decode(errors="replace"), pmids[i : i + 100])
        for pmid in pmids[i : i + 100]:
            rec = summ["result"].get(pmid, {})
            out[pmid] = {
                "exists": "error" not in rec and bool(rec),
                "title": rec.get("title", ""),
                "journal": rec.get("source", ""),
                "pubdate": rec.get("pubdate", ""),
                "epubdate": rec.get("epubdate", ""),
                "abstract": abstracts.get(pmid, ""),
            }
    return out


def _split_abstracts(text: str, pmids: list[str]) -> dict[str, str]:
    """efetch text output is one record per block ending in 'PMID: <id>'."""
    import re

    out: dict[str, str] = {}
    start = 0
    for m in re.finditer(r"PMID:\s*(\d+)", text):
        if m.group(1) in pmids and m.group(1) not in out:
            out[m.group(1)] = text[start : m.end()]
            start = m.end()
    return out


def discovery_year(rec: dict) -> int:
    """Earliest of the print and electronic publication years.

    The discovery is public from whichever came first.
    """
    years = [int(d[:4]) for d in (rec["pubdate"], rec["epubdate"]) if d[:4].isdigit()]
    if not years:
        raise ValueError(f"no parsable date in {rec}")
    return min(years)


def mentions_gene(rec: dict, names: list[str], full_names: tuple[str, ...] = ()) -> bool:
    """True if any symbol/alias appears as a whole word, or any full gene name appears
    (case-insensitive), in the title or abstract."""
    import re

    text = rec["title"] + " " + rec["abstract"]
    if any(n.lower() in text.lower() for n in full_names if n):
        return True
    return any(re.search(rf"(?<![A-Za-z0-9]){re.escape(n)}(?![A-Za-z0-9])", text) for n in names if n)


PN_TERMS = (
    'neuropathy[tiab] OR neuropathies[tiab] OR "Charcot-Marie-Tooth"[tiab] OR CMT[tiab] '
    'OR HMSN[tiab] OR HSAN[tiab] OR HSN[tiab] OR dHMN[tiab] OR "distal spinal muscular atrophy"[tiab] '
    'OR "hereditary motor"[tiab] OR "hereditary sensory"[tiab] OR neuronopathy[tiab] '
    'OR polyneuropathy[tiab]'
)


def earlier_hits(names: list[str], before_year: int, limit: int = 5) -> list[str]:
    """PMIDs mentioning the gene and a neuropathy term, published before before_year.

    A screen only: hits are candidates for manual review, not evidence.
    """
    import urllib.parse

    gene = " OR ".join(f'"{n}"[tiab]' for n in names if len(n) > 2)
    # Exclude NCBI Bookshelf records: GeneReviews chapters carry the series date (1993).
    term = (
        f"({gene}) AND ({PN_TERMS}) AND (mutation*[tiab] OR variant*[tiab]) "
        "NOT pubmed books[filter]"
    )
    q = urllib.parse.urlencode(
        {"db": "pubmed", "retmode": "json", "retmax": limit, "sort": "pub_date",
         "datetype": "pdat", "mindate": "1900", "maxdate": str(before_year - 1), "term": term}
    )
    res = json.loads(_get(f"{EUTILS}esearch.fcgi?{q}"))
    return res["esearchresult"].get("idlist", [])


def assemble_truth(
    curation: list[dict], candidates: list[dict[str, str]], force_tier2: frozenset[str] = frozenset()
) -> tuple[list[dict], list[dict]]:
    """Check curated discovery papers against PubMed and build truth rows.

    Returns (truth_rows, review_rows). Only status=verified records whose PMID
    exists and mentions the gene become truth rows; every problem is listed in
    review_rows for manual checking.
    """
    cand = {c["hgnc_symbol"]: c for c in candidates}
    missing = set(cand) - {c["hgnc_symbol"] for c in curation}
    review = [{"hgnc_symbol": s, "flag": "not_curated", "detail": ""} for s in sorted(missing)]
    verified = [c for c in curation if c["status"] == "verified"]
    recs = pubmed_records(sorted({str(c["pmid"]) for c in verified}))
    truth = []
    for c in curation:
        sym = c["hgnc_symbol"]
        if c["status"] != "verified":
            review.append({"hgnc_symbol": sym, "flag": c["status"], "detail": c.get("notes", "")})
            continue
        pmid = str(c["pmid"])
        rec = recs[pmid]
        names = [sym] + [a for a in cand[sym]["aliases"].split("|") if a]
        if not rec["exists"]:
            review.append({"hgnc_symbol": sym, "flag": "pmid_not_found", "detail": pmid})
            continue
        year = discovery_year(rec)
        flags = []
        if not mentions_gene(rec, names, tuple(cand[sym]["gene_names"].split("|"))):
            flags.append("gene_not_in_title_abstract")
        if c.get("confidence") != "high":
            flags.append(f"confidence_{c.get('confidence')}")
        earlier = earlier_hits(names, year)
        if earlier:
            flags.append("earlier_hits:" + ",".join(earlier))
        for f in flags:
            review.append({"hgnc_symbol": sym, "flag": f.split(":")[0], "detail": f"{pmid} {year} {f}"})
        truth.append(
            {
                "hgnc_symbol": sym,
                "ensembl_gene_id": cand[sym]["ensembl_gene_id"],
                "disease_category": c["disease_category"],
                "tier": "2" if sym in force_tier2 else str(c["tier"]),
                "also_hsp": cand[sym]["also_hsp"],
                "inheritance": c["inheritance"],
                "discovery_year": str(year),
                "pmid": pmid,
                "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "notes": c.get("phenotype_in_paper", ""),
            }
        )
    return sorted(truth, key=lambda r: r["hgnc_symbol"]), review
