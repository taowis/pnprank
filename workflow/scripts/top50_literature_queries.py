"""Re-run the literature searches behind results/tables/top50_literature.tsv.

For the top-N candidates of the released ranking (known tier-1/2 genes excluded, as in
`pnprank top`), this script records:

- PubMed hit counts for five fixed queries per gene (Q1-Q5 below), built from the HGNC
  approved symbol plus HGNC alias and previous symbols (HGNC REST API, fetched live);
- Genomics England PanelApp entries for the gene (all panels, with confidence level);
- ClinGen gene-disease validity curations for the gene (live CSV download).

The per-gene classification (A-E) is manual and is documented in
docs/case_studies/top50_literature.md; this script only reproduces the searches.
PubMed counts change as records are added, so re-runs on a later date will differ a little.

Queries (all fields [tiab]; G = symbol OR aliases; aliases shorter than 4
characters, non-ASCII, generic words, or another gene's approved symbol are dropped
because they match unrelated records):
  Q1  symbol AND NEURO
  Q2  G AND NEURO
  Q3  G AND NEURO AND GENETIC AND PATIENT    (human genetic reports)
  Q4  G AND RELATED AND GENETIC AND PATIENT  (related neurological diseases)
  Q5  N AND PNS AND MODEL                    (functional / model-organism work)
where N = G OR the HGNC approved name OR a common protein name (PROTEIN_NAMES).

Usage:
  python workflow/scripts/top50_literature_queries.py [-n 50] [--dump-dir DIR]

--dump-dir also writes, per gene, the PMIDs and titles returned by Q3, Q4 and Q5 (up to 200
each, PubMed relevance order) for manual review.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import subprocess
import time
import urllib.parse
from pathlib import Path

import pandas as pd

from pnprank.ranking import load_ranking, top_genes

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HGNC = "https://rest.genenames.org/fetch/symbol/"
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"
CLINGEN = "https://search.clinicalgenome.org/kb/gene-validity/download"

NEURO_TERMS = [
    "neuropathy", "neuropathies", "polyneuropathy", "Charcot-Marie-Tooth",
    "hereditary motor neuropathy", "distal hereditary motor neuronopathy",
    "hereditary sensory and autonomic neuropathy", "hereditary sensory neuropathy",
    "spastic paraplegia", "motor neuron disease", "motor neurone disease",
    "spinal muscular atrophy",
]
RELATED_TERMS = [
    "amyotrophic lateral sclerosis", "ataxia", "leukodystrophy", "hypomyelination",
    "hypomyelinating", "myopathy", "muscular dystrophy", "arthrogryposis", "epilepsy",
    "encephalopathy", "neurodevelopmental", "intellectual disability", "dystonia",
    "lissencephaly", "cortical malformation", "neurodegeneration",
]
PNS_TERMS = ["Schwann", "myelin*", "axon*", "sensory neuron*", "motor neuron*",
             "dorsal root ganglion", "peripheral nerve*", "sciatic", "neuromuscular junction"]
MODEL_TERMS = ["knockout", "knock-out", "knockdown", "deficient", "deficiency", "null",
               "mutant", "mutants", "zebrafish", "Drosophila", "conditional"]
# Common protein names not in HGNC (used in Q5 only, together with the HGNC approved name).
PROTEIN_NAMES = {
    "SPP1": ["osteopontin"], "CALCA": ["calcitonin gene-related peptide"],
    "LGALS1": ["galectin-1"], "SNCG": ["gamma-synuclein"], "TARDBP": ["TDP-43"],
    "S100B": ["S100B"], "MBP": ["myelin basic protein"], "PKM": ["pyruvate kinase M2"],
    "ATP5F1B": ["ATP synthase beta"], "KCNA1": ["Kv1.1"],
}
GENETIC_TERMS = ["mutation", "mutations", "variant", "variants", "pathogenic", "de novo",
                 "biallelic", "homozygous", "heterozygous"]
PATIENT_TERMS = ["patient", "patients", "family", "families", "proband", "individuals",
                 "pedigree"]


def _or(terms: list[str]) -> str:
    # Truncated terms (ending in *) cannot be quoted in PubMed.
    return "(" + " OR ".join(f"{t}[tiab]" if t.endswith("*") else f'"{t}"[tiab]'
                             for t in terms) + ")"


def curl(url: str, ua: str | None = None, sleep: float = 0.0) -> bytes:
    if sleep:
        time.sleep(sleep)
    cmd = ["curl", "-sSfL", "-m", "120", "--retry", "4", "--retry-delay", "3"]
    if ua:
        cmd += ["-A", ua]
    return subprocess.run(cmd + [url], check=True, capture_output=True).stdout


def ncbi(path: str, params: dict) -> bytes:
    # 0.4 s between calls keeps us under NCBI's 3 requests/second without an API key.
    return curl(EUTILS + path + "?" + urllib.parse.urlencode(params), sleep=0.4)


def esearch(term: str, retmax: int = 0) -> tuple[int, list[str]]:
    for attempt in range(5):
        d = json.loads(ncbi("esearch.fcgi", {"db": "pubmed", "term": term, "retmode": "json",
                                             "retmax": retmax, "sort": "relevance"}))
        r = d.get("esearchresult", {})
        if "count" in r:  # NCBI occasionally returns an error body with HTTP 200
            return int(r["count"]), r.get("idlist", [])
        time.sleep(2 ** attempt)
    raise RuntimeError(f"esearch failed: {term}\n{d}")


def titles(pmids: list[str]) -> dict[str, dict]:
    out = {}
    for i in range(0, len(pmids), 200):
        d = json.loads(ncbi("esummary.fcgi", {"db": "pubmed", "retmode": "json",
                                              "id": ",".join(pmids[i:i + 200])}))
        for p in d["result"].get("uids", []):
            r = d["result"][p]
            out[p] = {"title": r.get("title", ""), "pubdate": r.get("pubdate", ""),
                      "source": r.get("source", "")}
    return out


# Aliases that are common words or abbreviations rather than gene names.
GENERIC_ALIASES = {"BETA2", "BETA-5"}


def _clean(text: str) -> str:
    return " ".join(str(text).split())


def hgnc_terms(symbol: str, other_symbols: set[str]) -> tuple[list[str], list[str], str]:
    """HGNC alias and previous symbols usable as search terms, and those dropped: shorter
    than 4 characters, non-ASCII, generic words, or the approved symbol of another gene."""
    time.sleep(0.2)
    doc = json.loads(subprocess.run(
        ["curl", "-sSfL", "-m", "60", "-H", "Accept: application/json",
         HGNC + urllib.parse.quote(symbol)], check=True, capture_output=True).stdout)
    docs = doc["response"]["docs"]
    if not docs:
        return [], [], ""
    d = docs[0]
    raw = list(d.get("alias_symbol") or []) + list(d.get("prev_symbol") or [])
    keep = [a for a in raw if len(a) >= 4 and a.isascii() and a.upper() != symbol.upper()
            and a.upper() not in GENERIC_ALIASES and a.upper() not in other_symbols]
    drop = [a for a in raw if a not in keep]
    return keep, drop, d.get("name", "")


def panelapp(symbol: str) -> list[str]:
    out, url = [], PANELAPP + urllib.parse.quote(symbol) + "/"
    while url:
        d = json.loads(curl(url, ua="Mozilla/5.0 (pnprank literature check)", sleep=0.2))
        for r in d["results"]:
            if r.get("entity_name") != symbol:
                continue
            col = {"3": "Green", "2": "Amber", "1": "Red", "0": "Red"}.get(
                str(r.get("confidence_level")), str(r.get("confidence_level")))
            ph = _clean("; ".join(r.get("phenotypes") or []))
            out.append(f"{r['panel']['name']} (id {r['panel']['id']}): {col}"
                       + (f" [{ph}]" if ph else ""))
        url = d.get("next")
    return out


def clingen_table() -> pd.DataFrame:
    text = curl(CLINGEN).decode("utf-8")
    rows = list(csv.reader(io.StringIO(text)))
    hdr = next(i for i, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
    body = [r for r in rows[hdr + 1:] if r and not r[0].startswith("+")]
    return pd.DataFrame(body, columns=rows[hdr])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=50)
    ap.add_argument("--ranking", default="results/tables/pnprank_ranking.tsv.gz")
    ap.add_argument("--out", default="results/tables/top50_literature_queries.tsv")
    ap.add_argument("--sources-out", default="results/tables/top50_literature_curated.tsv")
    ap.add_argument("--dump-dir", default=None)
    a = ap.parse_args()

    today = dt.date.today().isoformat()
    ranking = load_ranking(a.ranking)
    genes = top_genes(ranking, a.n)
    approved = {str(s).upper() for s in ranking.hgnc_symbol}
    clingen = clingen_table()
    neuro, related = _or(NEURO_TERMS), _or(RELATED_TERMS)
    genetic, patient = _or(GENETIC_TERMS), _or(PATIENT_TERMS)
    pns, model = _or(PNS_TERMS), _or(MODEL_TERMS)
    qrows, srows = [], []
    dump = Path(a.dump_dir) if a.dump_dir else None
    if dump:
        dump.mkdir(parents=True, exist_ok=True)

    for _, g in genes.iterrows():
        sym = g.hgnc_symbol
        aliases, dropped, name = hgnc_terms(sym, approved - {sym.upper()})
        gterm = "(" + " OR ".join(f'"{t}"[tiab]' for t in [sym] + aliases) + ")"
        names = [sym] + aliases + [name] + PROTEIN_NAMES.get(sym, [])
        nterm = "(" + " OR ".join(f'"{t}"[tiab]' for t in dict.fromkeys(names) if t) + ")"
        queries = {
            "Q1": f'"{sym}"[tiab] AND {neuro}',
            "Q2": f"{gterm} AND {neuro}",
            "Q3": f"{gterm} AND {neuro} AND {genetic} AND {patient}",
            "Q4": f"{gterm} AND {related} AND {genetic} AND {patient}",
            "Q5": f"{nterm} AND {pns} AND {model}",
        }
        found = {}
        for qid, q in queries.items():
            n, ids = esearch(q, retmax=200 if qid in ("Q3", "Q4", "Q5") else 0)
            found[qid] = ids
            qrows.append({"rank": g["rank"], "hgnc_symbol": sym, "query_id": qid,
                          "hits": n, "date_searched": today, "query": q})
        cg = clingen[clingen["GENE SYMBOL"] == sym]
        srows.append({
            "rank": g["rank"], "hgnc_symbol": sym, "hgnc_name": name,
            "aliases_used": "|".join(aliases), "aliases_dropped": "|".join(dropped),
            "panelapp": " || ".join(panelapp(sym)) or "none",
            "clingen": " || ".join(_clean(f"{r['DISEASE LABEL']} ({r['MOI']}): {r['CLASSIFICATION']}"
                                   f" <{r['ONLINE REPORT']}>") for _, r in cg.iterrows())
                       or "none",
            "date_searched": today,
        })
        if dump:
            ids = sorted(set(found["Q3"]) | set(found["Q4"]) | set(found["Q5"]), key=int)
            meta = titles(ids) if ids else {}
            (dump / f"{sym}.json").write_text(json.dumps(
                {**{k: found[k] for k in ("Q3", "Q4", "Q5")}, "records": meta}, indent=1))
        print(sym, [r["hits"] for r in qrows[-len(queries):]], flush=True)

    pd.DataFrame(qrows).to_csv(a.out, sep="\t", index=False)
    pd.DataFrame(srows).to_csv(a.sources_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
