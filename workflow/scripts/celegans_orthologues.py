"""C. elegans orthologues of top-ranked candidates and of tier-1 neuropathy genes.

Outputs
- results/tables/celegans_shortlist.tsv: candidates from the top 200 unlabelled genes
  (`pnprank top -n 200`, known tier-1/2 genes excluded) that have a clear worm
  orthologue, with ranking columns, orthology evidence and a testability flag.
- results/tables/tier1_celegans_orthologues.tsv: the 100 tier-1 genes and their best
  worm orthologue (if any).

Sources (queried live, responses cached as JSON under <data_root>/external/celegans/)
- HGNC REST (Ensembl gene ID -> HGNC ID, needed by the Alliance API):
  https://rest.genenames.org/fetch/ensembl_gene_id/{ENSG}
- Ensembl REST, pan-taxonomic Compara (the default vertebrate Compara of release 116
  returns no C. elegans homologies):
  https://rest.ensembl.org/homology/id/human/{ENSG}?target_species=caenorhabditis_elegans
  &type=orthologues&compara=pan_homology&sequence=none
- Alliance of Genome Resources orthology (DIOPT-style method counts):
  https://www.alliancegenome.org/api/gene/{HGNC_ID}/orthologs?filter.stringency=all
- Optional enrichment of the chosen worm gene, from Alliance only:
  alleles  https://www.alliancegenome.org/api/gene/{WB_ID}/alleles?filter.alleleCategory=allele
  expression ribbon (POST) https://www.alliancegenome.org/api/expression/ribbon-summary

"Clear orthologue" rule (docs/decisions.md, 2026-10-08): a human-worm pair is clear if
(a) Ensembl pan-Compara calls it ortholog_one2one, or
(b) Alliance reports it as best score in both directions and at least
    MIN_METHODS of the prediction methods support it.

Usage: python workflow/scripts/celegans_orthologues.py
"""

from __future__ import annotations

import json
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

from pnprank.ranking import load_ranking, top_genes

N_TOP = 200
MIN_METHODS = 5          # of the 9 methods Alliance runs for human-worm pairs
MIN_INTERVAL_S = 0.25    # at most 4 requests per second
NERVOUS_SYSTEM = "UBERON:0001016"
ALL_ANATOMY = "UBERON:0001062"

ENSEMBL = "https://rest.ensembl.org"
ALLIANCE = "https://www.alliancegenome.org/api"
HGNC = "https://rest.genenames.org"

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
cache = Path(paths["data_root"]) / "external" / "celegans"
out_dir = Path("results/tables")

_lock = threading.Lock()
_last = [0.0]


def _request(url: str, body: bytes | None = None) -> object:
    """GET (or POST `body`) returning parsed JSON; throttled, retried with backoff."""
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(6):
        with _lock:
            wait = MIN_INTERVAL_S - (time.monotonic() - _last[0])
            if wait > 0:
                time.sleep(wait)
            _last[0] = time.monotonic()
        try:
            req = urllib.request.Request(url, data=body, headers=headers,
                                         method="POST" if body else "GET")
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                raise
            retry = e.headers.get("Retry-After")
            time.sleep(float(retry) if retry else 2 ** attempt)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed after retries: {url}")


def cached(kind: str, key: str, url: str, body: bytes | None = None) -> object:
    f = cache / kind / f"{key.replace(':', '_')}.json"
    if f.exists():
        return json.loads(f.read_text())
    f.parent.mkdir(parents=True, exist_ok=True)
    data = _request(url, body)
    f.write_text(json.dumps(data))
    return data


def hgnc_id(ensg: str) -> str | None:
    docs = cached("hgnc", ensg, f"{HGNC}/fetch/ensembl_gene_id/{ensg}")["response"]["docs"]
    return docs[0]["hgnc_id"] if len(docs) == 1 else None


def ensembl_pairs(ensg: str) -> pd.DataFrame:
    d = cached("ensembl", ensg,
               f"{ENSEMBL}/homology/id/human/{ensg}?target_species=caenorhabditis_elegans"
               "&type=orthologues&compara=pan_homology&sequence=none")
    rows = [{"wb_id": h["target"]["id"], "ens_type": h["type"].replace("ortholog_", ""),
             "ens_pct_id_human": round(h["source"]["perc_id"], 1),
             "ens_pct_id_worm": round(h["target"]["perc_id"], 1)}
            for x in d["data"] for h in x["homologies"]
            if h["target"]["species"] == "caenorhabditis_elegans"]
    return pd.DataFrame(rows, columns=["wb_id", "ens_type", "ens_pct_id_human",
                                       "ens_pct_id_worm"]).drop_duplicates("wb_id")


def alliance_pairs(hgnc: str | None) -> pd.DataFrame:
    cols = ["wb_id", "worm_gene", "all_n_methods", "all_n_methods_called", "all_best",
            "all_best_reverse", "all_methods"]
    if hgnc is None:
        return pd.DataFrame(columns=cols)
    d = cached("alliance_orthologs", hgnc,
               f"{ALLIANCE}/gene/{hgnc}/orthologs?filter.stringency=all&limit=1000")
    if d.get("total", 0) > d.get("returnedRecords", 0):
        raise RuntimeError(f"Alliance orthologs paginated for {hgnc}")
    rows = []
    for r in d["results"]:
        g = r["geneToGeneOrthologyGenerated"]
        o = g["objectGene"]
        if o["taxon"]["curie"] != "NCBITaxon:6239":
            continue
        matched = [m["name"] for m in g.get("predictionMethodsMatched", [])]
        not_matched = g.get("predictionMethodsNotMatched", [])
        rows.append({"wb_id": o["primaryExternalId"].removeprefix("WB:"),
                     "worm_gene": o["geneSymbol"]["displayText"],
                     "all_n_methods": len(matched),
                     "all_n_methods_called": len(matched) + len(not_matched),
                     "all_best": g["isBestScore"]["name"] == "Yes",
                     "all_best_reverse": g["isBestScoreReverse"]["name"] == "Yes",
                     "all_methods": ",".join(sorted(matched))})
    return pd.DataFrame(rows, columns=cols)


def orthologues(ensg: str) -> tuple[dict, list[str]]:
    """Best worm orthologue of one human gene, and a list of the other worm genes."""
    h = hgnc_id(ensg)
    ens, al = ensembl_pairs(ensg), alliance_pairs(h)
    m = ens.merge(al, on="wb_id", how="outer")
    base = {"hgnc_id": h, "n_worm_orthologues_any": len(m)}
    if m.empty:
        return base | {"clear_orthologue": False, "orthology_class": "none"}, []
    m["all_n_methods"] = pd.to_numeric(m.all_n_methods).fillna(0).astype(int)
    m["all_best"] = m.all_best.astype("boolean").fillna(False).astype(bool)
    m["all_best_reverse"] = m.all_best_reverse.astype("boolean").fillna(False).astype(bool)
    m["rule_ens_one2one"] = m.ens_type.eq("one2one")
    m["rule_alliance"] = m.all_best & m.all_best_reverse & (m.all_n_methods >= MIN_METHODS)
    m["clear"] = m.rule_ens_one2one | m.rule_alliance
    m = m.sort_values(["clear", "all_n_methods", "ens_pct_id_human"],
                      ascending=False, na_position="last")
    best = m.iloc[0].to_dict()
    others = [f"{(r.worm_gene if isinstance(r.worm_gene, str) else r.wb_id)}"
              f"({r.all_n_methods})" for r in m.iloc[1:].itertuples()
              if r.all_n_methods >= 3]
    rule = "+".join(n for n, k in [("ensembl_one2one", "rule_ens_one2one"),
                                    ("alliance_reciprocal_best", "rule_alliance")]
                    if best[k])
    return base | {
        "worm_gene": best.get("worm_gene") if isinstance(best.get("worm_gene"), str)
        else None,
        "wormbase_id": best["wb_id"],
        "ensembl_type": best.get("ens_type"),
        "ensembl_pct_id_human": best.get("ens_pct_id_human"),
        "ensembl_pct_id_worm": best.get("ens_pct_id_worm"),
        "alliance_n_methods": int(best["all_n_methods"]),
        "alliance_n_methods_called": best.get("all_n_methods_called"),
        "alliance_best_both_ways": bool(best["all_best"] and best["all_best_reverse"]),
        "alliance_methods": best.get("all_methods"),
        "clear_orthologue": bool(best["clear"]),
        "clear_rule": rule or None,
        "orthology_class": orthology_class(bool(best["clear"]), int(best["all_n_methods"])),
        "n_worm_genes_clear": int(m.clear.sum()),
    }, others


def orthology_class(clear: bool, n_methods: int) -> str:
    """clear: the rule above. supported: >= MIN_METHODS Alliance methods but not reciprocal
    best (typically a worm gene shared with a human paralogue that is the better match).
    weak: any other call."""
    if clear:
        return "clear"
    return "supported" if n_methods >= MIN_METHODS else "weak"


def human_genes_of_worm_gene(wb: str) -> list[str]:
    """Human genes that Alliance links to this worm gene with >= MIN_METHODS methods."""
    d = cached("alliance_orthologs", f"WB_{wb}",
               f"{ALLIANCE}/gene/WB:{wb}/orthologs?filter.stringency=all&limit=1000")
    if d.get("total", 0) > d.get("returnedRecords", 0):
        raise RuntimeError(f"Alliance orthologs paginated for {wb}")
    out = []
    for r in d["results"]:
        g = r["geneToGeneOrthologyGenerated"]
        o = g["objectGene"]
        if (o["taxon"]["curie"] == "NCBITaxon:9606"
                and len(g.get("predictionMethodsMatched", [])) >= MIN_METHODS):
            out.append(o["geneSymbol"]["displayText"])
    return sorted(set(out))


def worm_gene_name(wb: str) -> str | None:
    d = cached("alliance_gene", wb, f"{ALLIANCE}/gene/WB:{wb}")
    return d.get("gene", {}).get("geneSymbol", {}).get("displayText")


def enrichment(wb: str) -> dict:
    a = cached("alliance_alleles", wb,
               f"{ALLIANCE}/gene/WB:{wb}/alleles?limit=1&filter.alleleCategory=allele")
    e = cached("alliance_expression", wb, f"{ALLIANCE}/expression/ribbon-summary",
               body=json.dumps([f"WB:{wb}"]).encode())
    groups = e["subjects"][0]["groups"] if e.get("subjects") else {}
    n = lambda k: groups.get(k, {}).get("ALL", {}).get("nb_annotations", 0)  # noqa: E731
    return {"worm_n_alleles": a.get("total", 0),
            "worm_expr_annotations": n(ALL_ANATOMY),
            "worm_expr_nervous_system": n(NERVOUS_SYSTEM)}


def cell_class(top: str) -> str:
    ct = top.split("|")[1].split(" (")[0] if isinstance(top, str) and "|" in top else ""
    if ct.startswith("neuron_"):
        return "neuronal"
    if ct.startswith("schwann") or ct == "oligodendrocyte":
        return "myelinating_glia"
    if ct in ("satellite_glia", "astrocyte", "microglia", "opc"):
        return "other_glia"
    return "non_neural" if ct else "unknown"


def annotate(genes: pd.DataFrame, enrich: bool) -> pd.DataFrame:
    rows = []
    for i, g in enumerate(genes.itertuples(), 1):
        o, others = orthologues(g.ensembl_gene_id)
        if o.get("wormbase_id") and not o.get("worm_gene"):
            o["worm_gene"] = worm_gene_name(o["wormbase_id"])
        o["other_worm_orthologues_ge3_methods"] = ";".join(others) or None
        if o.get("wormbase_id") and o["orthology_class"] != "weak":
            hs = human_genes_of_worm_gene(o["wormbase_id"])
            o["human_genes_sharing_worm_gene"] = ";".join(
                h for h in hs if h != g.hgnc_symbol) or None
        if enrich and o.get("clear_orthologue"):
            o |= enrichment(o["wormbase_id"])
        rows.append({"ensembl_gene_id": g.ensembl_gene_id} | o)
        if i % 25 == 0:
            print(f"  {i}/{len(genes)}", flush=True)
    return genes.merge(pd.DataFrame(rows), on="ensembl_gene_id", how="left")


def main() -> None:
    release = cached("info", "ensembl_release", f"{ENSEMBL}/info/data")["releases"][0]
    eg = cached("info", "eg_version", f"{ENSEMBL}/info/eg_version")["version"]
    al = cached("info", "alliance_release", f"{ALLIANCE}/releaseInfo")["releaseVersion"]
    print(f"Ensembl REST release {release} (EG version {eg}); Alliance {al}")

    tab = load_ranking("results/tables/pnprank_ranking.tsv.gz")
    cand = top_genes(tab, N_TOP)
    assert len(cand) == N_TOP
    keep = ["rank", "hgnc_symbol", "ensembl_gene_id", "status", "percentile",
            "percentile_expr_matched", "percentile_non_atlas", "gtex_mean_decile",
            "top_cell_type"]
    print(f"Candidates: top {N_TOP} unlabelled genes")
    cand = annotate(cand[keep].reset_index(drop=True), enrich=True)
    cand["cell_class"] = cand.top_cell_type.map(cell_class)
    cand["worm_testable"] = cand.clear_orthologue & cand.cell_class.ne("myelinating_glia")
    short = cand[cand.clear_orthologue.fillna(False).astype(bool)].copy()
    short["note"] = short.cell_class.map(
        {"myelinating_glia": "top cell type is Schwann/myelinating glia: less testable "
                             "in worms (no myelin)"}).fillna("")
    short.to_csv(out_dir / "celegans_shortlist.tsv", sep="\t", index=False)
    cand.to_csv(cache / "top200_all_orthology.tsv", sep="\t", index=False)

    truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
    t1 = truth[truth.tier == "1"][["hgnc_symbol", "ensembl_gene_id", "disease_category"]]
    print(f"Tier-1 genes: {len(t1)}")
    t1 = annotate(t1.reset_index(drop=True), enrich=False)
    t1_cols = ["hgnc_symbol", "ensembl_gene_id", "disease_category", "hgnc_id",
               "worm_gene", "wormbase_id", "ensembl_type", "ensembl_pct_id_human",
               "alliance_n_methods", "alliance_best_both_ways", "clear_orthologue",
               "clear_rule", "orthology_class", "human_genes_sharing_worm_gene",
               "n_worm_orthologues_any"]
    t1[t1_cols].to_csv(out_dir / "tier1_celegans_orthologues.tsv", sep="\t", index=False)

    any1 = (t1.n_worm_orthologues_any > 0).sum()
    clear1 = t1.clear_orthologue.fillna(False).astype(bool).sum()
    print(f"Tier-1 with any worm orthologue call: {any1}/{len(t1)}; clear: {clear1}")
    print("Tier-1 orthology classes:", t1.orthology_class.value_counts().to_dict())
    print("Top-200 orthology classes:", cand.orthology_class.value_counts().to_dict())
    print(f"Top {N_TOP}: clear orthologue {len(short)}; "
          f"of which myelinating-glia top cell type {(short.cell_class == 'myelinating_glia').sum()}")


if __name__ == "__main__":
    main()
