"""Existing gene-prioritisation approaches, re-implemented on dated data for the
temporal back-test (docs/stage_reports/comparators.md).

Two seed-based comparators use only genes known before the cutoff (`known`) and only
annotation or interaction releases dated before the cutoff:

- `go_similarity_scores`: guilt by association on Gene Ontology, in the style of
  ToppGene/Endeavour functional-annotation ranking. Each gene is a vector of its GO
  terms (propagated to ancestors with the ontology of the same era, weighted by
  inverse document frequency); the score is the cosine similarity to the centroid of
  the known genes' vectors.
- `rwr_scores`: network propagation (random walk with restart) from the known genes on
  a BioGRID physical-interaction release.

Neither uses disease annotations of any kind. Genes with no annotation (GO) or not in
the network (PPI) cannot be ranked by the method and get score 0, i.e. they tie at the
bottom; coverage is reported alongside every result.
"""

from __future__ import annotations

import gzip
import io
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# --- Gene Ontology ------------------------------------------------------------------


def parse_obo(path: str | Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse an OBO file. Returns (parents, alt) where parents[term] are its is_a and
    part_of parents and alt maps alt_id -> primary id. Obsolete terms are left out."""
    parents: dict[str, set[str]] = {}
    alt: dict[str, str] = {}
    term = None
    in_term = False
    obsolete = False
    cur_par: set[str] = set()
    cur_alt: list[str] = []

    def close():
        if term is not None and not obsolete:
            parents[term] = cur_par
            for a in cur_alt:
                alt[a] = term

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("["):
                close()
                term, obsolete, cur_par, cur_alt = None, False, set(), []
                in_term = line == "[Term]"
                continue
            if not line or ":" not in line:
                continue
            key, val = line.split(":", 1)
            val = val.split("!")[0].strip()
            if key == "id" and in_term and term is None:
                term = val
            elif term is None:
                continue
            elif key == "is_a":
                cur_par.add(val.split()[0])
            elif key == "relationship" and val.startswith("part_of "):
                cur_par.add(val.split()[1])
            elif key == "alt_id":
                cur_alt.append(val)
            elif key == "is_obsolete" and val == "true":
                obsolete = True
        close()
    return parents, alt


def ancestors_closure(parents: dict[str, set[str]]) -> dict[str, frozenset[str]]:
    """Each term mapped to itself plus all its ancestors (is_a / part_of)."""
    memo: dict[str, frozenset[str]] = {}

    def anc(t, stack=()):
        if t in memo:
            return memo[t]
        out = {t}
        for p in parents.get(t, ()):
            if p in parents and p not in stack:
                out |= anc(p, stack + (t,))
        memo[t] = frozenset(out)
        return memo[t]

    for t in parents:
        anc(t)
    return memo


def read_gaf(path: str | Path, before: str | None = None,
             exclude_evidence: tuple[str, ...] = ("ND",)) -> pd.DataFrame:
    """Read a GO annotation file (GAF 1.0/2.x). NOT-qualified annotations and the given
    evidence codes are dropped. With `before` ("YYYYMMDD"), annotations dated on or after
    it are dropped as well (a guard; a release dated before the cutoff should have none).
    Returns columns db_id, symbol, go_id, evidence, date."""
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 14 or "NOT" in f[3].split("|"):
                continue
            if f[6] in exclude_evidence:
                continue
            rows.append((f[1], f[2], f[4], f[6], f[13]))
    df = pd.DataFrame(rows, columns=["db_id", "symbol", "go_id", "evidence", "date"])
    if before is not None:
        df = df[df.date.str[:8] < before]
    return df


def map_to_ensembl(ids: pd.Series, symbols: pd.Series, uniprot: dict[str, str],
                   sym: dict[str, str]) -> pd.Series:
    """Ensembl gene ID per row: UniProt accession via HGNC first, then the symbol via
    HGNC approved / unique previous / unique alias symbols (`atlases.symbol_to_ensembl`)."""
    by_acc = ids.str.split("-").str[0].map(uniprot)
    return by_acc.fillna(symbols.map(sym))


def hgnc_uniprot(hgnc_tsv: str | Path) -> dict[str, str]:
    """UniProt accession -> Ensembl gene ID (accessions shared by several genes dropped)."""
    h = pd.read_csv(hgnc_tsv, sep="\t", dtype=str, usecols=["uniprot_ids", "ensembl_gene_id", "status"])
    h = h[(h.status == "Approved") & h.ensembl_gene_id.notna() & h.uniprot_ids.notna()]
    m: dict[str, set[str]] = defaultdict(set)
    for acc, ens in zip(h.uniprot_ids, h.ensembl_gene_id):
        for a in acc.split("|"):
            m[a].add(ens)
    return {a: next(iter(e)) for a, e in m.items() if len(e) == 1}


def uniprot_to_ensembl(hgnc_tsv: str | Path, idmapping: str | Path | None = None,
                       sec_ac: str | Path | None = None) -> dict[str, str]:
    """UniProt accession -> Ensembl gene ID.

    1. HGNC's own UniProt cross-references (`hgnc_uniprot`).
    2. UniProt's human idmapping (accession -> HGNC ID -> Ensembl via HGNC), which also
       covers unreviewed (TrEMBL) entries.
    3. Secondary accessions (UniProt sec_ac.txt) resolved to their primary accession,
       then 1-2. Only unambiguous mappings are kept at each step.
    """
    out = hgnc_uniprot(hgnc_tsv)
    if idmapping is not None:
        hg = hgnc_column_map(hgnc_tsv, "hgnc_id")
        m: dict[str, set[str]] = defaultdict(set)
        with gzip.open(idmapping, "rt") as fh:
            for line in fh:
                acc, kind, val = line.rstrip("\n").split("\t")
                if kind == "HGNC" and val in hg:
                    m[acc.split("-")[0]].add(hg[val])
        for acc, ens in m.items():
            if acc not in out and len(ens) == 1:
                out[acc] = next(iter(ens))
    if sec_ac is not None:
        prim: dict[str, set[str]] = defaultdict(set)
        started = False
        with open(sec_ac) as fh:
            for line in fh:
                if line.startswith("_____"):
                    started = True
                    continue
                f = line.split()
                if started and len(f) == 2:
                    prim[f[0]].add(f[1])
        for sec, ps in prim.items():
            ens = {out[p] for p in ps if p in out}
            if sec not in out and len(ens) == 1:
                out[sec] = next(iter(ens))
    return out


def go_gene_matrix(gene_terms: pd.DataFrame, parents: dict[str, set[str]], alt: dict[str, str],
                   universe: list[str], min_genes: int = 3, max_genes: int = 500
                   ) -> tuple[sp.csr_matrix, list[str]]:
    """Binary gene x term matrix (rows follow `universe`), with annotations propagated to
    ancestors. `gene_terms` has columns ensembl_gene_id, go_id. Terms annotated to fewer
    than `min_genes` or more than `max_genes` universe genes are dropped."""
    closure = ancestors_closure(parents)
    pos = {g: i for i, g in enumerate(universe)}
    pairs = set()
    for g, t in zip(gene_terms.ensembl_gene_id, gene_terms.go_id):
        if g not in pos:
            continue
        t = alt.get(t, t)
        for a in closure.get(t, ()):
            pairs.add((pos[g], a))
    terms = sorted({t for _, t in pairs})
    tpos = {t: j for j, t in enumerate(terms)}
    r = np.fromiter((i for i, _ in pairs), dtype=np.int64, count=len(pairs))
    c = np.fromiter((tpos[t] for _, t in pairs), dtype=np.int64, count=len(pairs))
    m = sp.csr_matrix((np.ones(len(pairs)), (r, c)), shape=(len(universe), len(terms)))
    n = np.asarray(m.sum(axis=0)).ravel()
    keep = np.flatnonzero((n >= min_genes) & (n <= max_genes))
    return m[:, keep].tocsr(), [terms[j] for j in keep]


def go_similarity_scores(m: sp.csr_matrix, universe: list[str], known: set[str]) -> pd.Series:
    """Cosine similarity of each gene's IDF-weighted GO vector to the centroid of the
    known genes' (unit-normalised) vectors. Unannotated genes, and every gene when no
    known gene is annotated, score 0. Depends on `known` only through the centroid."""
    n_genes = np.asarray(m.sum(axis=1) > 0).ravel().sum()
    df = np.asarray(m.sum(axis=0)).ravel()
    idf = np.log(max(n_genes, 1) / np.maximum(df, 1))
    x = m.multiply(idf[None, :]).tocsr()
    norm = np.sqrt(np.asarray(x.multiply(x).sum(axis=1)).ravel())
    x = sp.diags(np.where(norm > 0, 1 / np.where(norm > 0, norm, 1), 0)) @ x
    idx = [i for i, g in enumerate(universe) if g in known and norm[i] > 0]
    if not idx:
        return pd.Series(0.0, index=universe)
    centroid = np.asarray(x[idx].mean(axis=0)).ravel()
    return pd.Series(x @ centroid, index=universe)


# --- Protein interactions -----------------------------------------------------------


# Experimental systems BioGRID classes as "genetic" (Experimental System Type in the
# 3.x tab2 files). Release 2.0 tab files have no type column; their systems are
# classified with this list, and the legacy "Invivo"/"Invitro" systems (protein
# interactions) count as physical. See docs/decisions.md.
GENETIC_SYSTEMS = frozenset({
    "Dosage Growth Defect", "Dosage Lethality", "Dosage Rescue", "Negative Genetic",
    "Phenotypic Enhancement", "Phenotypic Suppression", "Positive Genetic",
    "Synthetic Growth Defect", "Synthetic Haploinsufficiency", "Synthetic Lethality",
    "Synthetic Rescue", "Epistatic MAP",
})


def read_biogrid(path: str | Path, entrez: dict[str, str], sym: dict[str, str],
                 physical_only: bool = True, hgnc_ids: dict[str, str] | None = None
                 ) -> pd.DataFrame:
    """Human-human interactions from a BioGRID ORGANISM tab (2.0) or tab2 (>= 3.0) zip.

    Returns unique undirected Ensembl pairs (a < b), self-interactions removed. tab2 rows
    are mapped by Entrez Gene ID via HGNC, tab (2.0, no Entrez column) by HGNC ID
    (`hgnc_ids`); either falls back to the official symbol via HGNC approved / unique
    previous / unique alias symbols."""
    with zipfile.ZipFile(path) as z:
        name = next(n for n in z.namelist() if "Homo_sapiens" in n)
        raw = z.read(name).decode("utf-8", errors="replace")
    lines = raw.splitlines()
    hdr = next(i for i, ln in enumerate(lines) if ln.startswith(("#", "INTERACTOR_A")))
    df = pd.read_csv(io.StringIO("\n".join(lines[hdr:])), sep="\t", dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    if "Entrez Gene Interactor A" in df.columns:  # tab2
        org = df["Organism Interactor A"].eq("9606") & df["Organism Interactor B"].eq("9606")
        phys = df["Experimental System Type"].str.lower().eq("physical")
        a = df["Entrez Gene Interactor A"].map(entrez).fillna(df["Official Symbol Interactor A"].map(sym))
        b = df["Entrez Gene Interactor B"].map(entrez).fillna(df["Official Symbol Interactor B"].map(sym))
    else:  # tab 2.0: no system-type column, HGNC IDs as interactor IDs
        org = df["ORGANISM_A_ID"].eq("9606") & df["ORGANISM_B_ID"].eq("9606")
        phys = ~df["EXPERIMENTAL_SYSTEM"].isin(GENETIC_SYSTEMS)
        hg = hgnc_ids or {}
        a = df["INTERACTOR_A"].map(hg).fillna(df["OFFICIAL_SYMBOL_A"].map(sym))
        b = df["INTERACTOR_B"].map(hg).fillna(df["OFFICIAL_SYMBOL_B"].map(sym))
    sel = org & (phys if physical_only else True) & a.notna() & b.notna() & (a != b)
    e = pd.DataFrame({"a": a[sel], "b": b[sel]})
    lo, hi = e.min(axis=1), e.max(axis=1)
    return pd.DataFrame({"a": lo, "b": hi}).drop_duplicates().reset_index(drop=True)


def hgnc_column_map(hgnc_tsv: str | Path, column: str) -> dict[str, str]:
    """`column` (e.g. entrez_id, hgnc_id) -> Ensembl gene ID for approved HGNC genes."""
    h = pd.read_csv(hgnc_tsv, sep="\t", dtype=str, usecols=[column, "ensembl_gene_id", "status"])
    h = h[(h.status == "Approved") & h.ensembl_gene_id.notna() & h[column].notna()]
    return dict(zip(h[column], h.ensembl_gene_id))


def adjacency(edges: pd.DataFrame, universe: list[str]) -> sp.csr_matrix:
    """Symmetric binary adjacency over `universe` (edges outside it are dropped)."""
    pos = {g: i for i, g in enumerate(universe)}
    e = edges[edges.a.isin(pos) & edges.b.isin(pos)]
    r = e.a.map(pos).to_numpy()
    c = e.b.map(pos).to_numpy()
    n = len(universe)
    m = sp.csr_matrix((np.ones(len(e)), (r, c)), shape=(n, n))
    m = ((m + m.T) > 0).astype(float)
    return m.tocsr()


def rwr_scores(adj: sp.csr_matrix, universe: list[str], known: set[str], restart: float = 0.5,
               tol: float = 1e-10, max_iter: int = 1000) -> pd.Series:
    """Random walk with restart from the known genes (uniform restart vector over those
    in the network), on the column-normalised adjacency. Genes outside the network and
    isolated genes score 0. Depends on `known` only through the restart vector."""
    deg = np.asarray(adj.sum(axis=0)).ravel()
    w = adj @ sp.diags(np.where(deg > 0, 1 / np.where(deg > 0, deg, 1), 0))
    seeds = [i for i, g in enumerate(universe) if g in known and deg[i] > 0]
    if not seeds:
        return pd.Series(0.0, index=universe)
    p0 = np.zeros(len(universe))
    p0[seeds] = 1 / len(seeds)
    p = p0.copy()
    for _ in range(max_iter):
        nxt = (1 - restart) * (w @ p) + restart * p0
        if np.abs(nxt - p).sum() < tol:
            p = nxt
            break
        p = nxt
    return pd.Series(p, index=universe)


# --- back-test wrapper --------------------------------------------------------------


def seed_scores_at_cutoff(score_fn, y: pd.Series) -> pd.Series:
    """Run a seed-based comparator with the training labels of a cutoff.

    `y` comes from `pnprank.backtest.training_labels`: the seeds are the genes labelled 1
    (known before the cutoff) and nothing else; genes labelled 0 are scored, positives
    are not (as for the PNPrank back-test scores)."""
    known = set(y.index[y == 1])
    return score_fn(known).reindex(y.index[y == 0])
