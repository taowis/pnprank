"""Comparators inside the Stage 5 temporal back-test (docs/stage_reports/comparators.md).

At each cutoff and positive set, the PNPrank back-test scores cached by
stage5_backtest.py are compared with:
- go_similarity: GO guilt-by-association on the GOA human release and GO ontology dated
  before the cutoff, seeded with the genes known before the cutoff;
- ppi_rwr: random walk with restart on the BioGRID physical-interaction release dated
  before the cutoff (no release predates 2005, so not run at 2005), same seeds;
- label-free gnomAD constraint scores (pLI, missense Z) next to the Stage 5 baselines.

Same evaluation labels, transforms (within-fold percentile; fold x expression decile)
and paired bootstrap as Stage 5.

Usage: python workflow/scripts/comparators_backtest.py [scores_dir] [n_jobs]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score

from pnprank.atlases import symbol_to_ensembl
from pnprank.backtest import (
    auroc,
    bootstrap_ci,
    evaluation_labels,
    paired_bootstrap_diff,
    training_labels,
    within_fold_percentile,
)
from pnprank.comparators import (
    adjacency,
    go_gene_matrix,
    go_similarity_scores,
    hgnc_column_map,
    map_to_ensembl,
    parse_obo,
    read_biogrid,
    read_gaf,
    rwr_scores,
    seed_scores_at_cutoff,
    uniprot_to_ensembl,
)

scores_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/backtest/scores")
n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 1
cfg = yaml.safe_load(Path("config/config.yaml").read_text())
ccfg = yaml.safe_load(Path("config/comparators.yaml").read_text())
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
root = Path(paths["data_root"])
cdir = root / "external" / "comparators"
bc, meth = cfg["backtest"], ccfg["methods"]
seed, nboot, ci = cfg["seeds"]["global"], cfg["bootstrap"]["n_resamples"], cfg["bootstrap"]["ci"]
out = Path("results/backtest/comparators")
out.mkdir(parents=True, exist_ok=True)

X_static = pd.read_csv(root / "processed" / "features_static.tsv.gz", sep="\t", index_col=0)
universe = list(X_static.index)
truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
negatives = set(universe) - set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id)
expr = X_static["gtex|mean_all|logtpm"]
TIERS = {"tier1": {"1"}, "tier1+2": {"1", "2"}}

hgnc = root / "external" / "truth_sources" / "hgnc_complete_set.txt"
sym, _ = symbol_to_ensembl(hgnc)
idm = ccfg["id_mapping"]
udir = cdir / "uniprot" / idm["release"]
uniprot = uniprot_to_ensembl(hgnc, udir / Path(idm["human_idmapping"]).name,
                             udir / Path(idm["sec_ac"]).name)
entrez = hgnc_column_map(hgnc, "entrez_id")
hgnc_ids = hgnc_column_map(hgnc, "hgnc_id")

# --- dated data per cutoff ----------------------------------------------------------
data, coverage = {}, []
for c in bc["cutoff_years"]:
    ga, go = ccfg["go_annotation"]["releases"][c], ccfg["go_ontology"]["releases"][c]
    gaf = read_gaf(cdir / "goa" / Path(ga["url"]).name, before=f"{c}0101",
                   exclude_evidence=tuple(meth["go_exclude_evidence"]))
    gaf_all = read_gaf(cdir / "goa" / Path(ga["url"]).name,
                       exclude_evidence=tuple(meth["go_exclude_evidence"]))
    gaf["ensembl_gene_id"] = map_to_ensembl(gaf.db_id, gaf.symbol, uniprot, sym)
    parents, alt = parse_obo(cdir / "go" / go["release"] / Path(go["url"]).name)
    n_terms_known = gaf.go_id.map(lambda t: alt.get(t, t)).isin(parents).mean()
    m, terms = go_gene_matrix(gaf.dropna(subset=["ensembl_gene_id"]), parents, alt, universe,
                              meth["go_min_genes"], meth["go_max_genes"])
    rec = {"cutoff": c, "go_release": ga["release"], "go_release_date": ga["date"],
           "go_ontology": go["release"], "gaf_rows": len(gaf_all),
           "gaf_rows_dated_after_cutoff": len(gaf_all) - len(gaf),
           "gaf_proteins": gaf.db_id.nunique(),
           "gaf_proteins_mapped": gaf.dropna(subset=["ensembl_gene_id"]).db_id.nunique(),
           "gaf_terms_in_ontology": round(float(n_terms_known), 4),
           "go_terms_used": len(terms),
           "go_genes_annotated": int((m.sum(axis=1) > 0).sum())}
    adj = None
    if c in ccfg["interactions"]["releases"]:
        r = ccfg["interactions"]["releases"][c]
        edges = read_biogrid(cdir / "biogrid" / r["file"], entrez, sym,
                             meth["ppi_physical_only"], hgnc_ids)
        adj = adjacency(edges, universe)
        rec |= {"ppi_release": r["release"], "ppi_release_date": r["date"],
                "ppi_edges": int(adj.nnz // 2),
                "ppi_genes": int((np.asarray(adj.sum(axis=0)).ravel() > 0).sum())}
    coverage.append(rec)
    data[c] = (m, adj)
    print(rec, flush=True)
pd.DataFrame(coverage).to_csv(out / "comparator_data_coverage.tsv", sep="\t", index=False)

LABEL_FREE = {"neg_loeuf": -X_static["gnomad|loeuf"], "pli": X_static["gnomad|pli"],
              "mis_z": X_static["gnomad|mis_z"]}
# go_term_count and ppi_degree are seed-free study-bias controls: how well does "has
# many GO annotations" or "has many interaction partners" alone rank future genes?
SCORES = ["full", "non_atlas", "gtex_nerve", "go_similarity", "ppi_rwr", "neg_loeuf", "pli",
          "mis_z", "go_term_count", "ppi_degree"]
COMPARATORS = [s for s in SCORES if s != "full"]


def scores_at(set_name: str, c: int) -> pd.DataFrame:
    df = pd.read_csv(scores_dir / f"{set_name}_{c}.tsv.gz", sep="\t", index_col=0)
    y = training_labels(universe, truth, c, TIERS[set_name])
    zeros = y.index[y == 0]
    # the cached scores must be exactly the genes labelled 0 under these labels
    assert set(df.index) == set(zeros), (set_name, c)
    m, adj = data[c]
    df["go_similarity"] = seed_scores_at_cutoff(
        lambda k: go_similarity_scores(m, universe, k), y).reindex(df.index)
    df["ppi_rwr"] = (seed_scores_at_cutoff(
        lambda k: rwr_scores(adj, universe, k, meth["rwr_restart"]), y).reindex(df.index)
        if adj is not None else np.nan)
    for k, v in LABEL_FREE.items():
        df[k] = v.reindex(df.index)
    df["go_term_count"] = pd.Series(np.asarray(m.sum(axis=1)).ravel(), index=universe
                                    ).reindex(df.index)
    df["ppi_degree"] = (pd.Series(np.asarray(adj.sum(axis=0)).ravel(), index=universe
                                  ).reindex(df.index) if adj is not None else np.nan)
    n_seed_go = len([g for g in y.index[y == 1]
                     if m[universe.index(g)].sum() > 0])
    n_seed_ppi = (len([g for g in y.index[y == 1] if adj[universe.index(g)].sum() > 0])
                  if adj is not None else 0)
    df.attrs["seeds"] = {"n_seeds": int(y.sum()), "seeds_annotated_go": n_seed_go,
                         "seeds_in_ppi": n_seed_ppi}
    return df


def ap(s, y):
    ok = ~np.isnan(s)
    return average_precision_score(y[ok], s[ok])


def evaluate(set_name: str, c: int, k0: int):
    df = scores_at(set_name, c)
    e = expr.reindex(df.index)
    tr = {}
    for s in SCORES:
        tr[("fold_rank", s)] = within_fold_percentile(df[s], df["fold"])
        tr[("fold_rank_expr_matched", s)] = within_fold_percentile(df[s], df["fold"], e,
                                                                   bc["expression_bins"])
    perf, pairs, cov = [], [], []
    k = k0
    for window in (None, bc["interval_window_years"]):
        ev = evaluation_labels(truth, negatives, c, TIERS[set_name], window)
        ev = ev[ev.index.isin(df.index)]
        lab = ev.to_numpy().astype(bool)
        wname = "all_after" if window is None else f"within_{window}y"
        for s in ("go_similarity", "ppi_rwr"):
            v = df[s].reindex(ev.index)
            if v.notna().any():
                cov.append({"positives": set_name, "cutoff": c, "eval": wname, "score": s,
                            **df.attrs["seeds"],
                            "pos_covered": int((v[lab] > 0).sum()), "n_pos": int(lab.sum()),
                            "neg_covered": int((v[~lab] > 0).sum()), "n_neg": int((~lab).sum())})
        for t in ("fold_rank", "fold_rank_expr_matched"):
            for s in SCORES:
                sc = tr[(t, s)].reindex(ev.index).to_numpy(float)
                if np.isnan(sc).all():
                    continue
                a = bootstrap_ci(sc, lab, auroc, nboot, ci, seed + k)
                k += 1
                perf.append({"positives": set_name, "cutoff": c, "eval": wname, "transform": t,
                             "score": s, "n_pos": int(lab[~np.isnan(sc)].sum()),
                             "n_neg": int((~lab[~np.isnan(sc)]).sum()),
                             "auroc": round(a[0], 3), "auroc_lo": round(a[1], 3),
                             "auroc_hi": round(a[2], 3), "avg_precision": round(ap(sc, lab), 4)})
                if s == "full":
                    continue
                d, lo, hi, p = paired_bootstrap_diff(
                    tr[(t, "full")].reindex(ev.index).to_numpy(float), sc, lab, auroc, nboot,
                    ci, seed)
                pairs.append({"positives": set_name, "cutoff": c, "eval": wname,
                              "transform": t, "a": "full", "b": s, "n_pos": int(lab.sum()),
                              "auroc_diff": round(d, 3), "ci_lo": round(lo, 3),
                              "ci_hi": round(hi, 3), "frac_le0": round(p, 3)})
    (root / "interim" / "comparators").mkdir(parents=True, exist_ok=True)
    df[["fold"] + COMPARATORS[2:]].to_csv(
        root / "interim" / "comparators" / f"{set_name}_{c}.tsv.gz", sep="\t", compression="gzip")
    return perf, pairs, cov


jobs = [(sn, c) for sn in TIERS for c in bc["cutoff_years"]]
res = Parallel(n_jobs=n_jobs)(delayed(evaluate)(sn, c, 1000 * i) for i, (sn, c) in enumerate(jobs))
pd.DataFrame([r for p, _, _ in res for r in p]).to_csv(out / "comparator_performance.tsv",
                                                       sep="\t", index=False)
pd.DataFrame([r for _, q, _ in res for r in q]).to_csv(out / "comparator_paired_differences.tsv",
                                                       sep="\t", index=False)
pd.DataFrame([r for _, _, v in res for r in v]).to_csv(out / "comparator_coverage.tsv",
                                                       sep="\t", index=False)


# --- seed-specificity null (tier 1, all genes after the cutoff) ----------------------
# The same comparators seeded with random non-candidate genes (as many as the real
# seeds), raw AUROC. If the real seeds do no better, the comparator's AUROC reflects
# study bias (annotation depth, network degree), not similarity to known genes.
N_NULL = 20


def null_at(c: int, i: int):
    rng = np.random.default_rng(seed + 7919 * (i + 1))
    m, adj = data[c]
    y = training_labels(universe, truth, c, TIERS["tier1"])
    ev = evaluation_labels(truth, negatives, c, TIERS["tier1"])
    ev = ev[ev.index.isin(y.index[y == 0])]
    pool = sorted(negatives & set(y.index[y == 0]))
    fns = {"go_similarity": lambda k: go_similarity_scores(m, universe, k)}
    if adj is not None:
        fns["ppi_rwr"] = lambda k: rwr_scores(adj, universe, k, meth["rwr_restart"])
    real = {n: auroc(seed_scores_at_cutoff(f, y).reindex(ev.index).to_numpy(float),
                     ev.to_numpy().astype(bool)) for n, f in fns.items()}
    rows = []
    for d in range(N_NULL):
        seeds = set(rng.choice(pool, int(y.sum()), replace=False))
        e2 = ev[~ev.index.isin(seeds)]
        lab = e2.to_numpy().astype(bool)
        for name, fn in fns.items():
            rows.append({"cutoff": c, "score": name, "draw": d, "real_seeds_auroc": real[name],
                         "auroc": auroc(fn(seeds).reindex(e2.index).to_numpy(float), lab)})
    return rows


null = pd.DataFrame([r for rr in Parallel(n_jobs=n_jobs)(
    delayed(null_at)(c, i) for i, c in enumerate(bc["cutoff_years"])) for r in rr])
null.groupby(["cutoff", "score", "real_seeds_auroc"]).auroc.agg(
    null_mean="mean", null_lo=lambda v: v.quantile(0.025), null_hi=lambda v: v.quantile(0.975),
    n_draws="size").reset_index().round(3).to_csv(out / "comparator_random_seed_null.tsv",
                                                  sep="\t", index=False)
print("done", flush=True)
